import sys, os
import logging
import warnings
import pydicom
import matplotlib.pyplot as plt
from finufft import Plan
import numpy as np
from scipy.fft import fft2
from scipy.fft import fftshift
from scipy.sparse.linalg import lsqr
from scipy.sparse.linalg import LinearOperator
import pandas as pd
from pathlib import Path
from .utilities import get_all_files, convert_cartesian_to_spherical, generate_legendre_basis, dicom_to_numpy
from .utilities import get_gradient_spherical_harmonics
from .utilities import printProgressBar
from time import perf_counter



logging.basicConfig(format='[%(filename)s: line %(lineno)d] %(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)


class KspaceDistortionCorrector:
    """
    This algorithm is based on `this work`_<https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25487>
    """

    def __init__(self, ImageDirectory, NufftLibrary='finufft', Gx_Harmonics=None, Gy_Harmonics=None,
                 Gz_Harmonics=None, ImExtension='.dcm', dicom_data=None, correct_through_plane=True):
        """
        :param ImageDirectory:
        :param NufftLibrary:
        :param Gx_Harmonics:
        :param Gy_Harmonics:
        :param Gz_Harmonics:
        :param ImExtension:
        :param dicom_data:
        :param correct_through_plane:
        """
        self.correct_through_plane = correct_through_plane
        self._n_zero_pad = 0  # n_pixels to add around each edge of volume. set to 0 for no zero padding
        self._dicom_data = dicom_data
        self._calculate_gradient_strength()
        self._Gx_Harmonics, self._Gy_Harmonics, self._Gz_Harmonics = \
            get_gradient_spherical_harmonics(Gx_Harmonics, Gy_Harmonics, Gz_Harmonics)
        self._Gx_Harmonics = self._Gx_Harmonics * self._gradient_strength[0]
        self._Gy_Harmonics = self._Gy_Harmonics * self._gradient_strength[1]
        self._Gz_Harmonics = self._Gz_Harmonics * self._gradient_strength[2] * -1

        self.ImageDirectory = Path(ImageDirectory)
        self._all_dicom_files = get_all_files(self.ImageDirectory, ImExtension)
        self._all_dicom_files = np.sort(self._all_dicom_files)
        warnings.warn('making assumption that sorting by dicom name is same as series number')
        # self._all_dicom_files = sort_dicom_slices(self._all_dicom_files)
        self._n_dicom_files = len(self._all_dicom_files)


        self.ImageArray,\
        self._dicom_affine,\
        (self._X, self._Y, self._Z) = dicom_to_numpy(self.ImageDirectory,
                                                     file_extension='dcm',
                                                     return_XYZ=True,
                                                     zero_pad=self._n_zero_pad,
                                                     enforce_increasing_coords=False)
        self._get_rows_and_cols()
        self.n_order = int(np.sqrt(self._Gx_Harmonics.size) - 1)
        self.Images = get_all_files(self.ImageDirectory, ImExtension)
        self.r_DSV = 150  # only for drawing on the plot
        # Select nufft algorithm to use:
        if not NufftLibrary in ['pynufft', 'finufft', 'torchnufft']:
            raise NotImplementedError(f'Unknown NUFFTlibrary entered: {NufftLibrary}.'
                         f'\nAllowed options are: pynufft, finufft, torchnufft')
        self.NUFFTlibrary = NufftLibrary  # pynufft, finufft, or torchnufft
        self._check_input_data()

    def _get_rows_and_cols(self):
        """
        just extract the number of rows and columns from a dicom header - we need it for zero padding later.
        I suspect I should be able to derive this from the dicom affine, but this way seems more foolproof
        """
        demo_header = pydicom.read_file(self.ImageDirectory / self._all_dicom_files[0])
        self._Rows, self._Cols = self.ImageArray.shape[0:2]

        self._ImageOrientationPatient = demo_header.ImageOrientationPatient

        self._PixelSpacing = self._dicom_data['pixel_spacing']


    def _check_input_data(self):
        """
        Put any tests of input data here
        """
        assert self._Gx_Harmonics.shape == self._Gy_Harmonics.shape == self._Gz_Harmonics.shape
        assert self._n_zero_pad >= 0

    def _correct_image(self):
        """
        List of steps required for each independant input slice.
        Under some circumstances some steps could be recycled for subsequent slices
        """

        self._image_shape = np.array(np.squeeze(self._X_slice).shape)
        coords_cartesian = np.array([self._X_slice.flatten(), self._Y_slice.flatten(), self._Z_slice.flatten()])
        coords_cartesian = pd.DataFrame(coords_cartesian.T, columns=['x', 'y', 'z'])
        self.coords = convert_cartesian_to_spherical(coords_cartesian)
        self._calculate_encoding_fields()
        # self._plot_encoding_fields()  # useful for debugging
        self._generate_Kspace_data()
        self._generate_distorted_indices()
        self._perform_least_squares_optimisation()

    def _unpad_image_arrays(self):
        """
        remove any zero padding from original and reconstructed image volumes
        """
        if self._n_zero_pad > 0:
            self.ImageArray = self.ImageArray[self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad]
            self._image_array_corrected = self._image_array_corrected[self._n_zero_pad:-self._n_zero_pad,
                                          self._n_zero_pad:-self._n_zero_pad,
                                          self._n_zero_pad:-self._n_zero_pad]
            # self._image_array_corrected2 = self._image_array_corrected2[self._n_zero_pad:-self._n_zero_pad,
            #                               self._n_zero_pad:-self._n_zero_pad,
            #                               self._n_zero_pad:-self._n_zero_pad]

    def _calculate_gradient_strength(self):
        """
        Calculate the gradient strengths from dicom header.
        Probably have to use these to scale the harmonics... or is it only the relative values that matter...
        """
        self._gradient_strength = [1e3, 1e3, 1e3]

    def _calculate_encoding_fields(self):
        """
        Based on the spherical harmonics and the coordinates derived from the dicom, estimate the encoding
        fields that have been applied to each voxel
        """
        # X_middle = np.zeros([size_m, self.num_order])

        legendre_basis = generate_legendre_basis(self.coords, self.n_order)
        self.Gx_encode = (legendre_basis @ self._Gx_Harmonics)
        self.Gy_encode = (legendre_basis @ self._Gy_Harmonics)
        self.Gz_encode = (legendre_basis @ self._Gz_Harmonics)

    def _plot_encoding_fields(self, vmin=None, vmax=None):
        """
        debug code; make sure the gradient fields look correct
        """
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=[15, 5])
        gah = np.squeeze(np.reshape(self.Gx_encode.to_numpy(), self._image_shape)) * 1e6

        Xim = axs[0].imshow(gah, vmin=vmin, vmax=vmax)
        axs[0].set_title('Gx')
        # axs[0].set_xlabel(self._row_label)
        # axs[0].set_ylabel(self._col_label)

        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[0].add_artist(draw_circle)
        axs[0].grid(False)

        gah = np.squeeze(np.reshape(self.Gy_encode.to_numpy(), self._image_shape)) * 1e6

        Yim = axs[1].imshow(gah, vmin=vmin, vmax=vmax)
        axs[1].set_title('Gy')
        # axs[1].set_xlabel(self._row_label)
        # axs[1].set_ylabel(self._col_label)

        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[1].add_artist(draw_circle)
        axs[1].grid(False)

        gah = np.squeeze(np.reshape(self.Gz_encode.to_numpy(), self._image_shape)) * 1e6

        Zim = axs[2].imshow(gah, vmin=vmin, vmax=vmax)
        axs[2].set_title('Gz')
        # axs[2].set_xlabel(self._row_label)
        # axs[2].set_ylabel(self._col_label)

        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[2].add_artist(draw_circle)
        axs[2].grid(False)

        cbar = fig.colorbar(Xim, ax=axs[0], fraction=0.046, pad=0.04)
        cbar = fig.colorbar(Yim, ax=axs[1], fraction=0.046, pad=0.04)
        cbar = fig.colorbar(Zim, ax=axs[2], fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.show()

    def _plot_indices(self):
        """
        handy debug plot routine to plot linear/distorted indices
        :return:
        """

        fig, axs = plt.subplots(nrows=2,ncols=2, figsize=[8,8])
        axs[0, 0].plot(self.xj); axs[0, 0].set_title('xj')
        axs[0, 1].plot(self.yj); axs[0, 1].set_title('yj')
        axs[1, 0].plot(self.sk); axs[1, 0].set_title('sk')
        axs[1, 1].plot(self.tk); axs[1, 1].set_title('tk')

    def _generate_Kspace_data(self):
        """
        Get the k space of the current slice data.
        :return:
        """
        self.k_space = fftshift(fft2(fftshift(self._image_to_correct)))

    def _generate_distorted_indices(self):
        """
        generates both the linear (i.e. assumed) indices and the distorted indices.
        These indices are passed to the NUFFT library.
        """

        self.nj = self._Rows * self._Cols
        self.nk = self._Rows * self._Cols
        # as currently codded it will always come out as 1/2..??
        temp1 = np.linspace(-1/2, 1/2, self._Rows)
        temp2 = np.linspace(-1/2, 1/2, self._Cols)
        [T1, T2] = np.meshgrid(temp2, temp1)
        self.sk = T2.flatten()
        self.tk = T1.flatten()

        if (np.round(self._ImageOrientationPatient) == [0, 1, 0, 0, 0, -1]).all():
            xn_dis = self.Gz_encode / (self._PixelSpacing[2])
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gy_encode / (self._PixelSpacing[1])
            self.yj = yn_dis * 2 * np.pi
        elif (np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 0, -1]).all():
            xn_dis = self.Gz_encode / (self._PixelSpacing[2])
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gx_encode / (self._PixelSpacing[0])
            self.yj = yn_dis * 2 * np.pi
        elif (np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 1, 0]).all():
            xn_dis = self.Gy_encode / (self._PixelSpacing[1])
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gx_encode / (self._PixelSpacing[0])
            self.yj = yn_dis * 2 * np.pi
        elif np.round(self._ImageOrientationPatient == [2, 2, 2, 2, 2, 2]).all():
            # this is for through plane correction where the real images are [0, 1, 0, 0, 0, -1]
            x_lin_size, y_lin_size = self._image_to_correct.shape
            xn_lin = np.linspace(-x_lin_size/2, x_lin_size/2, x_lin_size)
            yn_lin = np.linspace(-y_lin_size/2, y_lin_size/2, y_lin_size)
            [xn_lin, yn_lin] = np.meshgrid(xn_lin, yn_lin, indexing='ij')
            xn_lin = xn_lin.flatten()
            yn_lin = yn_lin.flatten()
            xn_dis = self.Gx_encode / (self._PixelSpacing[0])
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            yn_dis = self.Gx_encode / (self._PixelSpacing[0])
            self.yj = yn_dis * 2 * np.pi
        elif np.round(self._ImageOrientationPatient == [3, 3, 3, 3, 3, 3]).all():
            # this is for through plane correction where the real images are [1, 0, 0, 0, 0, -1]
            x_lin_size, y_lin_size = self._image_to_correct.shape
            xn_lin = np.linspace(-x_lin_size / 2, x_lin_size / 2, x_lin_size)
            yn_lin = np.linspace(-y_lin_size / 2, y_lin_size / 2, y_lin_size)
            [xn_lin, yn_lin] = np.meshgrid(xn_lin, yn_lin, indexing='ij')
            xn_lin = xn_lin.flatten()
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            yn_dis = self.Gy_encode / (self._PixelSpacing[1])
            self.yj = yn_dis * 2 * np.pi
        elif (np.round(self._ImageOrientationPatient) == [1, 1, 1, 1, 1, 1]).all():
            # this is for through plane correction where the real images are [1, 0, 0, 0, 1, 0]
            x_lin_size, y_lin_size = self._image_to_correct.shape
            xn_lin = np.linspace(-x_lin_size/2, x_lin_size/2, x_lin_size)
            yn_lin = np.linspace(-y_lin_size/2, y_lin_size/2, y_lin_size)
            [xn_lin, yn_lin] = np.meshgrid(xn_lin, yn_lin, indexing='ij')
            xn_lin = xn_lin.flatten()
            #xn_dis should match the image indices
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            # self.xj = pd.Series(self.sk)
            yn_dis = -1*self.Gz_encode / (self._PixelSpacing[2])
            self.yj = yn_dis * 2 * np.pi


        else:
            raise NotImplementedError('this slice orientation is not handled yet sorry')
        try:
            self.xj = self.xj.to_numpy()
            self.yj = self.yj.to_numpy()
        except:
            print('fucks ache')


        if (np.round(self._ImageOrientationPatient) == [1, 1, 1, 1, 1, 1]).all() or \
                (np.round(self._ImageOrientationPatient) == [2,2,2,2,2,2]).all():
            try:
                self.dodgy_ind = self.dodgy_ind + 1
            except:
                self.dodgy_ind = 0

    def _fiNufft_Ax(self, x):
        """
        flatron instiute nufft
        Returns A*x
        equivalent to the 'notranpose' option in shanshans code
        self.xj and yj are non uniform nonuniform source points. they are essentially the encoding signals.
        self.sk and tk are uniform target points
        # """

        if x.dtype is not np.dtype('complex128'):
            x = x.astype('complex128')
        # y = nufft2d3(self.xj, self.yj, x, self.sk, self.tk, eps=1e-06, isign=-1)
        y = self.Nufft_Ax_Plan.execute(x, None)
        return y.flatten()

    def _fiNufft_Atb(self, x):
        """
        flatron instiute nufft
        This is to define the Nufft as a scipy.sparse.linalg.LinearOperator which can be used by the lsqr algorithm
        see here for explanation:
        https://stackoverflow.com/questions/48621407/python-equivalent-of-matlabs-lsqr-with-first-argument-a-function
        Returns A'*x
        equivalent to the 'tranpose' option in shanshans code
        """
        # y = nufft2d3(self.sk, self.tk, x, self.xj, self.yj, eps=1e-06, isign=1)
        y = self.Nufft_Atb_Plan.execute(x, None)
        return y.flatten()

    def _perform_least_squares_optimisation(self):
        """
        From
        `Integrated Image Reconstruction and Gradient Nonlinearity Correction <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4390402/pdf/nihms629785.pdf>`_
        "A common strategy for reconstructing both fully- and undersampled MRI data is  penalized regression,
        which seeks to **produce the image** most likely to have produced the set of noisy measurements while potentially
        also satisfying some other expected properties (e.g., sparsity) (17). Since noise in MRI is Gaussian
        distributed, penalized least squares regression of the following general form is often employed.
        So:
        - the image is what we want to get back.
        - The noisy measurements are self.k_space
        - We use the NUFFT to compute the image most likely to have produced k_space, given the encoding fields that we computed.
        """


        fk1 = np.reshape(self.k_space, [self._Rows * self._Cols])

        StartingImage = None  # x0 for lsqr. can be overwritten for each option below.

        if self.NUFFTlibrary == 'finufft':
            self.Nufft_Ax_Plan = Plan(3, 2, 1, 1e-06, -1)
            self.Nufft_Ax_Plan.setpts(self.xj, self.yj, None, self.sk, self.tk)
            self.Nufft_Atb_Plan = Plan(3, 2, 1, 1e-06, 1)
            self.Nufft_Atb_Plan.setpts(self.sk, self.tk, None, self.xj, self.yj)
            A = LinearOperator((fk1.shape[0], fk1.shape[0]), matvec=self._fiNufft_Ax, rmatvec=self._fiNufft_Atb)
            StartingImage = self._image_to_correct.flatten().astype(complex)

        if False:
            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[10,5])
            axs[0].imshow(self._image_to_correct)
            axs[1].imshow(self.outputImage)

        maxit = 20
        x1 = lsqr(A, fk1, iter_lim=maxit, x0=StartingImage)
        self.outputImage = abs(np.reshape(x1[0], [self._Rows, self._Cols]))

    def _force_linear_harmonics(self, dims_to_force):

        if 'x' in dims_to_force:
            _Gx_Harmonics_tmp = np.zeros(self._Gx_Harmonics.shape)
            _Gx_Harmonics_tmp[2] = self._Gx_Harmonics[2]
            self._Gx_Harmonics = _Gx_Harmonics_tmp

        if 'y' in dims_to_force:
            _Gy_Harmonics_tmp = np.zeros(self._Gy_Harmonics.shape)
            _Gy_Harmonics_tmp[3] = self._Gy_Harmonics[3]
            self._Gy_Harmonics = _Gy_Harmonics_tmp

        if 'z' in dims_to_force:
            _Gz_Harmonics_tmp = np.zeros(self._Gz_Harmonics.shape)
            _Gz_Harmonics_tmp[1] = self._Gz_Harmonics[1]
            self._Gz_Harmonics = _Gz_Harmonics_tmp

    # public methods

    def correct_all_images(self):
        """
        This will loop through and correct all IMA files in the input directory
        :return:
        """
        if self.correct_through_plane:
            n_images_to_correct = self._n_dicom_files + (self.ImageArray.shape[1] - 2 * self._n_zero_pad)
        else:
            n_images_to_correct = self._n_dicom_files
        loop_axis = 2  # ImageArray always has slice last
        # nb: the below code allows us to wrap all the data into one 'itterable'; it also changes the view
        # such that we are looping over the slice direction
        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._X, loop_axis),
                          np.rollaxis(self._Y, loop_axis),
                          np.rollaxis(self._Z, loop_axis))

        i = 0
        self._image_array_corrected = np.zeros(self.ImageArray.shape)
        for array_slice, X, Y, Z in zipped_data:
            if i < self._n_zero_pad or i > self._n_dicom_files + self._n_zero_pad:
                # skip these files, they have no meaning anyway and we can't (easily) write them to dicom
                i += 1
                continue

            t_start = perf_counter()
            print(f'2D correction: {i-self._n_zero_pad} of {self._n_dicom_files}')

            print(printProgressBar(i-self._n_zero_pad, n_images_to_correct))
            self._image_to_correct = array_slice
            self._X_slice = X
            self._Y_slice = Y
            self._Z_slice = Z
            self._correct_image()
            self._image_array_corrected[:, :, i] = self.outputImage
            t_stop = perf_counter()
            print(f"Elapsed time {t_stop - t_start}")
            i += 1

        if self.correct_through_plane:

            if (np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 1, 0]).all():
                self._ImageOrientationPatient = [1, 1, 1, 1, 1, 1]
            elif np.round(self._ImageOrientationPatient == [0, 1, 0, 0, 0, -1]).all():
                self._ImageOrientationPatient = [2, 2, 2, 2, 2, 2]
            elif np.round(self._ImageOrientationPatient == [1, 0, 0, 0, 0, -1]).all():
                self._ImageOrientationPatient = [3, 3, 3, 3, 3, 3]
            else:
                raise NotImplementedError
            # which directions are already corrected:
            corrected_dims = np.array(['x', 'y', 'z'])[([self._dicom_data['slice_direction'] not in
                                                         axis for axis in ['x', 'y', 'z']])]

            # corrected_dims = np.array(['x', 'y', 'z'])
            # self._force_linear_harmonics(corrected_dims)  # this forces already corrected dimensions to be linear
            loop_axis = 1
            zipped_data = zip(np.rollaxis(self._image_array_corrected, loop_axis),
                              np.rollaxis(self._X, loop_axis),
                              np.rollaxis(self._Y, loop_axis),
                              np.rollaxis(self._Z, loop_axis))
            self._image_array_corrected2 = self._image_array_corrected.copy()
            j = 0
            self._Rows = np.rollaxis(self._image_array_corrected, loop_axis).shape[1]
            self._Cols = np.rollaxis(self._image_array_corrected, loop_axis).shape[2]
            n_slices = np.rollaxis(self._image_array_corrected, loop_axis).shape[0]
            for array_slice, X, Y, Z in zipped_data:
                if j < self._n_zero_pad or j > n_slices + self._n_zero_pad:
                    # skip these files, they have no meaning anyway and we can't (easily) write them to dicom
                    j += 1
                    continue
                t_start = perf_counter()
                print(f'Through Plane correction: {j - self._n_zero_pad} of {self.ImageArray.shape[loop_axis] - (2*self._n_zero_pad)}')
                print(printProgressBar(i+j-(self._n_zero_pad*4), n_images_to_correct))
                self._image_to_correct = array_slice
                self._X_slice = X
                self._Y_slice = Y
                self._Z_slice = Z
                self._correct_image()
                self._image_array_corrected[:, j, :] = self.outputImage
                if False:
                    fig, axs = plt.subplots(nrows=1, ncols=3)
                    axs[0].imshow(self._image_array_corrected[:, j, :] )
                    axs[1].imshow(self._image_to_correct)
                    axs[2].imshow(self.outputImage)
                t_stop = perf_counter()
                print(f"Elapsed time {t_stop - t_start}")

                j += 1


        self._unpad_image_arrays()

    def save_all_images(self):
        if not os.path.isdir(self.ImageDirectory / 'Corrected'):
            os.mkdir(self.ImageDirectory / 'Corrected')
        loop_axis = np.where([self._dicom_data['slice_direction'] in axis for axis in ['x', 'y', 'z']])[0][0]
        loop_axis = 2

        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._image_array_corrected, loop_axis))
        i = 0
        for original_image, corrected_image in zipped_data:
            fig, axs = plt.subplots(nrows=1, ncols=2)

            # axs[0].imshow(self.pixel_array, extent=self._extent)
            im1 = axs[0].imshow(original_image)
            axs[0].set_title('Original Image')
            axs[0].set_aspect('equal')
            # axs[0].set_xlabel(self._row_label)
            # axs[0].set_ylabel(self._col_label)
            axs[0].grid(True)
            # fig.colorbar(im1, ax=axs[0])

            # axs[1].imshow(self.outputImage, extent=self._extent)
            im2 = axs[1].imshow(corrected_image)
            axs[1].set_title('Corrected Image')
            axs[1].set_aspect('equal')
            # axs[1].set_xlabel(self._row_label)
            # axs[1].set_ylabel(self._col_label)
            axs[1].grid(True)
            # fig.colorbar(im2, ax=axs[1])

            plt.tight_layout()
            plt.savefig(self.ImageDirectory / 'Corrected' / (str(i) + '.png'), format='png')
            plt.close(fig)
            i += 1

    def save_all_images_as_dicom(self):
        if not os.path.isdir(self.ImageDirectory / 'Corrected_dcm'):
            os.mkdir(self.ImageDirectory / 'Corrected_dcm')
        loop_axis = np.where([self._dicom_data['slice_direction'] in axis for axis in ['x', 'y', 'z']])[0][0]
        loop_axis = 2
        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._image_array_corrected, loop_axis))
        i = 0
        for original_image, corrected_image in zipped_data:
            current_slice = pydicom.read_file(self.ImageDirectory / self._all_dicom_files[i])
            temp_dcm = current_slice.copy()
            temp_dcm.PixelData = np.uint16(corrected_image).tobytes()
            temp_dcm.save_as(self.ImageDirectory / 'Corrected_dcm' / (str(i) + '.dcm'))
            i += 1

