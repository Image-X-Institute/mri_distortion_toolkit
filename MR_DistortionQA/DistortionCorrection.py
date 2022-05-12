import sys, os
import logging
import pydicom
import matplotlib.pyplot as plt
import multiprocessing as mp
from finufft import Plan
import numpy as np
from scipy.fft import fft2
from scipy.fft import fftshift
from scipy.sparse.linalg import lsqr
from scipy.sparse.linalg import LinearOperator
import pandas as pd
from pathlib import Path
from .utilities import get_all_files, convert_cartesian_to_spherical, generate_legendre_basis, dicom_to_numpy
from .utilities import  get_spherical_harmonics
import seaborn as sns

logging.basicConfig(format='[%(filename)s: line %(lineno)d] %(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)

sns.set_theme(style="whitegrid")

class geo:
    """
    extract relevant geometric parameters from input ImageHeader.
    Image header is a pydicom object
    """

    def __init__(self, ImageHeader):
        self.Nrows = int(ImageHeader.Rows)
        self.Ncolumns = int(ImageHeader.Columns)
        self.PixelSizeColumns = float(ImageHeader.PixelSpacing[0])
        self.PixelSizeRows = float(ImageHeader.PixelSpacing[1])
        self.Orientation = np.int_(ImageHeader.ImageOrientationPatient)
        self.Position = np.float_(ImageHeader.ImagePositionPatient)

        if (np.round(self.Orientation) == [0, 1, 0, 0, 0, -1]).all:  # Originally axial images
            # %Axial images isocenter position
            self.IsoRows = self.Nrows - self.Position[2] / self.PixelSizeRows  # isocenter calculated from dicom
            self.IsoColumns = self.Ncolumns + self.Position[1] / self.PixelSizeColumns
            self.Xpixels = 1
            self.Ypixels = self.Nrows
            self.Zpixels = self.Ncolumns
            # BW: this one defintely works.
        elif (np.round(self.Orientation) == [1, 0, 0, 0, 0, -1]).all:  # Originally coronal images
            logger.warning('please check the that the number of pixels in each direction is being calculated correctly.'
                           '\nRemove this warning when checked')
            # %Coronal images isocenter position
            self.IsoRows = self.Nrows - self.Position[2] / self.PixelSizeRows  # isocenter calculated from dicom
            self.IsoColumns = self.Ncolumns + self.Position[0] / self.PixelSizeColumns
            self.Xpixels = self.Nrows
            self.Ypixels = 1
            self.Zpixels = self.Ncolumns
        elif (np.round(self.Orientation) == [1, 0, 0, 0, 1, 0]).all:  # % Originally sagittal images
            # %Sagittal images isocenter position
            logger.warning('please check the that the number of pixels in each direction is being calculated correctly.'
                           '\nRemove this warning when checked')
            self.IsoRows = self.Nrows + self.Position[1] / self.PixelSizeRows  # isocenter calculated from dicom
            self.IsoColumns = self.Ncolumns + self.Position[0] / self.PixelSizeColumns
            self.Xpixels = self.Nrows
            self.Ypixels = self.Ncolumns
            self.Zpixels = 1


class KspaceDistortionCorrector:
    """
    This algorithm is based on `this work`_<https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25487>

    """

    def __init__(self, ImageDirectory, NufftLibrary='finufft', Gx_Harmonics=None, Gy_Harmonics=None,
                 Gz_Harmonics=None, ImExtension='.dcm', r_DSV=None):
        """

        """
        self.Gx_Harmonics, self.Gy_Harmonics, self.Gz_Harmonics = \
            get_spherical_harmonics(Gx_Harmonics, Gy_Harmonics, Gz_Harmonics)
        self.ImageDirectory = Path(ImageDirectory)
        self.n_order = int(np.sqrt(self.Gx_Harmonics.size) -1)
        self.Images = get_all_files(self.ImageDirectory, ImExtension)
        self.r_DSV = r_DSV

        # Select nufft algorithm to use:
        if not NufftLibrary in ['pynufft', 'finufft', 'torchnufft']:
            logger.error(f'Unknown NUFFTlibrary entered: {NufftLibrary}.'
                         f'\nAllowed options are: pynufft, finufft, torchnufft')
            sys.exit(1)
        self.NUFFTlibrary = NufftLibrary  # pynufft, finufft, or torchnufft
        self._check_input_data()
        '''
        I got sick of switching code branches so I decided to just put this all in one. At some point we can revert
        back to just using one.
        '''

    def _check_input_data(self):
        """
        Put any tests of input data here
        """
        assert self.Gx_Harmonics.shape == self.Gy_Harmonics.shape == self.Gz_Harmonics.shape

    def _correct_image(self, image):
        """
        List of steps required for each independant input slice.
        Under some circumstances some steps could be recycled for subsequent slices
        """
        self.CurrentImageName = image
        self.CurrentSlice = pydicom.read_file(self.ImageDirectory / image)
        self.geo = geo(self.CurrentSlice)  # want to get rid of this dependnency
        (X, Y, Z) = dicom_to_numpy(self.ImageDirectory, FilesToReadIn=image, file_extension='dcm', return_XYZ=True)[2]
        '''
        the ordering matters downstream.
        could be meant to be Z Y X - we need to test with images in different slice orientations.
        '''
        self._row_label = np.array(['x [mm]','y [mm]','z [mm]'])[np.equal(X.shape, self.CurrentSlice.Rows)][0]
        self._col_label = np.array(['x [mm]','y [mm]','z [mm]'])[np.equal(X.shape, self.CurrentSlice.Columns)][0]
        self._extent = [Z.min(), Z.max(), Y.min(), Y.max()]  # same for all plots
        self._extent = [Y.min(), Y.max(), Z.min(), Z.max()]  # same for all plots
        self._image_shape = np.array(np.squeeze(X).shape)
        coords = np.array([X.flatten(), Y.flatten(), Z.flatten()])
        self.coords = pd.DataFrame(coords.T, columns=['x', 'y', 'z'])
        self.coords = convert_cartesian_to_spherical(self.coords)
        self._calculate_gradient_strength()
        self._calculate_encoding_fields()
        # self._plot_encoding_fields()  # useful for debugging
        self.GenerateKspaceData()
        self.GenerateDistortedKspace()
        self.PerformLeastSquaresOptimisation()
        self.SaveCorrectedImage()
        self.SaveCorrectedDicom()
    
    def _calculate_gradient_strength(self):
        """
        Calculate the gradient strengths from dicom header.
        Probably have to use these to scale the harmonics... or is it only the relative values that matter...
        """
        bandwidth = float(self.CurrentSlice.PixelBandwidth)
        image_size = np.array([self.coords.x.unique().shape[0], self.coords.y.unique().shape[0], self.coords.z.unique().shape[0]])
        gama = self.CurrentSlice.MagneticFieldStrength * 42.57747851

        # get pixel spacing
        x_pixel_spacing = y_pixel_spacing = z_pixel_spacing = np.nan # overwrite below if possible
        if np.unique(self.coords.x).shape[0] > 1:
            x_pixel_spacing = np.mean(np.diff(np.unique(self.coords.x)))
        if np.unique(self.coords.y).shape[0] > 1:
            y_pixel_spacing = np.mean(np.diff(np.unique(self.coords.y)))
        if np.unique(self.coords.z).shape[0] > 1:
            z_pixel_spacing = np.mean(np.diff(np.unique(self.coords.z)))
        x_fov = self.coords.x.max() - self.coords.x.min() + x_pixel_spacing
        y_fov = self.coords.y.max() - self.coords.y.min() + y_pixel_spacing
        z_fov = self.coords.z.max() - self.coords.z.min() + z_pixel_spacing
        FOV = np.array([x_fov, y_fov, z_fov])
        self._gradient_strength = bandwidth * image_size / (gama*1e6 * FOV*1e-3) # in T/m
        self._scale_factor = 3.33
        logger.warning('scaling gradient strength by hard coded fudge factor')

    def im2double(self, im):
        """
        used to normalise image, replicates matlab function of the same name
        """
        info = np.iinfo(im.dtype)  # Get the data type of the input image
        return im.astype(float) / info.max  # Divide all values by the largest possible

    def ShowImages(self):
        """
        Plot the input and corrected images. If the correction isn't calculated yet will just show the input image
        :return:
        """
        fig, axs = plt.subplots(nrows=1, ncols=2)
        axs[0].imshow(self.CurrentSlice.pixel_array, extent=self._extent)
        axs[0].set_title('Original Image')
        axs[0].set_aspect('equal')
        axs[0].set_xlabel(self._row_label)
        axs[0].set_ylabel(self._col_label)
        axs[0].grid(False)
        if self.r_DSV:
            draw_circle = plt.Circle((0, 0), self.r_DSV, fill=False)
            axs[0].add_artist(draw_circle)

        axs[1].imshow(self.outputImage, extent=self._extent)
        axs[1].set_title('Corrected Image')
        axs[1].set_aspect('equal')
        axs[1].set_xlabel(self._row_label)
        axs[1].set_ylabel(self._col_label)
        axs[1].grid(False)
        if self.r_DSV:
            draw_circle = plt.Circle((0, 0), self.r_DSV, fill=False)
            axs[1].add_artist(draw_circle)

        plt.tight_layout()
        plt.show()

    def _calculate_encoding_fields(self):
        """
        Based on the spherical harmonics and the coordinates derived from the dicom, estimate the encoding
        fields that have been applied to each voxel
        """
        # X_middle = np.zeros([size_m, self.num_order])
        S_v = generate_legendre_basis(self.coords, self.n_order)

        self.Gx_encode = (S_v @ self.Gx_Harmonics) * self._scale_factor
        self.Gy_encode = (S_v @ self.Gy_Harmonics) * self._scale_factor
        self.Gz_encode = (S_v @ self.Gz_Harmonics) * self._scale_factor
        # nb [:,None] is a sneaky way of adding a singleton dimension

    def _plot_encoding_fields(self, vmin=None, vmax=None):
        """
        debug code; make sure the gradient fields look correct
        """

        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=[15, 5])
        gah = np.squeeze(np.reshape(self.Gx_encode.to_numpy(), self._image_shape)) * 1e6
        Xim = axs[0].imshow(gah, extent=self._extent, vmin=vmin, vmax=vmax)
        axs[0].set_title('Gx')
        axs[0].set_xlabel(self._row_label)
        axs[0].set_ylabel(self._col_label)
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[0].add_artist(draw_circle)
        axs[0].grid(False)

        gah = np.squeeze(np.reshape(self.Gy_encode.to_numpy(), self._image_shape)) * 1e6
        Yim = axs[1].imshow(gah, extent=self._extent, vmin=vmin, vmax=vmax)
        axs[1].set_title('Gy')
        axs[1].set_xlabel(self._row_label)
        axs[1].set_ylabel(self._col_label)
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[1].add_artist(draw_circle)
        axs[1].grid(False)

        gah = np.squeeze(np.reshape(self.Gz_encode.to_numpy(), self._image_shape)) * 1e6
        Zim = axs[2].imshow(gah, extent=self._extent, vmin=vmin, vmax=vmax)
        axs[2].set_title('Gz')
        axs[2].set_xlabel(self._row_label)
        axs[2].set_ylabel(self._col_label)
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[2].add_artist(draw_circle)
        axs[2].grid(False)

        cbar = fig.colorbar(Xim, ax=axs[0], fraction=0.046, pad=0.04)
        cbar = fig.colorbar(Yim, ax=axs[1], fraction=0.046, pad=0.04)
        cbar = fig.colorbar(Zim, ax=axs[2], fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.show()

        Bz_compare = Path('../../\MRILinac_DistortionCorrection/gah_zim.csv').resolve()
        gah_shanshan = np.loadtxt(Bz_compare, delimiter=',')
        new_gah = abs(gah-gah_shanshan)
        plt.figure()
        Zim = plt.imshow(new_gah, extent=self._extent, vmin=vmin, vmax=vmax)
        plt.colorbar(Zim, fraction=0.046, pad=0.04)
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs[2].add_artist(draw_circle)
        plt.title('subtracted data')
        plt.show()

        print('hello')

    def GenerateKspaceData(self):
        """
        Get the k space of the current slice data.

        :return:
        """
        self.k_space = fftshift(fft2(fftshift(self.CurrentSlice.pixel_array)))

    def GenerateDistortedKspace(self):
        """
        I don't entirely know what is happening here...
        we now reshape our encoding signals xn_dis into xj and yj plus more scaling factors...
         k_xy_dis represent the k-space indices of non-uniform points.
         ToDo:: is it correct to always have 0.5 as k-space indices? this is what is currently coded...
        :return:
        """

        # reshape to 2D and scale....:
        if (np.round(self.CurrentSlice.ImageOrientationPatient) == [0, 1, 0, 0, 0, -1]).all():
            self.xj = np.array((((self.Gz_encode * 1e5 / self.geo.PixelSizeRows) - (self.geo.IsoRows - self.geo.Nrows / 2) - 1) * 2 * np.pi))
            self.yj = np.array(((self.Gy_encode * 1e5 / self.geo.PixelSizeColumns) - (self.geo.IsoColumns - self.geo.Ncolumns / 2) - 1) * 2 * np.pi)
            # ToDo: this is the last line using geo.
        elif (np.round(self.CurrentSlice.ImageOrientationPatient) == [1, 0, 0, 0, 1, 0]).all():
            self.xj = np.array((((self.Gz_encode * 1e5 / self.geo.PixelSizeRows) - (self.geo.IsoRows - self.geo.Nrows / 2) - 1) * 2 * np.pi))
            self.yj = np.array((((self.Gx_encode * 1e5 / self.geo.PixelSizeColumns) - (self.geo.IsoColumns - self.geo.Ncolumns / 2) - 1) * 2 * np.pi))
        else:
            logger.error('not coded yet')
            sys.exit(1)

        self.nj = self.CurrentSlice.Rows * self.CurrentSlice.Columns
        self.nk = self.CurrentSlice.Rows * self.CurrentSlice.Columns
        # as currently codded it will always come out as 1/2..??
        temp1 = np.linspace(-1 / 2, 1 / 2 - (1 / self.CurrentSlice.Rows), self.CurrentSlice.Rows)
        temp2 = np.linspace(-1 / 2, 1 / 2 - (1 / self.CurrentSlice.Columns), self.CurrentSlice.Columns)
        [T1, T2] = np.meshgrid(temp2, temp1)
        self.sk = T2.flatten()
        self.tk = T1.flatten()

    def fiNufft_Ax(self, x):
        """
        flatron instiute nufft
        Returns A*x
        equivalent to the 'notranpose' option in shanshans code
        self.xj and yj are non uniform nonuniform source points. they are essentially the encoding signals.
        self.sk and tk are uniform target points
        # """
        y = self.Nufft_Ax_Plan.execute(x, None)
        return y.flatten()

    def fiNufft_Atb(self, x):
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

    def PerformLeastSquaresOptimisation(self):
        """
        From
        `Integrated Image Reconstruction and Gradient Nonlinearity Correction <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4390402/pdf/nihms629785.pdf>`_

        "A common strategy for reconstructing both fully- and undersampled MRI data is  penalized regression,
        which seeks to **produce the image** most likely to have produced the set of noisy measurements while potentially
        also satisfying some other expected properties (e.g., sparsity) (17). Since noise in MRI is Gaussian
        distributed, penalized least squares regression of the following general form is often employed.

        So, the image is what we want to get back.
        The noisy measurements are self.k_space
        We use the NUFFT to compute the image most likely to have produced k_space, given the encoding fields that we
        computed.
        :return:
        """
        fk1 = np.reshape(self.k_space, [self.CurrentSlice.Rows * self.CurrentSlice.Columns])
        StartingImage = None  # x0 for lsqr. can be overwritten for each option below.


        if self.NUFFTlibrary == 'finufft':
            self.Nufft_Ax_Plan = Plan(3, 2, 1, 1e-06, -1)
            self.Nufft_Ax_Plan.setpts(self.xj, self.yj, None, self.sk, self.tk)
            self.Nufft_Atb_Plan = Plan(3, 2, 1, 1e-06, 1)
            self.Nufft_Atb_Plan.setpts(self.sk, self.tk, None, self.xj, self.yj)
            A = LinearOperator((fk1.shape[0], fk1.shape[0]), matvec=self.fiNufft_Ax, rmatvec=self.fiNufft_Atb)
            StartingImage = self.CurrentSlice.pixel_array.flatten().astype(complex)

        maxit = 20
        x1 = lsqr(A, fk1, iter_lim=maxit, x0=StartingImage)

        self.outputImage = abs(np.reshape(x1[0], [self.CurrentSlice.Rows, self.CurrentSlice.Columns]))

    def SaveCorrectedImage(self):
        """
        Quick and dirty image save
        :return:
        """

        if not os.path.isdir(self.ImageDirectory / 'Corrected'):
            os.mkdir(self.ImageDirectory / 'Corrected')

        fig, axs = plt.subplots(nrows=1, ncols=2)
        # axs[0].imshow(self.CurrentSlice.pixel_array, extent=self._extent)
        axs[0].imshow(self.CurrentSlice.pixel_array)
        axs[0].set_title('Original Image')
        axs[0].set_aspect('equal')
        axs[0].set_xlabel(self._row_label)
        axs[0].set_ylabel(self._col_label)
        axs[0].grid(False)
        if self.r_DSV:
            draw_circle = plt.Circle((0, 0), self.r_DSV, fill=False)
            axs[0].add_artist(draw_circle)

        # axs[1].imshow(self.outputImage, extent=self._extent)
        axs[1].imshow(self.outputImage)
        axs[1].set_title('Corrected Image')
        axs[1].set_aspect('equal')
        axs[1].set_xlabel(self._row_label)
        axs[1].set_ylabel(self._col_label)
        axs[1].grid(False)
        if self.r_DSV:
            draw_circle = plt.Circle((0, 0), self.r_DSV, fill=False)
            axs[1].add_artist(draw_circle)

        plt.tight_layout()
        plt.savefig(self.ImageDirectory / 'Corrected' / (self.CurrentImageName + '.png'), format='png', dpi=600)
        plt.close(fig)

    def SaveCorrectedDicom(self):


        if not os.path.isdir(self.ImageDirectory / 'Corrected_dcm'):
            os.mkdir(self.ImageDirectory / 'Corrected_dcm')



        try:
            self._output_dcm_num =  self._output_dcm_num + 1
        except AttributeError:
            self._output_dcm_num = 0
        dcm_name = str(self._output_dcm_num) + '.dcm'
        if dcm_name == '25.dcm':
            print('hello')
            plt.figure()
            plt.imshow(self.CurrentSlice.pixel_array)

        temp_dcm = self.CurrentSlice.copy()
        temp_dcm.PixelData = np.uint16(self.outputImage).tobytes()
        temp_dcm.save_as(self.ImageDirectory / 'Corrected_dcm' / dcm_name)




    def CorrectAllImaFiles(self, ParralelProcessing=False):
        """
        This will loop through and correct all IMA files in the input directory
        :return:
        """
        if ParralelProcessing:
            try:
                processes = []
                for image in self.Images:
                    p = mp.Process(target=self._correct_image, args=(image,))
                    processes.append(p)
                [x.start() for x in processes]
            except:
                logger.error('Multi processing failed with the following error:')
                raise
        else:
            for image in self.Images:
                self._correct_image(image)
