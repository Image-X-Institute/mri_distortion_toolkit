import math
import os
import logging
import sys
import warnings
import pydicom
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from finufft import Plan
from finufft import nufft2d1, nufft2d2
import numpy as np
from scipy.fft import fft2
from scipy.fft import fftshift
from scipy.sparse.linalg import lsqr
from scipy.sparse.linalg import LinearOperator
from scipy.interpolate import RectBivariateSpline
import pandas as pd
from pathlib import Path
from .utilities import get_all_files, convert_cartesian_to_spherical, generate_legendre_basis, dicom_to_numpy
from .utilities import sort_dicom_slices
from .utilities import get_harmonics
from .utilities import printProgressBar
from time import perf_counter
from abc import abstractmethod
from .utilities import combine_harmonics

logging.basicConfig(format='[%(filename)s: line %(lineno)d] %(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)


class DistortionCorrectorBase:
    """
    This is the base class which other distortion correction methods can inherit from, and includes the core required
    behavior.

    :param ImageDirectory: Directory containing images to correct. Note that all images matching ImExtension will be
        read in; therefore only one series should be inside ImageDirectory (we dont check this!)
    :type ImageDirectory: pathlib.Path or str
    :param gradient_harmonics: a list of three harmonics, which should either be csv files exported from an instance of
        SphericalHarmonicFit, or SphericalHarmonicFit.harmonics
    :type gradient_harmonics: list
    :param B0_harmonics: Instance of SphericalHarmonicFit.harmonics describing B0, or else an equivalent csv file
    :type B0_harmonics: SphericalHarmonicFit.harmonics or csv, optional
    :param ImExtension: extension of files to read in ImageDirectory, e.g. 'dcm', or 'IMA'
    :type ImExtension: str, optional
    :param correct_through_plane: if True, through plane (3D) correction is carried out, which is roughly twice as slow
        as 2D only
    :type correct_through_plane: bool. optional
    :param correct_B0: if True, and if B0_harmonics is supplied, will also correct for B0 effects. Probably only works
        for standard (non EPI) sequences.
    :type correct_B0: bool, optional
    :param B0_direction: 'forward' or 'back'. Determines whether the B0 harmonics should be added (forward) or
        subtracted (backward) from the gradient harmonics
    :type B0_direction: str, optional
    :param pad: pixels to 0 pad image volume. This can be useful in the k-space domain
    :type pad: int, optional
    """

    def __init__(self, ImageDirectory, gradient_harmonics, B0_harmonics=None,
                 dicom_data=None,
                 ImExtension='.dcm', correct_through_plane=True, correct_B0=False,
                 B0_direction='forward', pad=0):
        """
        init method
        """

        if 'DistortionCorrectorBase' in str(self.__class__):
            raise TypeError('DistortionCorrectorBase should not be called directly; it only exists for other correctors'
                            'to inherit. Quitting')

        self.correct_through_plane = correct_through_plane
        self.correct_B0 = correct_B0
        self.B0_direction = B0_direction
        self._n_zero_pad = pad  # n_pixels to add around each edge of volume. set to 0 for no zero padding
        self._Gx_Harmonics, self._Gy_Harmonics, self._Gz_Harmonics, self._B0_Harmonics = \
            get_harmonics(gradient_harmonics[0], gradient_harmonics[1], gradient_harmonics[2], B0_harmonics)
        self._Gx_Harmonics = self._Gx_Harmonics * 1
        self._Gy_Harmonics = self._Gy_Harmonics * 1
        self._Gz_Harmonics = self._Gz_Harmonics * 1
        self._dicom_data = dicom_data

        if self.correct_B0:
            if self._B0_Harmonics is None:
                warnings.warn('cannot correct B0 as no B0 harmonics supplied, continuing')
            else:
                self._add_B0_to_gradient()

        self.ImageDirectory = Path(ImageDirectory)
        self._all_dicom_files = get_all_files(self.ImageDirectory, ImExtension)
        CompletePathFiles = [self.ImageDirectory / file for file in self._all_dicom_files]
        dicom_slices = [pydicom.read_file(f) for f in CompletePathFiles]
        sort_ind = sort_dicom_slices(dicom_slices)[1]
        self._all_dicom_files = list(np.array(self._all_dicom_files)[sort_ind])
        self._n_dicom_files = len(self._all_dicom_files)

        self.ImageArray, \
        self._dicom_affine, \
        (self._X, self._Y, self._Z) = dicom_to_numpy(self.ImageDirectory,
                                                     file_extension=ImExtension,
                                                     return_XYZ=True,
                                                     zero_pad=self._n_zero_pad,
                                                     enforce_increasing_coords=False)
        self._get_rows_and_cols()
        self._get_iso_offset()
        self.n_order = int(np.sqrt(self._Gx_Harmonics.size) - 1)
        self.Images = get_all_files(self.ImageDirectory, ImExtension)
        self.r_DSV = 150  # only for drawing on the plot
        self._check_input_data()

    @abstractmethod
    def _correct_image(self):
        """
        This is the method to correct a 2D slice of data; each inheriting corrector must supply its own _correct_image
        method.
        """
        pass

    def _get_rows_and_cols(self):
        """
        just extract the number of rows and columns from a dicom header - we need it for zero padding later.
        I suspect I should be able to derive this from the dicom affine, but this way seems more foolproof
        """
        demo_header = pydicom.read_file(self.ImageDirectory / self._all_dicom_files[0])
        self._Rows, self._Cols, self._Slices = self.ImageArray.shape
        self._ImageOrientationPatient = demo_header.ImageOrientationPatient
        self._pixel_spacing = self._dicom_affine[0:3, 0:3].sum(1)
        self._pixel_spacing = [ps if not ps == 0 else np.nan for ps in self._pixel_spacing]

    def _get_iso_offset(self):
        """
        calculate the offset between the center of the volume and scanner isocenter.
        This offset is used in _calculate_encoding_indices
        """
        self._image_position_patient_start = [self._X[0, 0, 0], self._Y[0, 0, 0], self._Z[0, 0, 0]]
        assert np.allclose(self._image_position_patient_start, self._dicom_affine[0:3, 3])
        self._image_position_patient_end = [self._X[-1, -1, -1], self._Y[-1, -1, -1], self._Z[-1, -1, -1]]
        self._iso_offset = -1 * np.add(self._image_position_patient_end, self._image_position_patient_start) / 2
        pixel_spacing = self._dicom_affine[0:3, 0:3].sum(1)
        self._iso_offset_pixels = np.divide(self._iso_offset, self._pixel_spacing)

    def _check_input_data(self):
        """
        Put any tests of input data here
        """
        assert self._Gx_Harmonics.shape == self._Gy_Harmonics.shape == self._Gz_Harmonics.shape
        assert self._n_zero_pad >= 0

    def _add_B0_to_gradient(self):
        """
        combine the B0 harmonics with the gradient harmonics in the frequency encode direciton
        :return:
        """

        assert self.B0_direction == 'forward' or self.B0_direction == 'back'
        if self.B0_direction == 'forward':
            _operation = 'add'
        elif self.B0_direction == 'back':
            _operation = 'subtract'
        else:
            raise TypeError('B0_direction must be "forward" or "back"')

        try:
            freq_encode_direction = self._dicom_data['freq_encode_direction']
            _test = self._dicom_data['gradient_strength']
        except (AttributeError, TypeError):
            raise AttributeError('Cannot combine gradient and B0 harmonics as no dicom_data was supplied. '
                           'At a minimum, we need the following: '
                           'dicom_data = {"freq_encode_direction": e.g. "x",'
                           '              "gradient_strength": e.g. [1,2,3]}'
                           'continuing with just the supplied gradient harmonics')
            
        if freq_encode_direction == 'x':
            gradient_order = int(np.sqrt(self._Gx_Harmonics.shape[0] - 1))
            scale = 1/self._dicom_data['gradient_strength'][0]
            self._Gx_Harmonics = combine_harmonics(self._Gx_Harmonics, self._B0_Harmonics*scale, operation=_operation,
                                                   n_order_return=gradient_order)
        elif freq_encode_direction == 'y':
            gradient_order = int(np.sqrt(self._Gy_Harmonics.shape[0] - 1))
            scale = 1 / self._dicom_data['gradient_strength'][1]
            self._Gy_Harmonics = combine_harmonics(self._Gy_Harmonics, self._B0_Harmonics*scale, operation=_operation,
                                                   n_order_return=gradient_order)
        elif freq_encode_direction == 'z':
            gradient_order = int(np.sqrt(self._Gz_Harmonics.shape[0] - 1))
            scale = 1 / self._dicom_data['gradient_strength'][2]
            self._Gz_Harmonics = combine_harmonics(self._Gz_Harmonics, self._B0_Harmonics*scale, operation=_operation,
                                                   n_order_return=gradient_order)

        # from .Harmonics import SphericalHarmonicFit
        # input_data = pd.DataFrame({'x': [0, 0, 0], 'y': [0, 0, 0], 'z': [0, 0, 0], 'Bz': [0, 1, 2]})
        # tempH = SphericalHarmonicFit(input_data)
        # tempH.harmonics = self._Gx_Harmonics
        # tempH._assess_harmonic_pk_pk()
        # tempH.plot_harmonics_pk_pk(cut_off=.005)

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
            self._X = self._X[self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad]
            self._Y = self._Y[self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad]
            self._Z = self._Z[self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad,
                              self._n_zero_pad:-self._n_zero_pad]

    def _calculate_encoding_fields(self):
        """
        Based on the spherical harmonics and the coordinates derived from the dicom, estimate the encoding
        fields that have been applied to each voxel
        """
        legendre_basis = generate_legendre_basis(self.coords, self.n_order)
        self.Gx_encode = (legendre_basis @ self._Gx_Harmonics) * 1e3
        self.Gy_encode = (legendre_basis @ self._Gy_Harmonics) * 1e3
        self.Gz_encode = (legendre_basis @ self._Gz_Harmonics) * 1e3

    def _zero_volume(self, volume_to_zero):
        """
        takes an input volume and replaces all values <0 with 0.
        required for the image domain corrector which tends to produce -ve pixels
        """
        negative_ind = volume_to_zero < 0
        logger.warning(
            f'{np.count_nonzero(negative_ind)} negative pixels detected; setting these to zero and continuing')
        volume_to_zero[negative_ind] = 0
        return volume_to_zero

    def _generate_Kspace_data(self):
        """
        Get the k space of the current slice data.
        :return:
        """
        self.k_space = fftshift(fft2(fftshift(self._image_to_correct)))

    def _calculate_encoding_indices(self):
        """
        generates both the linear (i.e. assumed) indices and the distorted indices.
        These indices are passed to the NUFFT library.
        """

        x_lin_size, y_lin_size = self._image_to_correct.shape
        xn_lin = np.linspace(-x_lin_size / 2, -x_lin_size / 2 + x_lin_size - 1, x_lin_size)
        yn_lin = np.linspace(-y_lin_size / 2, -y_lin_size / 2 + y_lin_size - 1, y_lin_size)
        [xn_lin, yn_lin] = np.meshgrid(xn_lin, yn_lin, indexing='ij')
        xn_lin = xn_lin.flatten()
        yn_lin = yn_lin.flatten()
        self.sk = xn_lin / x_lin_size
        self.tk = yn_lin / y_lin_size

        if (np.round(self._ImageOrientationPatient) == [0, 1, 0, 0, 0, -1]).all():
            xn_dis = self.Gz_encode / (self._pixel_spacing[2]) + self._iso_offset_pixels[2]
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gy_encode / (self._pixel_spacing[1]) + self._iso_offset_pixels[1]
            self.yj = yn_dis * 2 * np.pi
        elif (np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 0, -1]).all():
            xn_dis = (self.Gz_encode / (self._pixel_spacing[2])) + self._iso_offset_pixels[2]
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gx_encode / (self._pixel_spacing[0]) + self._iso_offset_pixels[0]
            self.yj = yn_dis * 2 * np.pi
        elif (np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 1, 0]).all():
            xn_dis = self.Gy_encode / (self._pixel_spacing[1]) + self._iso_offset_pixels[1]
            self.xj = xn_dis * 2 * np.pi
            yn_dis = self.Gx_encode / (self._pixel_spacing[0]) + self._iso_offset_pixels[0]
            self.yj = yn_dis * 2 * np.pi
        elif np.round(self._ImageOrientationPatient == [2, 2, 2, 2, 2, 2]).all():
            # this is for through plane correction where the real images are [0, 1, 0, 0, 0, -1]
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            yn_dis = self.Gx_encode / (self._pixel_spacing[0]) + self._iso_offset_pixels[0]
            self.yj = pd.Series(yn_dis * 2 * np.pi)
            xn_dis = pd.Series(np.copy(xn_lin))
        elif np.round(self._ImageOrientationPatient == [3, 3, 3, 3, 3, 3]).all():
            # this is for through plane correction where the real images are [1, 0, 0, 0, 0, -1]
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            yn_dis = self.Gy_encode / (self._pixel_spacing[1]) + self._iso_offset_pixels[1]
            self.yj = yn_dis * 2 * np.pi
            xn_dis = pd.Series(np.copy(xn_lin))
        elif (np.round(self._ImageOrientationPatient) == [1, 1, 1, 1, 1, 1]).all():
            # this is for through plane correction where the real images are [1, 0, 0, 0, 1, 0]
            self.xj = pd.Series(xn_lin * 2 * np.pi)
            yn_dis = -1 * self.Gz_encode / (self._pixel_spacing[2]) + self._iso_offset_pixels[2]
            self.yj = yn_dis * 2 * np.pi
            xn_dis = pd.Series(np.copy(xn_lin))
        else:
            raise NotImplementedError('this slice orientation is not handled yet sorry')

        self.xj = self.xj.to_numpy()
        self.yj = self.yj.to_numpy()
        self.xn_dis_pixel = np.reshape((xn_dis - xn_lin).values, self._image_to_correct.shape)
        self.yn_dis_pixel = np.reshape((yn_dis - yn_lin).values, self._image_to_correct.shape)

    def _plot_encoding_fields(self, vmin=None, vmax=None):  # pragma: no cover
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

    def _plot_indices(self):  # pragma: no cover
        """
        handy debug plot routine to plot linear/distorted indices
        :return:
        """

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[8, 8])
        axs[0, 0].plot(self.xj)
        axs[0, 0].set_title('xj')
        axs[0, 1].plot(self.yj)
        axs[0, 1].set_title('yj')
        axs[1, 0].plot(self.sk)
        axs[1, 0].set_title('sk')
        axs[1, 1].plot(self.tk)
        axs[1, 1].set_title('tk')

    def _plot_coords(self):  # pragma: no cover

        fig, axs = plt.subplots(1, 3)
        axs[0].plot(self._X_slice.flatten())
        axs[1].plot(self._Y_slice.flatten())
        axs[2].plot(self._Z_slice.flatten())
        plt.show()

    def _plot_images(self):  # pragma: no cover
        """
        debugging method to show input, output, and diff images
        :return:
        """
        vmin = 20
        vmax = 100
        fig, axs = plt.subplots(nrows=1, ncols=3)
        axs[0].imshow(self._image_to_correct, vmin=vmin, vmax=vmax)
        axs[1].imshow(self.outputImage, vmin=vmin, vmax=vmax)
        axs[2].imshow((abs(np.subtract(self._image_to_correct, self.outputImage))), vmin=vmin, vmax=vmax)
        axs[2].set_title('abs_difference')
        plt.tight_layout()

    # public methods

    def correct_all_images(self):
        """
        This loops through the image array, and calls _correct_image on each slice.
        If correct_through_plane is True, this process occurs twice, with the 'slice direction' being changed
        in the second loop.
        :return:
        """
        start_time = perf_counter()

        update_every_n_slices = 10
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

            if i % update_every_n_slices == 0:
                print(printProgressBar(i - self._n_zero_pad, n_images_to_correct))
            self._image_to_correct = array_slice
            self._X_slice = X
            self._Y_slice = Y
            self._Z_slice = Z
            self._correct_image()
            self._image_array_corrected[:, :, i] = self.outputImage
            i += 1

        if self.correct_through_plane:

            if np.allclose(np.round(self._ImageOrientationPatient), [1, 0, 0, 0, 1, 0]):
                self._ImageOrientationPatient = [1, 1, 1, 1, 1, 1]
            elif np.allclose(np.round(self._ImageOrientationPatient), [0, 1, -0, -0, 0, -1]):
                self._ImageOrientationPatient = [2, 2, 2, 2, 2, 2]
            elif np.allclose(np.round(self._ImageOrientationPatient), [1, 0, 0, 0, 0, -1]):
                self._ImageOrientationPatient = [3, 3, 3, 3, 3, 3]
            else:
                raise NotImplementedError

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

                if j % update_every_n_slices == 0:
                    print(printProgressBar(i + j - (self._n_zero_pad * 4), n_images_to_correct))
                self._image_to_correct = array_slice
                self._X_slice = X
                self._Y_slice = Y
                self._Z_slice = Z
                self._correct_image()
                self._image_array_corrected[:, j, :] = self.outputImage

                j += 1

        execution_time = perf_counter() - start_time
        print(f'\ntotal time: {execution_time: 1.1f}s')
        print(f'mean time per slice: {execution_time / n_images_to_correct: 1.1}s')
        self._unpad_image_arrays()


    def save_all_images(self, save_loc=None, DSV_radius=None, grid=True):
        """
        save corrected data as png

        :param save_loc: path to save data at.
        :type save_loc: string or path
        """

        import cProfile, pstats, io
        from pstats import SortKey
        pr = cProfile.Profile()
        pr.enable()

        _start_time = perf_counter()
        print('saving all data as png...')
        plt.ioff()
        if save_loc is None:
            save_loc = self.ImageDirectory / 'corrected'
        save_loc = Path(save_loc)
        if not (self.ImageDirectory / save_loc).is_dir():
            save_loc.mkdir()

        loop_axis = 2
        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._image_array_corrected, loop_axis))

        # generate extent and labels
        extent = None
        row_label = None
        col_label = None
        FOV = None
        slice_coords = None

        if ((np.round(self._ImageOrientationPatient) == [0, 1, 0, 0, 0, -1]).all() or
                np.round(self._ImageOrientationPatient == [2, 2, 2, 2, 2, 2]).all()):
            extent = (self._X.min(), self._X.max(), self._Z.min(), self._Z.max())
            slice_coords = np.unique(self._X)
            row_label = 'X [mm]'
            col_label = 'Z [mm]'

        elif ((np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 0, -1]).all() or
              np.round(self._ImageOrientationPatient == [3, 3, 3, 3, 3, 3]).all()):

            extent = (self._Y.min(), self._Y.max(), self._Z.min(), self._Z.max())
            slice_coords = np.unique(self._Y)
            row_label = 'Y [mm]'
            col_label = 'Z [mm]'

        elif ((np.round(self._ImageOrientationPatient) == [1, 0, 0, 0, 1, 0]).all() or
              np.round(self._ImageOrientationPatient == [1, 1, 1, 1, 1, 1]).all()):
            extent = (self._X.min(), self._X.max(), self._Y.min(), self._Y.max())
            slice_coords = np.unique(self._Z)
            row_label = 'X [mm]'
            col_label = 'Y [mm]'

        if slice_coords is None:
            slice_coords = np.ones(np.rollaxis(self.ImageArray, loop_axis).shape[0])
        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._image_array_corrected, loop_axis),
                          slice_coords)

        i = 0
        for original_image, corrected_image, slice_coord in zipped_data:
            fig, axs = plt.subplots(nrows=1, ncols=2)

            # axs[0].imshow(self.pixel_array, extent=self._extent)
            im1 = axs[0].imshow(original_image, extent=extent)
            axs[0].set_title('Original Image')
            axs[0].set_aspect('equal')
            axs[0].set_xlabel(row_label)
            axs[0].set_ylabel(col_label)
            if grid:
                axs[0].grid(True)
            if DSV_radius:
                h = 150 - slice_coord
                with warnings.catch_warnings(record=True) as who_cares:
                    warnings.simplefilter("always")
                    FOV_radius = np.sqrt((2 * h * DSV_radius) - (h ** 2))
                circ = Circle((0, 0), FOV_radius, facecolor='none', edgecolor='white')
                axs[0].add_patch(circ)

            im2 = axs[1].imshow(corrected_image, extent=extent)
            axs[1].set_title('Corrected Image')
            axs[1].set_aspect('equal')
            axs[1].set_xlabel(row_label)
            axs[1].set_ylabel(col_label)
            if grid:
                axs[1].grid(True)
            if DSV_radius:
                circ = Circle((0, 0), FOV_radius, facecolor='none', edgecolor='white')
                axs[1].add_patch(circ)
            plt.tight_layout()
            plt.savefig(save_loc / (str(i) + '.png'), format='png')
            plt.close('all')
            i += 1
        print(f'images export to png successful in {perf_counter() - _start_time} s')

        pr.disable()
        pr.dump_stats('save_images_stats')


    def save_all_images_as_dicom(self, save_loc=None):
        """
        save corrected data as dicom

        :param save_loc: path to save data at.
        :type save_loc: string or path
        """
        _start_time = perf_counter()
        print('saving all data as dcm...')
        if self._image_array_corrected.min() < 0:
            self._image_array_corrected = self._zero_volume(self._image_array_corrected)
        if save_loc is None:
            save_loc = self.ImageDirectory / 'corrected_dcm'
        save_loc = Path(save_loc)
        if not save_loc.is_dir():
            save_loc.mkdir()
        loop_axis = 2
        zipped_data = zip(np.rollaxis(self.ImageArray, loop_axis),
                          np.rollaxis(self._image_array_corrected, loop_axis))
        i = 0
        for original_image, corrected_image in zipped_data:
            current_slice = pydicom.read_file(self.ImageDirectory / self._all_dicom_files[i])
            temp_dcm = current_slice.copy()
            temp_dcm.PixelData = np.uint16(corrected_image).tobytes()
            temp_dcm.save_as(save_loc / (str(i) + '.dcm'))
            i += 1
        print(f'images exported to dicom in {perf_counter() - _start_time: 1.1f}s')


class KspaceDistortionCorrector(DistortionCorrectorBase):
    """
    This performs correction using a least squares based approach in the k-space domain,
    and is based on `this work <https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25487>`_
    """

    def __init__(self, pad=10, **kwds):

        super().__init__(pad=pad, **kwds)
        self.NUFFTlibrary = 'finufft'  # pynufft, finufft, or torchnufft
        self._check_input_data()

    def _generate_Kspace_data(self):
        """
        Get the k space of the current slice data.
        :return:
        """
        self.k_space = fftshift(fft2(fftshift(self._image_to_correct)))

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
        self._generate_Kspace_data()
        self._calculate_encoding_indices()
        self._perform_least_squares_optimisation()

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

        maxit = 20
        x1 = lsqr(A, fk1, iter_lim=maxit, x0=StartingImage)
        self.outputImage = abs(np.reshape(x1[0], [self._Rows, self._Cols]))


class ImageDomainDistortionCorrector(DistortionCorrectorBase):
    """
    This performs image correction in the image domain by constructing a spline between the linear and distorted
    indices.
    """

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
        self._calculate_encoding_indices()

        x_n = np.arange(0, self._image_shape[0], 1)
        y_n = np.arange(0, self._image_shape[1], 1)
        pixel_interp = RectBivariateSpline(x_n, y_n, self._image_to_correct)

        yy, xx = np.meshgrid(y_n, x_n)
        xvector = (self.xn_dis_pixel + xx).flatten()
        yvector = (self.yn_dis_pixel + yy).flatten()

        _output_image = pixel_interp.ev(xvector, yvector)

        _output_image = np.where((xvector < 0) | (yvector < 0) | (xvector > self._image_shape[0]) |
                                 (yvector > self._image_shape[1]), 0, _output_image)

        self.outputImage = _output_image.reshape(self._image_shape)
