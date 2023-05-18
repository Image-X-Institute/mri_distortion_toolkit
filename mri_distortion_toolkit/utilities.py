"""
General purpose reusable functions
"""

import os
import pathlib
import logging
from pathlib import Path
import numpy as np
from scipy.special import lpmv
import pandas as pd
from matplotlib import pyplot as plt
import pydicom
import json
from itertools import compress
import plotly.graph_objects as go
import warnings
from get_all_files import get_all_files
from datetime import datetime

ch = logging.StreamHandler()
formatter = logging.Formatter('[%(filename)s: line %(lineno)d %(levelname)8s] %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(ch)
logger.setLevel(logging.INFO)  # This toggles all the logging in your app
logger.propagate = False


class bcolors:
    """
    This is just here to enable me to print pretty colors to the linux terminal
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def build_dicom_affine(Dicomfiles):
    """
    Build the dicom affine matrix.
    The reference I used to create the affine is `NIbabel <https://nipy.org/nibabel/dicom/dicom_orientation.html>`_
    (we are using the r/c defined under 'DICOM affines again')

    :param Dicomfiles: list of pydicom.dataset.FileDataset instances (created using pydicom.read_file)
    :type Dicomfiles: list
    :return: CoordinateMatrix, a 4x4 matrix which transforms the image based coordinates (i,j,k) to world coordinates (x,y,z)
    """

    assert isinstance(Dicomfiles, list)
    assert isinstance(Dicomfiles[0], pydicom.dataset.FileDataset)

    ds = Dicomfiles[0]  # just use first one, they all give the same result
    if not np.logical_or(np.isclose(np.abs(ds.ImageOrientationPatient), 1),
                         np.isclose(np.abs(ds.ImageOrientationPatient), 0)).all():
        '''
        This just checks that the ImageOrientationPatient tag only contains 0 or 1, i.e. that no weird oblique
        orientation is in use. For now I don't want to support or think about non cardinal encoding directions
        '''
        raise NotImplementedError('at the moment we only support cardinal imaging planes')

    CoordinateMatrix = np.zeros([4, 4])
    F = np.array([ds.ImageOrientationPatient[3:], ds.ImageOrientationPatient[:3]]).T
    F = np.round(F)
    first_col = np.multiply(F[:, 0], ds.PixelSpacing[0])
    sec_col = np.multiply(F[:, 1], ds.PixelSpacing[1])
    CoordinateMatrix[:, 0][0:3] = first_col
    CoordinateMatrix[:, 1][0:3] = sec_col
    CoordinateMatrix[:, 3][0:3] = [float(i) for i in ds.ImagePositionPatient]
    if Dicomfiles.__len__() > 1:
        # figure out slice direction
        SliceLoc = []
        image_position_patient_list = []
        for file in Dicomfiles:
            # need to extract the slice coordinates one at a time (can't think of a better way anyway...)
            image_position_patient_list.append(file.ImagePositionPatient)
            SliceLoc.append(file.SliceLocation)
        SliceLoc = np.array(SliceLoc)
        image_position_patient_list = np.array(image_position_patient_list)
        if not np.allclose(np.diff(SliceLoc), np.diff(SliceLoc)[0], atol=.1, rtol=.1):
            logger.error(f'it appears there are missing or inconsistently spaced slices. Continuing')
        k = (image_position_patient_list[-1] - image_position_patient_list[0])/(image_position_patient_list.shape[0]-1)
        CoordinateMatrix[0:3, 2] = k

    return CoordinateMatrix


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar. `credit here <https://gist.github.com/greenstick/b23e475d2bfdc3a82e34eaa1f6781ee4>`_

    :param iteration: current iteration (Int)
    :param total: total iterations (Int)
    :param prefix:  prefix string (Str)
    :param suffix: suffix string (Str)
    :param decimals: positive number of decimals in percent complete (Int)
    :param length: character length of bar (Int)
    :param fill: bar fill character (Str)
    :return: progress_string: the string to print
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    progress_string = f'\r{prefix} |{bar}| {percent}% {suffix}'
    return progress_string


def dicom_to_numpy(path_to_dicoms, FilesToReadIn=None, file_extension='dcm', return_XYZ=False,
                   zero_pad=0, enforce_increasing_coords=False):
    """
    This function does two things:

    1. Loads all the dicom files in path_to_dicoms in a numpy array. data is ordered according to [row, column, slice]
    2. Builds the matrix that transforms the [r,c,s] indices into [x,y,z] coordinates

    :param path_to_dicoms: folder where all dicoms images are stored
    :type path_to_dicoms: pathlib.Path or string
    :param FilesToReadIn: By default, this function will read in all files of file type file_extension. If you only want
        to read in a subset of files, you can specify them with FilesToReadIn
    :type FilesToReadIn: list, optional
    :param file_extension: extension of dicom files, e.g. 'dcm' or 'IMA'
    :type file_extension: string
    :param return_XYZ: if True, will build and return three matrics containing the X, Y, and Z coords at each i,j,k
        pixel.
    :type return_XYZ: bool, optional
    :param zero_pad: this many zeros will be placed around the returned volume and the coordinates updated accordingly.
        main use case was for distortion correction. Note that all returned objects will be  consistently affected
        by this parameter
    :type zero_pad: int, optional
    :returns: ImageArray: a numpy array of voxels
    :returns: dicom_affine: a matrix that transforms the [r,c,s] indices into [x,y,z] coordinates
    :returns: (X, Y, Z): Coordinate arrays the same size as ImageArray. Optional
    """
    if not FilesToReadIn:
        FilesToReadIn = get_all_files(path_to_dicoms, file_extensions=file_extension)
    if not FilesToReadIn:
        raise FileNotFoundError(f'could not find any {file_extension} files in input directory {path_to_dicoms}'
                                f'\nTry using the ImExtenstion= input option to look for different image types')
    if isinstance(FilesToReadIn, list):
        CompletePathFiles = [str(Path(path_to_dicoms) / file) for file in FilesToReadIn]
    elif isinstance(FilesToReadIn, (str, pathlib.Path)):
        CompletePathFiles = [str(Path(path_to_dicoms) / FilesToReadIn)]

    dicom_slices = [pydicom.read_file(f) for f in CompletePathFiles]

    dicom_slices = sort_dicom_slices(dicom_slices)[0]
    dicom_affine = build_dicom_affine(dicom_slices)
    dicom_affine[0:3, 3] = dicom_affine[0:3, 3] + (-1*zero_pad*dicom_affine[0:3, 0:3].sum(axis=1))  # update start point for zero padding
    n_rows = dicom_slices[0].Rows + zero_pad*2
    n_cols = dicom_slices[0].Columns + zero_pad*2
    n_slices = len(dicom_slices) + zero_pad*2

    ImageArray = np.zeros([n_rows, n_cols, n_slices])
    if enforce_increasing_coords:
        # this option exists for the distortion correction code, which increasing indices
        new_affine = dicom_affine.copy()
        coordinates_end = dicom_affine @ [n_rows-1, n_cols-1, n_slices-1, 1]
        for col in [0, 1]:  # check first two cols
            for i, pixel_increment in enumerate(dicom_affine[:, col]):
                if pixel_increment < 0:
                    '''
                    then we need to force this to positive, and switch the start position to the 
                    current end position
                    '''
                    new_affine[i, col] = -1 * pixel_increment
                    new_start = coordinates_end[i]
                    new_affine[i, 3] = new_start
        dicom_affine = new_affine
    for i, file in enumerate(dicom_slices):
        if zero_pad > 0:
            ImageArray[zero_pad:-zero_pad, zero_pad:-zero_pad, i+zero_pad] = file.pixel_array
        else:
            ImageArray[:, :, i] = file.pixel_array
    if return_XYZ:
        rcsMatrix = np.zeros([4, n_rows * n_cols * n_slices])
        r_indices = np.arange(0, n_rows)
        c_indices = np.arange(0, n_cols)
        s_indices = np.arange(0, n_slices)
        [RR, CC, SS] = np.meshgrid(r_indices, c_indices, s_indices, indexing='ij')
        rcsMatrix[0, :] = RR.flatten()
        rcsMatrix[1, :] = CC.flatten()
        rcsMatrix[2, :] = SS.flatten()
        rcsMatrix[3, :] = 1
        XYZtemp = np.matmul(dicom_affine, rcsMatrix)
        X = np.reshape(XYZtemp[0, :], ImageArray.shape)
        Y = np.reshape(XYZtemp[1, :], ImageArray.shape)
        Z = np.reshape(XYZtemp[2, :], ImageArray.shape)

        assert ImageArray.shape == X.shape

        return ImageArray, dicom_affine, (X, Y, Z)
    else:
        return ImageArray, dicom_affine


def sort_dicom_slices(dicom_datasets, filter_non_image_data=True):
    """
    sort slices by instance number

    :param dicom_datasets: list of pydicom files
    :return: sorted list of pydicom files
    :param filter_non_image_data: if True, any dicom file not containing a pixel_data attribute is removed
    :type filter_non_image_data: boolean, optional
    """
    if filter_non_image_data:
        contains_pixel_data_ind = [hasattr(dicom, 'pixel_array') for dicom in dicom_datasets]
        dicom_datasets = list(compress(dicom_datasets, contains_pixel_data_ind))
    instance_number = [el.InstanceNumber for el in dicom_datasets]
    sort_ind = np.argsort(instance_number)
    return list(np.array(dicom_datasets)[sort_ind]), sort_ind

def _get_MR_acquisition_data(data_loc, X, Y, Z, file_extension='.dcm'):
    """
    If we have MR images, extract the information from thall_dicom_files = get_all_files(data_loc)e
    dicom header that is needed to convert marker positions to magnetic fields downstream::

        InputVolume, dicom_affine, (X, Y, Z) = dicom_to_numpy(distorted_data_loc,
                                                      file_extension='dcm',
                                                      return_XYZ=True)

        dicom_data = get_MR_acquisition_data(distorted_data_loc, X, Y, Z)

    :param data_loc: location of dicom files
    :param X, Y, Z: full 3D coordinate matrices as generated by dicom_to_numpy
    :param file_extension: extension of files, e.g. 'dcm', 'IMA'.
    :return: dicom_data: a dictionary of useful information about the acquisition
    """
    all_dicom_files = get_all_files(data_loc, file_extensions=file_extension)
    # use first one that has pixel data:
    j = 0
    while True:
        example_dicom_file = pydicom.read_file(data_loc / all_dicom_files[j])
        try:
            test = example_dicom_file.pixel_array
            break
        except AttributeError:
            j += 1

    if example_dicom_file.Modality == 'MR':
        if not example_dicom_file.Manufacturer == 'SIEMENS':
            print(
                'this code has not been tested with non-siemens scanners. If dicom standards work properly'
                'this code will too...')

        x = np.unique(X)
        y = np.unique(Y)
        z = np.unique(Z)
        x_pixel_spacing = np.mean(np.diff(x))
        y_pixel_spacing = np.mean(np.diff(y))
        z_pixel_spacing = np.mean(np.diff(z))

        x_fov = x.max() - x.min() + x_pixel_spacing
        y_fov = y.max() - y.min() + y_pixel_spacing
        z_fov = z.max() - z.min() + z_pixel_spacing

        _dicom_header = example_dicom_file

        dicom_data = {}  # fill up below
        dicom_data['FOV'] = [x_fov, y_fov, z_fov]
        dicom_data['bandwidth'] = float(example_dicom_file.PixelBandwidth)

        dicom_data['gama'] = float(example_dicom_file.MagneticFieldStrength * 42.57747851)
        dicom_data['pixel_spacing'] = [x_pixel_spacing, y_pixel_spacing, z_pixel_spacing]
        dicom_data['image_size'] = [x.size, y.size, z.size]
        dicom_data['magnetic_field_strength'] = example_dicom_file.MagneticFieldStrength
        dicom_data['imaging_frequency'] = example_dicom_file.ImagingFrequency
        try:
            acquisition_date = datetime.strptime(example_dicom_file.AcquisitionDate, '%Y%m%d')
            dicom_data['acquisition_date'] = acquisition_date.strftime("%d_%B_%Y")
        except ValueError as e:
            print(e)
            print('failed to extract datetime, probably didnt understand date format, continuing...')
        dicom_data['manufacturer'] = example_dicom_file.Manufacturer
        try:
            dicom_data['sequence_name'] = example_dicom_file.SequenceName
        except AttributeError:
            pass
        try:
            dicom_data['sequence_name'] = example_dicom_file.ScanningSequence
        except AttributeError:
            pass
        # the above works on siemens scanner, not sure if it will prove consistent on others
        fat_water_delta_f = example_dicom_file.MagneticFieldStrength * 42.57747851 * 3.5
        dicom_data['chem_shift_magnitude'] = (example_dicom_file.PixelBandwidth / fat_water_delta_f) / 2
        dicom_data['InPlanePhaseEncodingDirection'] = example_dicom_file.InPlanePhaseEncodingDirection
        # extract the frequency, phase, and slice directions (B0 effects occur in freq and slice)
        directions = ['x', 'y', 'z']
        phase_encode_direction = None
        freq_encode_direction = None
        if dicom_data['InPlanePhaseEncodingDirection'] == 'ROW':
            phase_cosines = example_dicom_file.ImageOrientationPatient[:3]
            row_cosines = example_dicom_file.ImageOrientationPatient[3:]
        elif dicom_data['InPlanePhaseEncodingDirection'] == 'COL':
            phase_cosines = example_dicom_file.ImageOrientationPatient[3:]
            row_cosines = example_dicom_file.ImageOrientationPatient[:3]
        for i, (phase_cosine, freq_cosine) in enumerate(zip(phase_cosines, row_cosines)):
            if abs(float(phase_cosine)) == 1:
                phase_encode_direction = directions[i]
            if abs(float(freq_cosine)) == 1:
                freq_encode_direction = directions[i]
        if phase_encode_direction is None:
            print('failed to extract phase encoding direction from dicom data, continuing')
        if freq_encode_direction is None:
            print('failed to extract frequency encoding direction from dicom data, continuing')
        dicom_data['phase_encode_direction'] = phase_encode_direction
        dicom_data['freq_encode_direction'] = freq_encode_direction
        # get slice direction
        slice_dir_ind = np.equal(example_dicom_file.ImageOrientationPatient[:3],
                                 example_dicom_file.ImageOrientationPatient[3:])
        dicom_data['slice_direction'] = np.array(directions)[slice_dir_ind][0]
        bandwidth = np.array(dicom_data['bandwidth'])
        image_size = np.array(dicom_data['image_size'])
        gama = np.array(dicom_data['gama'])
        FOV = np.array(dicom_data['FOV'])
        gradient_strength = bandwidth * image_size / (gama * 1e6 * FOV * 1e-3)  # unit(T / m)
        dicom_data['gradient_strength'] = list(gradient_strength)

    else:
        dicom_data = None

    return dicom_data


def convert_cartesian_to_spherical(InputCoords):
    """
    Converts cartesian coordinates [x,y,z] to spherical coordinates [r, azimuth, elevation].
    If r, azimuth, or elevation already exist as column names they are overwritten.
    The convention used here matches the convention
    `matlab uses <https://www.mathworks.com/help/matlab/ref/cart2sph.html>`_, except that we add pi to the azimuth and
    pi/2 to the elevation such that

    - 0 <= azimuth <= 2*pi
    - 0 <= elevation <= pi

    There was once a very good reason for doing this; I can't remember what it was anymore but I like to stick
    to tradition.

    :param InputCoords: a data frame containing (at least) columns x, y, z.
    :returns InputCoords: same data frame with additional columns [r, azimuth, elevation]
    """
    if 'r' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='r')
    if 'azimuth' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='azimuth')
    if 'elevation' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='elevation')

    InputCoords.insert(3, "r", np.sqrt(InputCoords.x ** 2 + InputCoords.y ** 2 + InputCoords.z ** 2))
    InputCoords.insert(4, "azimuth", np.arctan2(InputCoords.y, InputCoords.x) + np.pi)
    InputCoords.insert(5, "elevation", np.arctan2(InputCoords.z, np.sqrt(InputCoords.x ** 2 + InputCoords.y ** 2)) + np.pi/2)


    return InputCoords


def generate_harmonic_names(n_order):
    """
    generate the names of each harmonic. Used to label the columns in the harmonics dataframe

    :param n_order: order of harmonics
    :type n_order: int

    :return: coeef_names: a list of names
    """
    k = 0
    coeff_names = []
    for n in range(0, n_order + 1):  # the plus 1 is because range stops at -1 for some reason
        coeff_names.append(f'A_{n}_0')
        k = k + 1
        for m in range(0, n):
            coeff_names.append(f'A_{n}_{m + 1}')
            k = k + 1
            coeff_names.append(f'B_{n}_{m + 1}')
            k = k + 1
    return coeff_names


def combine_harmonics(harmonics_1, harmonics_2, operation='add', n_order_return=None):
    """
    This is a general purpose function to combine two harmonic series.
    If the series are of the same order, this is the default order of the returned series.
    If they are different, by default the returned series will match the higher order.
    In either case, you can specify the return order to be lower than the default.

    :param harmonics_1: first series to combine
    :type harmonics_1: array-like
    :param harmonics_2: second series to combine
    :type harmonics_2: array-like
    :param operation: 'add' or 'subtract'
    :type operation: str, optional
    :param n_order_return: order of returned series
    :type n_order_return: int, optional

    :return: combined_harmonics
    """

    assert isinstance(n_order_return, int) or n_order_return is None
    assert operation=='add' or operation=='subtract'

    if (len(harmonics_1) == len(harmonics_2)):
        if operation == 'add':
            combined_harmonics = np.add(harmonics_1, harmonics_2)
        elif operation == 'subtract':
            combined_harmonics = np.subtract(harmonics_1, harmonics_2)
        harmonics_zeros = harmonics_1.copy()
        harmonics_zeros[:] = 0
        new_index = generate_harmonic_names(len(harmonics_2))
    elif len(harmonics_1) > len(harmonics_2):
        combined_harmonics = list(np.zeros(len(harmonics_1)))
        lower_ind = len(harmonics_2)
        if operation == 'add':
            combined_harmonics[:lower_ind] = np.add(harmonics_2, harmonics_1[:lower_ind])
        elif operation == 'subtract':
            combined_harmonics[:lower_ind] = np.subtract(harmonics_2, harmonics_1[:lower_ind])
        combined_harmonics[lower_ind:] = harmonics_1[lower_ind:]
    else:
        combined_harmonics = list(np.zeros(len(harmonics_1)))
        lower_ind = len(harmonics_1)
        if operation == 'add':
            combined_harmonics[:lower_ind] = np.add(harmonics_2[:lower_ind], harmonics_1)
        elif operation == 'subtract':
            combined_harmonics[:lower_ind] = np.subtract(harmonics_2[:lower_ind], harmonics_1)
        combined_harmonics[lower_ind:] = harmonics_2[lower_ind:]

    n_order_in = int(np.sqrt(len(combined_harmonics)) - 1)
    new_index = generate_harmonic_names(n_order_in)
    assert len(combined_harmonics) == np.max([len(harmonics_2), len(harmonics_1)])
    if n_order_return is not None:

        if n_order_return > n_order_in:
            logger.warning(f'return order of {n_order_return} specified, but maximum input order is {n_order_in} '
                          f'so this has no effect')
            clip_ind = len(combined_harmonics)
        else:
            clip_ind = int((n_order_return + 1) ** 2)
        combined_harmonics = combined_harmonics[:clip_ind]
        new_index = new_index[:clip_ind]

    combined_harmonics = pd.Series(combined_harmonics, index=new_index)

    return combined_harmonics


def convert_spherical_to_cartesian(InputCoords):
    """
    Converts spherical coordinates [r, azimuth, elevation] to cartesian coordinates [x,y,z].
    if the cartesian columns 'x', 'y', 'z' already exist, they will be overwritten

    :param InputCoords: a data frame containing columns r, azimuth, elevation
    :returns InputCoords:   same data frame with additional columns [x,y,z]
    """
    if 'x' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='x')
    if 'y' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='y')
    if 'z' in InputCoords.columns:
        InputCoords = InputCoords.drop(columns='z')

    # InputCoords.insert(3, "x", InputCoords.r * np.sin(InputCoords.azimuth - np.pi) * np.cos(InputCoords.elevation - np.pi/2))
    # InputCoords.insert(4, "y", InputCoords.r * np.sin(InputCoords.azimuth - np.pi) * np.sin(InputCoords.elevation - np.pi/2))
    # InputCoords.insert(5, "z", InputCoords.r * np.cos(InputCoords.azimuth))

    InputCoords.insert(3, "x", InputCoords.r * np.cos(InputCoords.azimuth - np.pi) * np.cos(InputCoords.elevation - np.pi/2))
    InputCoords.insert(4, "y", InputCoords.r * np.sin(InputCoords.azimuth - np.pi) * np.cos(InputCoords.elevation - np.pi/2))
    InputCoords.insert(5, "z", InputCoords.r * np.sin(InputCoords.elevation - np.pi/2))

    return InputCoords


def generate_legendre_basis(InputCoords, n_order):
    """
    Produce spherical base at all coordinate values.
    This utilises the 'fully normalised' convention described here <https://au.mathworks.com/help/matlab/ref/legendre.html>`_
    Note that we scale all r values by 1000, which is to help ensure numerical stability. This is based on the assumption
    that this input data is in mm. The code will throw a warning if this does not appear to the case.

    :param InputCoords: a data frame which can contain whatever you want, but MUST contain columns called
        elevation, azimuth, and r
    :type InputCoords: Dataframe
    :param n_order: order of harmonic fit
    :type n_order: int
    :returns legendre_basis: a pandas data frame of size [Ncoords, (n_order +1)**2)]
    """

    if InputCoords.r.mean() < 10:
        logger.warning('it appears that your input data is in m, not mm - please use mm!')

    legendre_basis = np.zeros([len(InputCoords.index), (n_order+1)**2])
    k = 0
    coeff_names = []
    for n in range(0, n_order+1):  # the plus 1 is because range stops at -1 for some reason
        Norm = np.sqrt(n+0.5)
        P = Norm * lpmv(0, n, np.cos(InputCoords.elevation))
        legendre_basis[:, k] = np.multiply((InputCoords.r/1000) ** n, P)  # m=0 terms
        coeff_names.append(f'A_{n}_0')
        k = k + 1
        for j in range(0, n):
            m = j+1  # this looks weird but it works to give values m=1:n
            Norm = ((-1)**m) * np.sqrt( (n+0.5)*(np.math.factorial(n-m)) / np.math.factorial(n+m))
            P = Norm * lpmv(m, n, np.cos(InputCoords.elevation))
            legendre_basis[:,k] = np.multiply(np.multiply((InputCoords.r/1000) ** n, P), np.cos(m * InputCoords.azimuth))
            coeff_names.append(f'A_{n}_{m}')
            k = k + 1
            legendre_basis[:, k] = np.multiply(np.multiply((InputCoords.r/1000) ** n, P), np.sin(m * InputCoords.azimuth))
            coeff_names.append(f'B_{n}_{m}')
            k = k + 1

    return pd.DataFrame(legendre_basis, columns=coeff_names)


def convert_spherical_harmonics(harmonics, input_format='full', output_format='none'):
    """
    function to convert between different harmonic formalisms.
    At the moment only works between 'fully normalised' and 'unnormalised' but we could add more later if
    we need to.
    example::

        g_x_harmonics_no_norm = convert_spherical_harmonics(g_x_harmonics, input_format='full', output_format='none')

    :param harmonics: a pandas series or numpy array of size (n_order+1)^2, where the harmonics are ordered as
        A00, A10, A11, B11, A20, A21, B21, A22, B22...etc.
    :param input_format: format of input harmonic normalisation; 'none' or 'full'
    :param output_format: format of output harmonic normalisation; 'none' or 'full'
    :return: converted harmonics as a pandas series
    """

    formats = ['none', 'full']
    if input_format not in formats:
        raise NotImplementedError(f'no method exists to convert {input_format} format')
    if output_format not in formats:
        raise NotImplementedError(f'no method exists to convert {output_format} format')

    case_string = input_format.lower() + '_' + output_format.lower()
    if input_format.lower == output_format.lower():
        case_string = 'same'  # will do nothing, except return the data as a series

    converted_harmonics = harmonics.copy()
    n_order = np.sqrt(converted_harmonics.shape[0]) - 1
    assert n_order.is_integer()  # if this is false then an invalid number of harmonics is entered
    n_order = int(n_order)
    # generate normalisation factors:
    k = 0
    coeff_names = []
    Norm = []
    for n in range(0, n_order + 1):  # the plus 1 is because range stops at -1 for some reason
        Norm.append(np.sqrt(n + 0.5))
        coeff_names.append(f'A_{n}_0')
        k = k + 1
        for j in range(0, n):
            m = j + 1  # this looks weird but it works to give values m=1:n
            Norm.append(((-1) ** m) * np.sqrt((n + 0.5) * (np.math.factorial(n - m)) / np.math.factorial(n + m)))
            coeff_names.append(f'A_{n}_{m}')
            coeff_names.append(f'B_{n}_{m}')
            Norm.append(((-1) ** m) * np.sqrt((n + 0.5) * (np.math.factorial(n - m)) / np.math.factorial(n + m)))

    if case_string == 'full_none':
        converted_harmonics = np.multiply(harmonics, Norm)
    if case_string == 'none_full':
        converted_harmonics = np.divide(harmonics, Norm)
    if isinstance(converted_harmonics, np.ndarray):
        converted_harmonics = pd.DataFrame(converted_harmonics, index=coeff_names).squeeze()

    return converted_harmonics


def get_harmonics(Gx_Harmonics, Gy_Harmonics, Gz_Harmonics, B0_Harmonics=None):
    """
    return the gradient spherical harmonics as pandas series. this function is simply a clean way to handle
    the different ways users can specify harmonics to other componets of this code

    :param Gx_Harmonics: either a pandas series or a path to a csv. If it is already a series we do nothing, but handle
        this option for cleaner code elsewhere. Gy_Harmonics and Gz_Harmonics are the same.
    :returns: Gx_Harmonics, Gy_Harmonics, Gz_Harmonics as Pandase series
    """

    if isinstance(Gx_Harmonics, pd.Series):
        Gx_Harmonics = Gx_Harmonics
    elif isinstance(Gx_Harmonics, (str,Path)):
        if not (Path(Gx_Harmonics).is_file()):
            raise FileNotFoundError(f'{Gx_Harmonics} does not exist')
        Gx_Harmonics = pd.read_csv(Gx_Harmonics, index_col=0).squeeze("columns")
    else:
        raise AttributeError('could not read in Gx harmonics...please input either a series or a '
                     'path to a saved csv file')

    if isinstance(Gy_Harmonics, pd.Series):
        Gy_Harmonics = Gy_Harmonics
    elif isinstance(Gy_Harmonics, (str,Path)):
        if not (Path(Gy_Harmonics).is_file()):
            raise FileNotFoundError(f'{Gy_Harmonics} does not exist')
        Gy_Harmonics = pd.read_csv(Gy_Harmonics, index_col=0).squeeze("columns")
    else:
        raise AttributeError('could not read in Gy harmonics...please input either a series or a '
                     'path to a saved csv file')

    if isinstance(Gz_Harmonics, pd.Series):
        Gz_Harmonics = Gz_Harmonics
    elif isinstance(Gz_Harmonics, (str,Path)):
        if not (Path(Gz_Harmonics).is_file()):
            raise FileNotFoundError(f'{Gz_Harmonics} does not exist')
        Gz_Harmonics = pd.read_csv(Gz_Harmonics, index_col=0).squeeze("columns")
    else:
        raise AttributeError('could not read in Gz harmonics...please input either a series or a '
                     'path to a saved csv file')

    if B0_Harmonics is not None:
        if isinstance(B0_Harmonics, pd.Series):
            B0_Harmonics = B0_Harmonics
        elif isinstance(B0_Harmonics, (str, Path)):
            if not (Path(B0_Harmonics).is_file()):
                raise FileNotFoundError(f'{B0_Harmonics} does not exist')
            B0_Harmonics = pd.read_csv(B0_Harmonics, index_col=0).squeeze("columns")
        else:
            raise AttributeError('could not read in B0 harmonics...please input either a series or a '
                                 'path to a saved csv file')

    return Gx_Harmonics, Gy_Harmonics, Gz_Harmonics, B0_Harmonics


def get_dicom_data(dicom_data):
    """
    figures out whether dicom data is a dict or a path to a json file
    if the latter, reads into  a dict and returns
    """

    if isinstance(dicom_data, dict):
        return dicom_data
    elif isinstance(dicom_data,(pathlib.Path, str)):
        with open(dicom_data, 'r') as f:
            dicom_data = json.load(f)
        return dicom_data
    else:
        raise AttributeError(f'could not read in dicom_data of type {type(dicom_data)}')


def reconstruct_Bz(harmonics, coords, quantity='uT', r_outer=None):
    """
    Reconstruct Bz at an arbitrary set of coordinates.
    n_order of reconstruction is derived from the shape of harmonics

    :param harmonics: a pandas series of size (n_order+1)^2, such as that defined within
        FieldAnalysis.SphericalHarmonicFit
    :param coords: a pandas Dataframe which can contain whatever you want, but MUST contain columns called:
        elevation, azimuth, and r
    :param quantity:
        - 'T': return raw field in T
        - 'uT': return field in uT, with A00 set to 0.
    :param r_outer: optional parameter; if you try to reconstruct points with r > r_outer a warning is thrown
    :type r_outer: float or None
    :returns Bz_recon: Pandas series with the value of Bz at each element in coords
    """

    if not quantity in ['uT', 'T']:
        logger.error(f'unsupported quantity "{quantity}" entered for recon.'
                       f'\n valid options are "field" or "uT". Using uT')
        quantity = 'uT'

    if any(abs(coords.r) > r_outer):
        logger.warning(f'you are trying to reconstruct {np.count_nonzero(abs(coords.r) > r_outer)} '
                       f'data points outside rOuter. You can do it, but the futher away'
                       'from rOuter you get the less you should trust this.')

    n_order = int(np.sqrt(harmonics.shape[0])-1)
    legendre_recon = generate_legendre_basis(coords, n_order)
    if quantity == 'uT':  # field in uT, with A00 set to 0 (e.g. only looking at perturbation)
        harmonics_temp = harmonics.copy()
        harmonics_temp.A_0_0 = 0
        Bz_recon = legendre_recon @ harmonics_temp
        Bz_recon = Bz_recon * 1e6
    elif quantity == 'T':  # raw field in T
        Bz_recon = legendre_recon @ harmonics
    else:
        raise NotImplementedError(f'unknown reconstruction option entered: {quantity}')

    return Bz_recon


def compare_recon_report_with_ground_truth_report(ground_truth_report, recon_report, tolerance=1):  # pragma: no cover
    """
    This is to compare the distortion data contained within two Reports.MRI_QA_recon_reporter
    objects. The original purpose was to compare a recon_report generated from direct data, and one
    generated via harmonics, although in principle it can take any two reports.

    :param ground_truth_report: Reports.MRI_QA_recon_reporter instance
    :param recon_report: Reports.MRI_QA_recon_reporter instance
    :param tolerance: tolerance in mm
    """
    # direct overlay of real versus recon
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    if ground_truth_report.r_outer == recon_report.r_outer:
        # then use r_outer to filter the data
        gt_data = ground_truth_report._extract_data_from_MatchedMarkerVolume(r_max=ground_truth_report.r_outer)
        recon_data = recon_report._extract_data_from_MatchedMarkerVolume(r_max=recon_report.r_outer)
    else:
        gt_data = ground_truth_report._extract_data_from_MatchedMarkerVolume()
        recon_data = recon_report._extract_data_from_MatchedMarkerVolume()

    axs[0, 0].scatter(gt_data.r_gt, gt_data.abs_dis)
    axs[0, 0].scatter(recon_data.r_gt, recon_data.abs_dis)
    axs[0, 0].legend(['original_data', 'recon_data'])
    axs[0, 0].set_xlabel('r [mm]')
    axs[0, 0].set_ylabel('absolute distortion [mm]')

    # scatter of real versus recon
    axs[0, 1].scatter(gt_data.abs_dis, recon_data.abs_dis)
    axs[0, 1].set_xlabel('original data')
    axs[0, 1].set_ylabel('recon data')

    # box plot of errors
    error = gt_data.abs_dis - recon_data.abs_dis
    axs[1, 0].boxplot(error)
    axs[1, 0].set_ylabel('recon error [mm]')


    # get index of outliers
    outliers_ind = error > tolerance
    print(f'{np.count_nonzero(outliers_ind)} outliers exist, at index {np.where(outliers_ind)}')
    outliers_ind = np.where(outliers_ind)
    for ind in outliers_ind[0]:
        print(f'index {ind} is at {gt_data.x_gt.iloc[ind]: 1.1f},'
              f' {gt_data.y_gt.iloc[ind]: 1.1f},'
              f' {gt_data.z_gt.iloc[ind]: 1.1f}  the value in the ground truth is'
              f' {gt_data.abs_dis.iloc[ind]: 1.1f}'
              f' and the value in the reconstructed is {recon_data.abs_dis.iloc[ind]: 1.1f}')


def enumerate_subfolders(data_loc):
    """
    A simple function that prints all the subfolders in data_loc as a dict, e.g::

        {'1': 'folder1',
         '2': 'folder2'} etc

    I find this is a useful to create a dict object that I copy to the start of my analysis scripts
    """
    data_dict = {}
    subfolders = [ f.path for f in os.scandir(data_loc) if f.is_dir() ]
    subfolders = np.sort(subfolders)
    for i, dir in enumerate(subfolders):
        folder_name = os.path.split(dir)[1]
        data_dict[str(i)] = folder_name
    for key in data_dict.keys():
        print(f"'{key}': '{data_dict[key]}',")


def plot_MarkerVolume_overlay(MarkerVolumeList, legend=None):  # pragma: no cover
    """
    Plot overlaid 3D scatter plots of the marker positions in each MarkerVolume

    :param MarkerVolumeList: list of MarkerVolume instances
    :param legend: legend to display on plot
    :return: None
    """

    assert isinstance(MarkerVolumeList, list)

    fig = plt.figure()
    axs = fig.add_subplot(111, projection='3d')
    for MarkerVolume in MarkerVolumeList:
        axs.scatter(MarkerVolume.MarkerCentroids.x, MarkerVolume.MarkerCentroids.y, MarkerVolume.MarkerCentroids.z)
        axs.set_xlabel('X [mm]')
        axs.set_ylabel('Y [mm]')
        axs.set_zlabel('Z [mm]')
        axs.set_title('3D marker positions')
    axs.set_box_aspect((np.ptp(MarkerVolume.MarkerCentroids.x),
                        np.ptp(MarkerVolume.MarkerCentroids.y),
                        np.ptp(MarkerVolume.MarkerCentroids.z)))
    if legend:
        plt.legend(legend)
    plt.show()


def plot_compressed_MarkerVolumes(MarkerVolumeList, p_max=20, p_min=-20, title=None, legend=None,
                                  projection_direction='z'): # pragma: no cover
    """
    plot overlay of compressed marker volumes (z coord is ignored)

    :param MarkerVolumeList:  list of MarkerVolume instances
    :type MarkerVolumeList: list
    :param p_max: max data to include in projection_direction
    :type p_max: float, optional
    :param p_min: min z data to include in projection_direction
    :type p_min: float, optional
    :param title: title of plot
    :type title: str, optional
    :param legend: legend of plot
    :type legend: list, optional
    :return: None
    """

    fig, axs = plt.subplots(figsize=[8, 8], ncols=1, nrows=1)
    if not isinstance(MarkerVolumeList, list):
        MarkerVolumeList = list(MarkerVolumeList)

    if not projection_direction in ['x', 'y', 'z']:
        warnings.warn(f'projection direction must be "x", "y", or "z", not {projection_direction}.'
                      f'setting to z and continuing')
        projection_direction = 'z'

    for MarkerVolume in MarkerVolumeList:
        if projection_direction == 'z':
            data_ind = np.logical_and(MarkerVolume.MarkerCentroids.z > p_min, MarkerVolume.MarkerCentroids.z < p_max)
            plot_data = MarkerVolume.MarkerCentroids[data_ind]
            axs.scatter(plot_data.x, plot_data.y)
            axs.set_xlabel('x [mm]', fontsize=15)
            axs.set_ylabel('y [mm]', fontsize=15)
        elif projection_direction == 'x':
            data_ind = np.logical_and(MarkerVolume.MarkerCentroids.x > p_min, MarkerVolume.MarkerCentroids.x < p_max)
            plot_data = MarkerVolume.MarkerCentroids[data_ind]
            axs.scatter(plot_data.z, plot_data.y)
            axs.set_xlabel('z [mm]', fontsize=15)
            axs.set_ylabel('y [mm]', fontsize=15)
        elif projection_direction == 'y':
            data_ind = np.logical_and(MarkerVolume.MarkerCentroids.y > p_min, MarkerVolume.MarkerCentroids.y < p_max)
            plot_data = MarkerVolume.MarkerCentroids[data_ind]
            axs.scatter(plot_data.x, plot_data.z)
            axs.set_xlabel('x [mm]', fontsize=15)
            axs.set_ylabel('z [mm]', fontsize=15)

        axs.set_xlim([-150, 150])
        axs.set_ylim([-150, 150])
        axs.tick_params(axis='both', which='major', labelsize=12)
        axs.tick_params(axis='both', which='minor', labelsize=12)
        axs.axis("equal")
        axs.grid()
        plt.subplots_adjust(bottom=0.25)
        if legend:
            plt.legend(legend)
        if title:
            plt.title(title)


def plot_matched_volume_hist(VolumeList, legend=None):
    """
    create a histogram of total distortion for all input volumes

    :param VolumeList: a list of [matched marker volumes](https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MatchedMarkerVolumes)
    :type VolumeList: list
    :param legend: legend to print
    :type legend: list, optional
    :return: None
    """
    bins = np.linspace(0, 10, 30)
    plt.figure()
    for marker_volume in VolumeList:
        plt.hist(marker_volume.MatchedCentroids.match_distance, bins=bins, alpha=0.5)
    if legend:
        plt.legend(['original', 'corrected'])
    plt.xlabel('distortion [mm]')
    plt.tight_layout()
    plt.show()


def plot_distortion_xyz_hist(MatchedMarkerVolume):
    """
    plot overlaid x, y, z distortion for an input instance of
    [MatchedMarkerVolume](https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MatchedMarkerVolumes)

    :param MatchedMarkerVolume: an instance of MatchedMarkerVolumes
    """

    x_dis = MatchedMarkerVolume.MatchedCentroids.x_gt - MatchedMarkerVolume.MatchedCentroids.x_gnl
    y_dis = MatchedMarkerVolume.MatchedCentroids.y_gt - MatchedMarkerVolume.MatchedCentroids.y_gnl
    z_dis = MatchedMarkerVolume.MatchedCentroids.z_gt - MatchedMarkerVolume.MatchedCentroids.z_gnl

    plt.figure()
    bins = np.linspace(0, 10, 30)
    plt.hist(abs(x_dis), bins=bins, alpha=0.5)
    plt.hist(abs(y_dis), bins=bins, alpha=0.5)
    plt.hist(abs(z_dis), bins=bins, alpha=0.5)
    plt.legend(['x', 'y', 'z'])
    plt.xlabel('distortion [mm]')
    plt.tight_layout()
    plt.show()


def print_dict(dict):
    """
    just provides a nicer display of the input dict than the default print function,
    got sick of writing the same code over and over...
    """
    for key in dict.keys():
        print(f'{key}: {dict[key]}')


def plot_harmonics_over_sphere(harmonics, radius, title=None):  # pragma: no cover
    """
    plot the reconstructed harmonics over the surface of a sphere
    I've used plotly instead of matplotlib for this as I wasnt happy with the matplotlib result.

    :param harmonics: either an instance of SphericalHarmonicFit or the equivalent data frame loaded from csv
    :param radius: radius of sphere
    :type radius: float
    :return: None
    """
    if title is None:
        title = 'DSV surface [uT]'
    # generate coordinates:
    _density = 100
    azimuth = np.linspace(0, 2*np.pi, _density)
    elevation = np.linspace(0, np.pi, _density)
    [AZ, EL] = np.meshgrid(azimuth, elevation)
    R = np.ones(AZ.shape) * radius
    coords = pd.DataFrame({'elevation': EL.flatten(), 'azimuth': AZ.flatten(), 'r': R.flatten()})
    coords = convert_spherical_to_cartesian(coords)
    X = coords.x.to_numpy().reshape(R.shape)
    Y = coords.y.to_numpy().reshape(R.shape)
    Z = coords.z.to_numpy().reshape(R.shape)

    Bz = reconstruct_Bz(harmonics, coords)
    Bz = Bz.to_numpy().reshape(R.shape)

    fig_harmonic_surface = go.Figure(data=[go.Surface(z=Z, x=X, y=Y, colorbar=dict(title='B(uT)',
                                                                                   titlefont={'size': 20},
                                                                                   tickfont=dict(size=40)),
                                                                    surfacecolor=Bz)])

    _available_plotly_themes = ["plotly", "plotly_white", "plotly_dark", "ggplot2",
                                "seaborn", "simple_white", "none"]
    _plotly_theme = 'seaborn'
    assert _plotly_theme in _available_plotly_themes
    fig_harmonic_surface.update_layout(title=title, autosize=True,
                                       template=_plotly_theme,
                                       font={'size': 20},
                                       scene=dict(xaxis=dict(title='X [mm]', titlefont=dict(size=20)),
                                                  yaxis=dict(title='Y [mm]', titlefont=dict(size=20)),
                                                  zaxis=dict(title='Z [mm]', titlefont=dict(size=20))))

    return fig_harmonic_surface


def read_opera_table_file(data_loc):
    """
    reads a table file from opera Tosca and returns the data as a data frame that can be
    read into spherical harmonic fit.
    Expected column order is  X [MM], Y [MM], Z [MM], BZ [T], with 6 header lines.
    No error handling is included for incorrect input

    :param data_loc: location of table file
    :type data_loc: str or pathlib.Path
    :return: Bz_data, a pandas dataframe
    """
    with open(data_loc, 'r') as f:
        data = np.loadtxt(f, skiprows=6)
    assert data.shape[1] == 4

    input_Bz_data = pd.DataFrame({'x': data[:, 0], 'y': data[:, 1], 'z': data[:, 2], 'Bz': data[:, 3]})
    return input_Bz_data