import pathlib
import numpy as np
import logging
import json
from .utilities import get_dicom_data
from .utilities import dicom_to_numpy
import pydicom
from pathlib import Path
from .utilities import get_all_files

ch = logging.StreamHandler()
formatter = logging.Formatter('[%(filename)s: line %(lineno)d %(levelname)8s] %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(ch)
logger.setLevel(logging.INFO)
logger.propagate = False


class ConvertMatchedMarkersToBz:
    """
    an object to calculate and store magnetic fields.
    Gradient fields are in T/m
    B0 fields are in T

    :param MatchedCentroids: a pandas data frame containing columns ['x_gt', 'y_gt', 'z_gt', 'x_gnl', 'y_gnl', 'z_gnl'].
        If you want to calculate B0, is should also have ['x_B0', 'y_B0', 'z_B0']. Such a dataframe is normally created
        within the MatchedMarkerVolumes class. All coordinates in mm.
    :param dicom_data: either a dictionary or a path to a .json file. The minimume required parameters are below.
        dicom_data is calculated automatically when MR data is used to create a MarkerVolume, and can be exported
        to json using::

            marker_volume.save_dicom_data().

        dicom_data example::

            dicom_data = {'FOV': [400, 400, 400],
            'bandwidth': 203,
            'gama': 42.6,
            'image_size': [256, 256, 256]}
    :type dicom_data: dictionary or path

    """

    def __init__(self, MatchedCentroids, dicom_data):
        """
        init
        """

        self.dicom_data = get_dicom_data(dicom_data)
        self.MatchedCentroids = MatchedCentroids
        self.MagneticFields = self.MatchedCentroids[['x_gt', 'y_gt', 'z_gt']].copy()
        self.MagneticFields = self.MagneticFields .rename(columns={"x_gt": "x", "y_gt": "y", "z_gt": "z"})
        # ^ aim is to calculate B field at the ground truth marker positions
        if self.MatchedCentroids.r_gt.mean() < 10 or self.MatchedCentroids.r_gnl.mean() < 10:
            logger.warning('it appears that your input data is in m, not mm - please use mm!')
        self._check_dicom_data()
        self._calculate_Bz()

    def _check_dicom_data(self):
        """
        check that the minium data we will be needing is there.
        """
        expected_keys = ['FOV', 'bandwidth', 'gama', 'pixel_spacing', 'image_size', 'gradient_strength']
        for key in expected_keys:
            if not key in self.dicom_data.keys():
                raise AttributeError(f'missing required dicom_data field: {key}')

    def _get_B0_calc_direction(self):
        """
        if we have information about the frequency encode direction and information about B0 distortion,
        we should find that the max B0 distortion occurs in the frequency direction.
        If not, something problably is wrong.
        in any case, the direction of maximum B0 distortion is returned and used to calcualte B0
        """

        x_B0_mean = np.mean(np.abs(self.MatchedCentroids.x_B0))
        y_B0_mean = np.mean(np.abs(self.MatchedCentroids.y_B0))
        z_B0_mean = np.mean(np.abs(self.MatchedCentroids.z_B0))
        max = np.max([x_B0_mean, y_B0_mean, z_B0_mean])
        directions = ['x', 'y', 'z']
        direction = [int(el == max) for el in [x_B0_mean, y_B0_mean, z_B0_mean]]
        self._B0_direction_bool = [bool(el) for el in direction]
        if 'freq_encode_direction' in self.dicom_data.keys():
            # check that the maximum distortion matches the frequence encode direction
            if not directions[np.where(direction)[0][0]] == self.dicom_data['freq_encode_direction']:
                logger.warning(f'\nYou have entered a frequency code direction in '
                               f'{self.dicom_data["freq_encode_direction"]}, but maximum B0 distortion is detected'
                               f'in {directions[np.where(direction)[0][0]]}. Proceeding using '
                               f'{self.dicom_data["freq_encode_direction"]} but you should check this...\n')
            self._B0_direction_bool = [self.dicom_data["freq_encode_direction"]==el for el in directions]
        else:
            logger.warning("\nyou have input data on B0 distortion, but dicom_data['freq_encode_direction'] does not exist."
                           "Therefore, B0 fields will be calculated based on the direction of maximum B0 distortion, "
                           "which is \n")
        self._B0_direction_string = np.array(['x_B0', 'y_B0', 'z_B0'])[self._B0_direction_bool][0]

    def _calculate_Bz(self):
        """
        Calculate Gradient Bz at each point. This is based on code from Shanshan.
        Shanshans notes:
        - When reversing gradient polarity, the readout (z) and phase encoding (y) directions are reversed,
         however, the slice selction direction (x) is not.
        - B0 inhomogeneity distortion appears at readout diretion and is neglectable at phase encoding direction.
        - This is distorted marker position caused by GNL-only along readout and phase encoding directions
        """
        gradient_strength = np.array(self.dicom_data['gradient_strength'])  # unit(T / m)
        # ^ this is a vector [gx, gy, gz]
        mm_to_m_factor = 1e-3
        self.MagneticFields['B_Gx'] = self.MatchedCentroids.x_gnl * gradient_strength[0] * mm_to_m_factor
        self.MagneticFields['B_Gy'] = self.MatchedCentroids.y_gnl * gradient_strength[1] * mm_to_m_factor
        self.MagneticFields['B_Gz'] = self.MatchedCentroids.z_gnl * gradient_strength[2] * mm_to_m_factor

        # if possible, add in B0 to data
        B0_fields = ['x_B0', 'y_B0', 'z_B0']
        B0_calc_possible = all([el in self.MatchedCentroids.columns for el in B0_fields])
        if B0_calc_possible:
            self._get_B0_calc_direction()
            self.MagneticFields['B0'] = self.MatchedCentroids[[self._B0_direction_string]] * \
                                        gradient_strength[self._B0_direction_bool][0] * mm_to_m_factor


class ConvertPhaseMapsToBz:
    """
    accepts two phase map images taken with different echo time, and returns the Bz field calculation
    """
    def __init__(self, image_echo1_loc, image_echo2_loc, te1=None, te2=None, file_extension='dcm'):

        self._image_echo1_loc = Path(image_echo1_loc)
        self._image_echo2_loc = Path(image_echo2_loc)
        self._file_extension = file_extension

        self._echo1, self._dicom_affine1, (self.X, self.Y, self.Z) = dicom_to_numpy(image_echo1_loc,
                                                                                    file_extension=file_extension,
                                                                                    return_XYZ=True)
        self._echo2, self._dicom_affine2, (self.X, self.Y, self.Z) = dicom_to_numpy(image_echo2_loc,
                                                                                    file_extension=file_extension,
                                                                                    return_XYZ=True)
        # assumption is image coordinates are the same; check:
        assert np.isclose(self._dicom_affine1, self._dicom_affine2).all()
        self._get_echo_time(te1, te2)
        self._unwrap_phase()

    def _unwrap_phase(self):
        """
        unwrap the two input phase images
        """
        from skimage.restoration import unwrap_phase
        self.phase1 = unwrap_phase(self._echo1)
        self.phase2 = unwrap_phase(self._echo2)

    def _calculate_field(self):
        """
        see this paper
        https://pubmed.ncbi.nlm.nih.gov/19810464/
        equation 2
        """
        print('ger')

    def _get_echo_time(self, te1, te2):
        """
        read echo time of the dicoms in data_loc and set it as a new attribute in self
        note that we only check one echo time; assumption is that folders are appropriately sorted
        """
        if te1 is None:
            dicom_files = get_all_files(self._image_echo1_loc, file_extension=self._file_extension)
            self._te1 = float(pydicom.read_file(self._image_echo1_loc / dicom_files[0]).EchoTime) *1e-3
        else:
            self._te1 = te1
        if te2 is None:
            dicom_files = get_all_files(self._image_echo2_loc, file_extension=self._file_extension)
            self._te2 = float(pydicom.read_file(self._image_echo2_loc / dicom_files[0]).EchoTime) * 1e-3
        else:
            self._te2 = te2

        self._delta_t = abs(te2 - te1)
        print(f'echo time 1: {self._te1: 1.2f} ms')
        print(f'echo time 2: {self._te2: 1.2f} ms')
