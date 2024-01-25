import time
from pathlib import Path
import numpy as np
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label
from scipy.interpolate import NearestNDInterpolator
from matplotlib import pyplot as plt
import logging
import os
import json
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import transform
from warnings import warn
from .utilities import _get_MR_acquisition_data
from .utilities import dicom_to_numpy

ch = logging.StreamHandler()
formatter = logging.Formatter('[%(filename)s: line %(lineno)d %(levelname)8s] %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(ch)
logger.setLevel(logging.INFO)
logger.propagate = False


class MarkerVolume:
    """
    This class excepts either a path to a folder of dicom files or a numpy array of marker positions.
    In the former case, it will attempt to automatically extract the marker positions. In the latter, it will just
    use the input marker positions. Once the instance is created, it contains a number of handy plots and functions.

    :param input_data_path: Either a path to a folder of dicoms, or a numpy array with three columns corresponding to the
        [x,y,z] marker positions
    :type input_data_path: path, string, np.ndarray
    :param ImExtension: change to whatever extension your dicom files, most commonly either dcm or IMA
    :type ImExtension: str, optional
    :param r_min: any markers at radii less than this will be discounted.
    :type r_min: float ,optional
    :param r_max: any markers at radii greater than this will be discounted.
    :type r_max: float, optional
    :param iterative_segmentation: If set to true, a slower iterative method will be used to find the threshold value
        used to segment the markers. Otherwise, otsu's method is used.
    :type iterative_segmentation: bool, optional
    :param threshold: Manually set the threshold value. If this is left as None otsu's method is used.
        This will replace the iterative segmentation method regardless of flag.
    :type threshold: float, optional
    :param n_markers_expected: if you know how many markers you expect, you can enter the value here. The code will
        then warn you if it finds a different number.
    :type n_markers_expected: int, optional
    :param gaussian_image_filter_sd: size of filter to blur input data with. Higher values are better at removing
        noise prior to thresholding, but may also blur out small/low signal markers!
    :type gaussian_image_filter_sd: float, optional
    :param correct_fat_water_shift: if True, will attempt to automatically correct marker positions for calculated
        fat water shift
    :type correct_fat_water_shift: bool, optional
    :param marker_size_lower_tol: lower acceptance window for a marker compared to the median size of markers.
        e.g. markermarker_size_tol=0.5 will result in markers being accepted within when
        (1-marker_size_lower_tol)*median_marker_size <= marker_size <= (1+markermarker_size_tol)*marker_size_upper_tol
    :type marker_size_lower_tol: float, optional
    :param marker_size_upper_tol: upper acceptance window for a marker compared to the median size of markers.
    :type marker_size_upper_tol: float, optional
    :param verbose: if True, prints a lot of info about the marker extraction process to the screen. mostly useful
        during debugging
    :type verbose: bool, optional
    :param fat_shift_direction: We aren't sure yet if we can confidently predict the direction that the fat water shift
        will be in. if you set correct_fat_water_shift=True, you can change the direction it shifts in here. it should
        take a value of -1, or 1. pull requests making this easier to use very welcome!
    :type fat_shift_direction: int, optional
    """

    def __init__(self, input_data, ImExtension='dcm', r_min=None, r_max=None, iterative_segmentation=False,
                 threshold=None, n_markers_expected=None, fat_shift_direction=None, verbose=False,
                 gaussian_image_filter_sd=1, correct_fat_water_shift=False, marker_size_lower_tol=0.9,
                 marker_size_upper_tol=1):

        self.verbose = verbose
        self._file_extension = ImExtension
        self._correct_fat_water_shift = correct_fat_water_shift
        self.dicom_data = None  # overwritten below if possible.
        self._n_markers_expected = n_markers_expected
        self._r_min = r_min
        self._r_max = r_max
        self._cutoffpoint = threshold
        self._iterative_segmentation = iterative_segmentation
        self._chemical_shift_vector = None
        self._gaussian_image_filter_sd = gaussian_image_filter_sd
        self._marker_size_lower_tol = marker_size_lower_tol
        self._marker_size_upper_tol = marker_size_upper_tol

        if isinstance(input_data, (Path, str)):
            self.input_data_path = Path(input_data)
            if self.input_data_path.is_dir():
                # dicom input
                self.InputVolume, self.dicom_affine, (self.X, self.Y, self.Z) = \
                    dicom_to_numpy(self.input_data_path, file_extension=self._file_extension, return_XYZ=True)
                self.dicom_data = _get_MR_acquisition_data(self.input_data_path, self.X, self.Y, self.Z,
                                                           file_extension=self._file_extension)
                # Segmenting markers
                self._filter_volume_by_r()
                start_time = time.perf_counter()
                self.ThresholdVolume, self.BlurredVolume = self._threshold_volume(self.InputVolume)
                centroids = self._find_contour_centroids()
                self.MarkerCentroids = pd.DataFrame(centroids, columns=['x', 'y', 'z'])
                end_time = time.perf_counter()
                print(f'total segmentation time: {end_time - start_time: 1.1f} s')

                # Correct for oil water shift
                if self._correct_fat_water_shift:
                    self._calculate_chemical_shift_vector(fat_shift_direction=fat_shift_direction)
                    self.MarkerCentroids = self.MarkerCentroids + self._chemical_shift_vector
            elif (os.path.isfile(self.input_data_path) and os.path.splitext(self.input_data_path)[1] == '.json'):
                # slicer input
                with open(self.input_data_path) as f:
                    data = json.load(f)
                arrays = [el for el in data['markups'][0]['controlPoints']]
                points = [el['position'] for el in arrays]
                self.MarkerCentroids = pd.DataFrame(points, columns=['x', 'y', 'z'])

                # look for dicom_data json file
                try:
                    with open(os.path.join(self.input_data_path.parent, 'dicom_data.json')) as f:
                        self.dicom_data = json.load(f)
                except:
                    logger.warning(
                        f'MR data file dicom_data.json not found at {self.input_data_path.parent}. Continuing')
            else:
                raise FileNotFoundError(f'could not find any data at {self.input_data_path}')
        elif isinstance(input_data, np.ndarray):
            # don't need to do anything then.
            self.MarkerCentroids = pd.DataFrame(input_data, columns=['x', 'y', 'z'])
            logger.warning("numpy input is deprecated, please use pandas input with ['x', 'y', 'z'] cols instead")
        elif isinstance(input_data, pd.core.frame.DataFrame):
            _expected_columns = ['x', 'y', 'z']
            if not all([col in input_data for col in _expected_columns]):
                raise AttributeError("input pandas data must have columns called ['x', 'y', 'z']")
            self.MarkerCentroids = input_data.copy()
        else:
            raise TypeError(f'Uknown data input: {type(input_data)}')

        # insert the radial value of each marker:
        self.MarkerCentroids['r'] = self.MarkerCentroids.apply(
            lambda row: np.sqrt(row[0] ** 2 + row[1] ** 2 + row[2] ** 2), axis=1)
        if self.MarkerCentroids.r.mean() < 10:
            logger.warning('it appears that your input data is in m, not mm - please use mm!')

        self._filter_markers_by_r()
        if self._n_markers_expected is not None:
            if not (self.MarkerCentroids.shape[0] == n_markers_expected):
                logger.warning(f'For data {self.input_data_path}\n'
                               f'You entered that you expected to find'
                               f' {n_markers_expected}, but actually found {self.MarkerCentroids.shape[0]}.')

    def _filter_volume_by_r(self):
        """
        Accelerates and improves segmenting by setting all voxels outside of specified radius to the background
        """
        if self._r_max or self._r_min:
            _max = 10000
            _min = -1
            if self._r_max:
                _max = self._r_max
            if self._r_min:
                _min = self._r_min
            background = np.median(self.InputVolume)

            it = np.nditer(self.InputVolume, flags=['multi_index'])
            for voxel in it:
                distance = np.sqrt(
                    self.X[it.multi_index] ** 2 + self.Y[it.multi_index] ** 2 + self.Z[it.multi_index] ** 2)
                if distance < _min or distance > _max:
                    self.InputVolume[it.multi_index] = background

    def _filter_markers_by_r(self):
        """
        Remove any markers that are less than r_min or more than r_max
        """
        if self._r_max:
            keep_ind = abs(self.MarkerCentroids.r) < self._r_max
            self.MarkerCentroids = self.MarkerCentroids[keep_ind].reset_index(drop=True)
        if self._r_min:
            keep_ind = abs(self.MarkerCentroids.r) > self._r_min
            self.MarkerCentroids = self.MarkerCentroids[keep_ind].reset_index(drop=True)

    def _calculate_chemical_shift_vector(self, fat_shift_direction=1):
        """
        Calculates a chemical shift vector (x,y,z) to be applied to each marker to account for the fat/water shift.

        :param fat_shift_direction: dictates whether the shift is positive or negative. we aren't sure how consisten
            this is between scanners so have left as used defined.
        """
        if self.dicom_data is None:
            logger.warning('No dicom data for chemical shift')
            return

        # TODO check these
        if fat_shift_direction == 'AP' or fat_shift_direction == 'RL' or fat_shift_direction == 'FH' or fat_shift_direction == -1:
            fat_shift_direction = -1
        elif fat_shift_direction == 'PA' or fat_shift_direction == 'LR' or fat_shift_direction == 'HF' or fat_shift_direction == 1:
            fat_shift_direction = 1
        else:
            fat_shift_direction = 0
            logger.warning('Chemical shift input was invalid')

        # Calculate magnitude (in pixels) of the chemical shift based on the shift frequency of a 1T MRI field (149Hz)
        directions = ['x', 'y', 'z']
        direction = np.array(([int(el == self.dicom_data['freq_encode_direction']) for el in directions]))
        # ^ this just gives a vector with a 1 in the correct direction e.g. [0,1,0]
        shift = direction * self.dicom_data['chem_shift_magnitude']
        self._chemical_shift_vector = shift * np.array(self.dicom_data['pixel_spacing'])[
            direction.astype(bool)] * fat_shift_direction

    def _find_iterative_cutoff(self, blurred_volume):
        """
        Replaces the otsu threshold with a slower, iterative method to find the best threshold for segmenting markers.
        A range of valid thresholds are found and the mean is returned as the best threshold.
        If no valid thresholds are found, the mean of candidate thresholds that will miss some markers is returned.
        """

        if self._n_markers_expected is None:
            logger.warning('Iterative segmentation method requires the expected number of markers to be input.')
            return None

        # finds a range of thresholds that give a number of segments near to the number of expected markers
        histogram_division = 100

        candidate_thresholds = []  # Thresholds that have fewer markers than expected, use when no valid thresholds
        candidate_n_points = []
        valid_thresholds = []  # Thresholds that should give the corrected number of markers

        # divide the range into segments and calculate how many volumes result from each segment
        background_value = np.median(np.round(blurred_volume))
        intensity_range = np.max(blurred_volume) - background_value
        histogram_segment = intensity_range / histogram_division

        # Testing each histogram segment
        for i in range(1, histogram_division - 1):

            n_segments = 0
            cutoff = background_value + i * histogram_segment
            threshold_volume = blurred_volume > cutoff
            labels = label(threshold_volume, background=0)
            unique_labels = np.unique(labels)[1:]  # first label is background so skip
            # Check number of segments falls within valid range:
            tol = 0.1
            if len(unique_labels) < self._n_markers_expected * (1 - tol) or len(
                    unique_labels) > self._n_markers_expected * (1 + tol):
                continue  # Too few or too many segments
            else:
                # This threshold is possibly valid
                candidate_thresholds.append(cutoff)
                # Remove small volumes and check if still valid
                for label_level in unique_labels:
                    RegionInd = labels == label_level
                    if np.count_nonzero(RegionInd) > 2:
                        n_segments += 1
                candidate_n_points.append(n_segments)
                if n_segments >= self._n_markers_expected:
                    valid_thresholds.append(cutoff)
                    if self.verbose:
                        print('Threshold: ' + str(format(cutoff, '.2f')) + ', ' + str(n_segments) + ' segments.')
                else:
                    if self.verbose:
                        print('Threshold: ' + str(format(cutoff, '.2f')) + ' is a candidate.')

            if len(unique_labels) < 10:
                # threshold is too high and we can stop
                break

        if valid_thresholds:
            print(f'Using thresold of {np.mean(valid_thresholds): 1.1f}')
            return np.mean(valid_thresholds)
        elif len(candidate_thresholds) > 0:
            ind = np.argmin(abs(np.subtract(self._n_markers_expected, candidate_n_points)))
            best_threshold = candidate_thresholds[ind]
            logger.warning(f'No valid thresholds were found. Using closest possible value of {best_threshold}')
            return best_threshold
        else:
            logger.warning('No valid thresholds were found. Try lowering gaussian blur level.')
            return None

    def _threshold_volume(self, VolumeToThreshold):
        """
        For a 3D numpy array, turn all elements < cutoff to zero, and all other elements to 1.
        If no cutoff is entered, otsu's method is used to auto-threshold.
        """

        BlurredVolume = gaussian(VolumeToThreshold, sigma=self._gaussian_image_filter_sd)

        if self._cutoffpoint is not None:
            # If cutoff point has been manually entered, go straight to thresholding
            self._iterative_segmentation = False
        if self._iterative_segmentation is True:
            self._cutoffpoint = self._find_iterative_cutoff(BlurredVolume)
        if self._cutoffpoint is None:
            self._cutoffpoint = threshold_otsu(BlurredVolume)

        ThresholdVolume = BlurredVolume > self._cutoffpoint

        return ThresholdVolume, BlurredVolume

    def _get_marker_max_min_volume(self):
        """
        figures out the range of voxels that a marker should have
        """
        self._n_voxels = []  # going to keep track of this so we can remove any very small regions if needed

        for label_level in self.unique_labels:  # first label is background so skip
            # extract x,y,z of each connected region
            RegionInd = self._labels == label_level
            self._n_voxels.append(np.count_nonzero(RegionInd))
        # Set up min and max marker volumes
        self._n_voxels_median = np.median(np.array(self._n_voxels))
        self._voxel_min = (1 - self._marker_size_lower_tol) * self._n_voxels_median
        self._voxel_max = (1 + self._marker_size_upper_tol) * self._n_voxels_median

        # Modify based on number of markers
        if self._n_markers_expected is not None:
            self._n_voxels_copy = np.array(self._n_voxels)
            voxel_difference = np.absolute(np.array(self._n_voxels) - self._n_voxels_median)
            # Loop that deletes the marker size with the largest difference from median until the expected # of markers
            while len(self._n_voxels_copy) > self._n_markers_expected:
                self._n_voxels_copy = np.delete(self._n_voxels_copy, np.argmax(voxel_difference))
                voxel_difference = np.delete(voxel_difference, np.argmax(voxel_difference))
            # Set min and max based on the remaining list
            if np.min(self._n_voxels_copy) > self._voxel_min:
                self._voxel_min = np.min(self._n_voxels_copy)
            if np.max(self._n_voxels_copy) < self._voxel_max:
                self._voxel_max = np.max(self._n_voxels_copy)

        # 3 voxels is the absolute floor
        if self._voxel_min < 3:
            self._voxel_min = 3

    def _remove_load(self):
        '''
        if it seems very likely that the only thing that has been segmented is the load
        remove it and try again. We may want to enable this in situations where any large object is detected...
        '''

        for label_level in self.unique_labels:
            RegionInd = self._labels == label_level
            self.InputVolume[RegionInd] = 0
        self._cutoffpoint = None
        self.ThresholdVolume, self.BlurredVolume = self._threshold_volume(self.InputVolume)
        self._labels = label(self.ThresholdVolume, background=0)
        self.unique_labels = np.unique(self._labels)[1:]  # first label is background so skip

    def _print_thresholding_information(self):
        print('Threshold value: ' + str(self._cutoffpoint))
        print('Median marker volume: ' + str(self._n_voxels_median))
        print('Expected marker volume: ' + str(self._voxel_min) + ' to ' + str(self._voxel_max))
        print(f'Regions to check: {self.unique_labels.shape}')

    def _print_thresholding_summary(self, n_found_markers, n_skipped_markers):
        print('----------')
        print(f'Total regions: {len(self.unique_labels)}')
        print(f'Total markers found:   {n_found_markers}')
        print(f'Total skipped:         {n_skipped_markers}')

    def _find_contour_centroids(self):
        """
        This code loops through all the found regions, extracts the cartesian coordiantes, and takes the
        intensity-weighted average as the centroid.
        It excludes volumes that are bigger/ smaller than the median using the params
        marker_size_upper_tol and marker_size_lower_tol
        """
        self._labels = label(self.ThresholdVolume, background=0)
        self.unique_labels = np.unique(self._labels)[1:]  # first label is background so skip
        if self.unique_labels.shape[0] < 3:
            logger.warning('automatic thresholding didnt work, trying to remove load and try again...')
            self._remove_load()
        self._get_marker_max_min_volume()

        if self.verbose:
            self._print_thresholding_information()

        x_centroids = []
        y_centroids = []
        z_centroids = []
        skipped = 0

        # extract x,y,z of each connected region
        for i, label_level in enumerate(self.unique_labels):
            RegionInd = np.equal(self._labels, label_level)
            voxels = np.count_nonzero(RegionInd)
            if voxels < self._voxel_min or voxels > self._voxel_max:
                skipped += 1
                if self.verbose:
                    print('\033[93mRegion ' + str(i + 1) + ': Skipped, v = ' + str(voxels) + '\033[0m')
                continue  # skip outliers

            region_sum = np.sum(self.InputVolume[RegionInd])
            weighted_x = np.sum(np.multiply(self.X[RegionInd], self.InputVolume[RegionInd]))
            weighted_y = np.sum(np.multiply(self.Y[RegionInd], self.InputVolume[RegionInd]))
            weighted_z = np.sum(np.multiply(self.Z[RegionInd], self.InputVolume[RegionInd]))
            x_centroids.append(weighted_x / region_sum)
            y_centroids.append(weighted_y / region_sum)
            z_centroids.append(weighted_z / region_sum)
            if self.verbose:
                print('Region ' + str(i + 1) + ': Marker, v = ' + str(voxels))

        if self.verbose:
            self._print_thresholding_summary(len(x_centroids), skipped)

        return np.array([x_centroids, y_centroids, z_centroids]).T

    # public methods

    def rotate_markers(self, xaxis_angle=0, yaxis_angle=0, zaxis_angle=0):
        """
        rotate markers using euler angles

        :param xaxis_angle: rotation around x angle in degrees
        :param yaxis_angle: rotation around y angle in degrees
        :param zaxis_angle: rotation around z angle in degrees
        """

        rotate = np.array([xaxis_angle, yaxis_angle, zaxis_angle])
        rotation_vector = transform.Rotation.from_euler('xyz', rotate, degrees=True)
        # Transform centroids
        rotated = rotation_vector.apply(self.MarkerCentroids[['x', 'y', 'z']])
        self.MarkerCentroids = pd.DataFrame(rotated, columns=['x', 'y', 'z'])
        self.MarkerCentroids['r'] = self.MarkerCentroids.apply(
            lambda row: np.sqrt(row[0] ** 2 + row[1] ** 2 + row[2] ** 2), axis=1)

    def translate_markers(self, x_shift=0, y_shift=0, z_shift=0):
        """
        translate markers

        :param x_shift: x translation in mm
        :param y_shift: y translation in mm
        :param z_shift: z translation in mm
        :return:
        """
        translate = np.array([x_shift, y_shift, z_shift])
        translated = self.MarkerCentroids[['x', 'y', 'z']] + translate
        self.MarkerCentroids = pd.DataFrame(translated, columns=['x', 'y', 'z'])
        self.MarkerCentroids['r'] = self.MarkerCentroids.apply(
            lambda row: np.sqrt(row[0] ** 2 + row[1] ** 2 + row[2] ** 2), axis=1)

    def plot_3D_markers(self, title='3D marker positions'):  # pragma: no cover
        """
        Just a quick plot of the marker positions; very useful sanity check!
        """

        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')
        axs.scatter(self.MarkerCentroids.x, self.MarkerCentroids.y, self.MarkerCentroids.z, s=10)
        axs.set_xlabel('X [mm]')
        axs.set_ylabel('Y [mm]')
        axs.set_zlabel('Z [mm]')
        axs.set_title(title)
        axs.set_box_aspect(
            (np.ptp(self.MarkerCentroids.x), np.ptp(self.MarkerCentroids.y), np.ptp(self.MarkerCentroids.z)))
        plt.show()

    def perturb_marker_positions(self, random_perturbation, systemic_perturbation=0):
        """
        add random noise to the marker positions. Useful to perform sensitivity analysis. operates in place.

        :param random_perturbation: upper and lower limits of random_perturbation to each marker coordinates. uniform
            distribution is used
        :type random_perturbation: float
        :param systemic_perturbation: constant offset added to all markers
        :type systemic_perturbation: float, optional
        """

        x_peturb = np.random.uniform(low=-random_perturbation, high=random_perturbation,
                                     size=self.MarkerCentroids.x.shape)
        y_peturb = np.random.uniform(low=-random_perturbation, high=random_perturbation,
                                     size=self.MarkerCentroids.x.shape)
        z_peturb = np.random.uniform(low=-random_perturbation, high=random_perturbation,
                                     size=self.MarkerCentroids.x.shape)
        self.MarkerCentroids.x = self.MarkerCentroids.x + x_peturb + systemic_perturbation
        self.MarkerCentroids.y = self.MarkerCentroids.y + y_peturb + systemic_perturbation
        self.MarkerCentroids.z = self.MarkerCentroids.z + z_peturb + systemic_perturbation

    def export_to_slicer(self, save_path=None, filename='slicer_centroids'):
        """
        export a json file that can be read in by slicer, allowing a good way to visualise marker segmentation
        performance. This file will be saved at the same spot as input data if dicom was input, otherwise it will be
        saved to save_path.
        """
        # sort out where to save:
        if save_path is None:
            try:
                _save_path = self.input_data_path
            except AttributeError:  # indicates numpy input
                logger.error('To save numpy data to slicer format, please enter a value for "save_path" when'
                             'calling MarkerVolume.export_to_slicer')
                return
        else:
            _save_path = Path(save_path)

        file_content = []  # will fill up one line at a time
        file_content.append('{"@schema": "https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/'
                            'Markups/Resources/Schema/markups-schema-v1.0.0.json#",\n')
        file_content.append('"markups": [{"type": "Fiducial", "coordinateSystem": "LPS", "controlPoints": [\n')
        for index, row in self.MarkerCentroids.iterrows():
            if index == self.MarkerCentroids.shape[0] - 1:
                # for last entry have to delete a comma
                file_content.append(
                    f'{{ "label": " ", "position": [{row.x: 1.2f}, {row.y: 1.2f}, {row.z: 1.2f}], "locked": true}}\n')
            else:
                file_content.append(
                    f'{{ "label": " ", "position": [{row.x: 1.2f}, {row.y: 1.2f}, {row.z: 1.2f}], "locked": true}},\n')

        file_content.append(']}]}')
        # write contents
        filename = filename + '.mrk.json'
        file_loc = _save_path / filename
        with open(file_loc, 'w') as f:
            f.writelines(file_content)

    def save_dicom_data(self, save_path=None, filename='dicom_data'):
        """
        save the dicom data as json.  This file will be saved at the same
        spot as input data if dicom was input, otherwise it will be saved
        to save_path.
        """

        if self.dicom_data is None:
            logger.warning('cannot save dicom data because there is none...')
            return
        if save_path is None:
            save_path = self.input_data_path
        save_path = Path(save_path)
        _file_name, _file_extension = os.path.splitext(filename)
        if _file_extension == '':
            _file_extension = '.json'
        filename = _file_name + _file_extension
        full_filename = save_path / filename

        with open(full_filename, 'w') as f:
            json.dump(self.dicom_data, f)


class MatchedMarkerVolumes:
    """
    :param GroundTruthData:
    :type GroundTruthData: MarkerVolume object, or path to dicom folder, or numpy array
    :param DistortedData:
    :type DistortedData: MarkerVolume object, or path to dicom folder, or numpy array
    :param reverse_gradient_data:
    :type reverse_gradient_data: None, or MarkerVolume object, or path to dicom folder, or numpy array
    :param WarpSearchData: if True, position of found markers is used to update expected positions for ground truth.
        Recomended is True if there is substantial distortion present.
    :type WarpSearchData: boolean
    :param AutomatchMarkers: if True, will automatically match the distorted data to the ground truth using a search
        algorith. If False, will simply use input data in the order it is entered. This is useful if you have sets
        of markers that were independently matched through some other process.
    :type AutomatchMarkers: boolean
    :param AllowDoubleMatching: when False, each ground truth marker is only allowed to match with one distorted
        marker. However, sometimes this can result in a chain reaction where one marker gets mismatched then every
        other marker gets mismatched too. In these cases you can set this parameter to True, and any double matched
        markers will simply be deleted.
    :type AllowDoubleMatching: bool, optional
    :param sorting_method: 'radial' or 'nearest'. This is only important if WarpSearchData is True; in that case,
        we are building a motion model on-the-fly based on previously seen markers, so the order we see them can be
        important. For data on a sphere, you may have more success by setting this to 'nearest', but this is a bit
        of an untested feature
    :type sorting_method: str, optional
    :param n_refernce_markers: the n inner_most Reference markers are used to align the ground truth to the MR volume
        before  matching. this can correct for setup errors. Note that we always move the reference. the best way
        for this to not matter either way is to not have any set up error!
    :type n_refernce_markers: int, optional
    :param skip_unmatchable_markers: if true, distorted markers where the nearest centroid marker is above the maximum
        distortion threshold are removed
    :type skip_unmatchable_markers: bool, optional
    """

    def __init__(self, GroundTruthData, DistortedData, reverse_gradient_data=None, WarpSearchData=True,
                 AutomatchMarkers=True, AllowDoubleMatching=False, sorting_method='radial', n_refernce_markers=0,
                 skip_unmatchable_markers=True):

        # warping parameters:
        self.WarpSearchData = WarpSearchData
        self.AutomatchMarkers = AutomatchMarkers
        self._motion_estimate_update_rate = 1  # 1 = every found marker, which is potentially more accurate but slower.
        self.sorting_method = sorting_method
        self.AllowDoubleMatching = AllowDoubleMatching

        # match check parameters:
        self._mean_match_tolerance = 20  # warning raised if mean larger than this
        self._max_match_tolerance = 30  # warning raised if max larger than this

        self._skip_unmatchable_markers = skip_unmatchable_markers
        self._n_reference_markers = n_refernce_markers

        # marker data:
        self.ground_truth_centroids = GroundTruthData.MarkerCentroids
        self.distorted_centroids = DistortedData.MarkerCentroids

        if reverse_gradient_data is None:
            self.distorted_centroidsRev = None
        else:
            self.distorted_centroidsRev = reverse_gradient_data.MarkerCentroids


        # run analysis:
        if self._n_reference_markers > 0:
            self._align_reference()
        if self._skip_unmatchable_markers:
            self._remove_unmatchable_distorted_centroids()
        self._check_input_data()



        if self.AutomatchMarkers:
            self.distorted_centroids = self._sort_distorted_centroids(self.distorted_centroids)
            self._CentroidMatch = self._match_distorted_markers_to_ground_truth(self.distorted_centroids)
            self._check_match()
            if self.distorted_centroidsRev is not None:
                self.distorted_centroidsRev = self._sort_distorted_centroids(self.distorted_centroidsRev)
                self._CentroidMatchRev = self._match_distorted_markers_to_ground_truth(self.distorted_centroidsRev)
        else:
            self._create_marker_data_prematched_markers()
        self._handle_double_matched_markers()

        # analyse marker positions
        self._generate_marker_position_data()

    def _check_input_data(self):
        """
        put any tests on the input data here to make sure the rest of the code can run.

        - tests that there at least as many ground truth centroids than distorted centroids
        - tests that data appears to be in mm (requirement)
        """
        if not self.ground_truth_centroids.shape[0] >= self.distorted_centroids.shape[0]:
            raise ValueError('There are fewer ground truth centroids than distorted centroids; this means the disorted'
                             'centroids can never be matched')

        if not self.AutomatchMarkers:
            if not self.ground_truth_centroids.shape[0] == self.distorted_centroids.shape[0]:
                raise ValueError(
                    f'Since you have set AutomatchMarkers=False, there must be equal numbers of ground truth'
                    f'markers and distorted markers, but ground truth has {self.ground_truth_centroids.shape[0]}'
                    f'and distorted has {self.distorted_centroids.shape[0]}.')

        if self.distorted_centroidsRev is not None:
            if not self.ground_truth_centroids.shape[0] >= self.distorted_centroidsRev.shape[0]:
                raise ValueError('There are fewer ground truth centroids than reversed gradient distorted centroids;'
                                 ' this means the disorted centroids can never be matched. Quitting')

        if np.max([self.ground_truth_centroids.x.max(), self.ground_truth_centroids.y.max(),
                   self.ground_truth_centroids.z.max()]) < 50:
            logger.warning('it appears that the ground truth centroids may not be defined in mm?'
                           'this will cause the algorithm to fail. Continuing but double check this!')

        if np.max([self.distorted_centroids.x.max(), self.distorted_centroids.y.max(),
                   self.distorted_centroids.z.max()]) < 50:
            logger.warning('it appears that the distorted centroids may not be defined in mm?'
                           'this will cause the algorithm to fail. Continuing but double check this!')
        if self.distorted_centroidsRev is not None:
            if np.max([self.distorted_centroidsRev.x.max(), self.distorted_centroidsRev.y.max(),
                       self.distorted_centroidsRev.z.max()]) < 50:
                logger.warning('it appears that the reverse gradient distorted centroids may not be defined in mm? '
                               'this will cause the algorithm to fail. Continuing but double check this!')

    def _get_grid_spacing(self):
        """
        Get the spacing between markers using the ground truth data. This is used to build the kernel for RBF
        interpolation during the marker matching step.
        """
        self._x_grid_spacing = np.mean(np.diff(np.sort(np.unique(self.ground_truth_centroids.x))))
        self._y_grid_spacing = np.mean(np.diff(np.sort(np.unique(self.ground_truth_centroids.y))))
        self._z_grid_spacing = np.mean(np.diff(np.sort(np.unique(self.ground_truth_centroids.z))))

    def _generate_extrapolant_models(self, CentroidsToSearch, distorted_centroids):
        """
        Use the found marker positions to extrapolate the distortion at other positions.

        :param CentroidsToSearch: Current search data. At this point this is based on the
            ground truth data with matched markers removed.

        returns X_dist, Y_dist, Z_dist, which is the distortion at every marker position in
        CentroidsToSearch.
        """
        points = self.ground_truth_centroids.iloc[self._MatchInd][['x', 'y', 'z']]
        distortion = self.ground_truth_centroids.iloc[self._MatchInd].to_numpy() - \
                     distorted_centroids.iloc[0:len(self._MatchInd)].to_numpy()

        # Nearest Neighbore interpolator
        self._X_interpolator = NearestNDInterpolator(points, distortion[:, 0])
        self._Y_interpolator = NearestNDInterpolator(points, distortion[:, 1])
        self._Z_interpolator = NearestNDInterpolator(points, distortion[:, 2])

        points = CentroidsToSearch[['x', 'y', 'z']].to_numpy().squeeze()
        X_dist = self._X_interpolator(points)
        Y_dist = self._Y_interpolator(points)
        Z_dist = self._Z_interpolator(points)
        return X_dist, Y_dist, Z_dist

    def _align_reference(self, rotation=True):
        """
        Moves the ground truth markers to be aligned with the distorted markers based on the central n refernce markers
        """
        if self.ground_truth_centroids is None or self.distorted_centroids is None:
            logger.warning('No centroids loaded')
            return

        def _calculate_radial_distance(centroids):
            """
            insert the radial value of each marker:
            """
            centroids['r'] = centroids.apply(
                lambda row: np.sqrt(row[0] ** 2 + row[1] ** 2 + row[2] ** 2), axis=1)
            if centroids.r.mean() < 10:
                logger.warning('it appears that your input data is in m, not mm - please use mm!')
            return centroids

        def _match_crosshair(search_centroids, reference_centroids):
            """
            match the distorted/ ground truth
            """

            _matched_reference = pd.DataFrame(columns=['x', 'y', 'z'])
            search_centroids.reset_index(inplace=True, drop=True)
            reference_centroids.reset_index(inplace=True, drop=True)
            for index, row in reference_centroids.iterrows():
                # Calculate distance from this marker to each item in search list
                Distances = cdist(search_centroids, np.atleast_2d(row), metric='euclidean')
                # Match marker with closest and save it
                ind_nearest = np.argmin(Distances)
                matched = search_centroids.iloc[ind_nearest]
                _matched_reference.loc[index] = matched

            return _matched_reference

        # Find the central position in the ground truth data
        ground_truth = self.ground_truth_centroids.copy()
        central_x = (ground_truth.x.min() + ground_truth.x.max()) / 2
        central_y = (ground_truth.y.min() + ground_truth.y.max()) / 2
        central_z = (ground_truth.z.min() + ground_truth.z.max()) / 2

        # Calculate the distance between each marker and central position
        ground_truth['r_centre'] = ground_truth.apply(
            lambda row: np.sqrt((row[0] - central_x) ** 2 + (row[1] - central_y) ** 2 + (row[2] - central_z) ** 2),
            axis=1)

        # Find the crosshair reference markers
        self._gt_reference = ground_truth.nsmallest(self._n_reference_markers, 'r_centre')[['x', 'y', 'z']]
        self._dist_reference = self.distorted_centroids.nsmallest(self._n_reference_markers, 'r')[['x', 'y', 'z']]

        # Calculate translation vector
        translation_vector = self._dist_reference.mean() - self._gt_reference.mean()

        # Calculate rotation vector
        if self._n_reference_markers < 3:
            rotation = False
        else:
            matched_reference = _match_crosshair(self._gt_reference + translation_vector, self._dist_reference)
            rotation_vector, error = transform.Rotation.align_vectors(self._dist_reference, matched_reference)

        # Transform centroids
        aligned = self.ground_truth_centroids[['x', 'y', 'z']] + translation_vector
        print(f'aligning ground truth centroids by \n{translation_vector}')
        if rotation is True:
            aligned = rotation_vector.apply(aligned)
            print(f'rotating ground truth centroids by \n{rotation_vector.as_matrix()}')


        aligned = pd.DataFrame(aligned, columns=['x', 'y', 'z'])
        self.ground_truth_centroids = _calculate_radial_distance(aligned)

    def _remove_unmatchable_distorted_centroids(self):
        """
            Calculates the distance between each distorted marker and each ground truth markers
            Distorted markers where the nearest ground truth marker exceeds the distance thresholds are deleted from
             the search list
        """

        unmatchable_index = []
        CentroidsToSearch = self.ground_truth_centroids.copy()
        DistortedCentroids = self.distorted_centroids.copy()

        for index, row in DistortedCentroids.iterrows():
            XA = CentroidsToSearch[["x", "y", "z"]].to_numpy()
            XB = row[0:3].values
            # Calculate distance from this marker to each item in search list
            Distances = cdist(XA, np.atleast_2d(XB), metric='euclidean')
            # Match marker with closest and save it
            distance_to_closest = np.min(Distances)
            if distance_to_closest > self._max_match_tolerance:
                unmatchable_index.append(index)
        if unmatchable_index:
            logger.warning(f'removing {len(unmatchable_index)} markers because no ground truth markers were within'
                           f'{self._max_match_tolerance} mm')

        self.distorted_centroids = self.distorted_centroids.drop(unmatchable_index, axis=0)

    def _match_distorted_markers_to_ground_truth(self, CentroidsToMatch):
        """
        find the closest ground truth marker to a given distortion. The logic of this algorithm is based on Paul
        Liu's matlab version

        1. Distorted data is sorted by r so that the smallest (least distorted) markers are processed first
        2. the match is simply based on the minimum distance in cartesian space between the distorted marker and the
            search markers.
        3. Once a match has been found, that ground truth is removed from the search space
        4. if self.WarpSearchData=True, the search markers are warped according to the previously found data. For cases
            where distortion is small this is probably unnecessary, but for extreme distortions this algorithm will fail
            without this
        """

        self._MatchInd = []
        distances = []
        CentroidsToSearch = self.ground_truth_centroids.copy()
        for index, row in CentroidsToMatch.iterrows():

            XA = CentroidsToSearch[["x", "y", "z"]].to_numpy()
            XB = row[0:3].values
            # Calculate distance from this marker to each item in search list
            Distances = cdist(XA, np.atleast_2d(XB), metric='euclidean')
            # Match marker with closest and save it
            ind_closest_point_in_ground_truth = np.argmin(Distances)
            self._MatchInd.append(ind_closest_point_in_ground_truth)
            distances.append(np.linalg.norm(
                self.ground_truth_centroids.iloc[ind_closest_point_in_ground_truth][['x', 'y', 'z']] - row[
                    ['x', 'y', 'z']]))
            if not self.AllowDoubleMatching:
                CentroidsToSearch.iloc[ind_closest_point_in_ground_truth] = pd.Series(
                    {'x': 1000, 'y': 1000, 'z': 1000, 'r': 1732})

            # warp search data?
            if self.WarpSearchData and (index % self._motion_estimate_update_rate == 0) and index > 0:
                # get extrapolated motion at all points:
                X_dist, Y_dist, Z_dist = self._generate_extrapolant_models(CentroidsToSearch, CentroidsToMatch)
                # regenerate centroids to search:
                CentroidsToSearch = self.ground_truth_centroids.copy()
                # remove found entries:
                if not self.AllowDoubleMatching:
                    CentroidsToSearch.iloc[self._MatchInd] = pd.Series({'x': 1000, 'y': 1000, 'z': 1000, 'r': 1732})
                # update the data with the extrapolated motion
                CentroidsToSearch[['x']] = CentroidsToSearch[['x']] - X_dist[:, None]
                CentroidsToSearch[['y']] = CentroidsToSearch[['y']] - Y_dist[:, None]
                CentroidsToSearch[['z']] = CentroidsToSearch[['z']] - Z_dist[:, None]

        # create a new data frame that has the distorted and ground truth markers
        MatchedMarkers = self.ground_truth_centroids.iloc[self._MatchInd]
        MatchedMarkers = MatchedMarkers.rename(columns={"x": "x_gt", "y": "y_gt", "z": "z_gt", "r": "r_gt"})
        MatchedMarkers.reset_index(inplace=True, drop=True)
        MatchedCentroids = pd.concat([CentroidsToMatch, MatchedMarkers], axis=1)
        MatchedCentroids['match_distance'] = distances

        return MatchedCentroids

    def _check_match(self):
        """
        just a simple test to detect where bad matches have been made
        could add more sophisticated things in here if needed
        :return:
        """

        _throw_warning = False
        if self._CentroidMatch.match_distance.mean() > self._mean_match_tolerance:
            _throw_warning = True
        if self._CentroidMatch.match_distance.max() > self._max_match_tolerance:
            _throw_warning = True

        if _throw_warning:
            warn(
                f'\n\nThe marker match may have failed.\n\nThe mean detected distortion is {self._CentroidMatch.match_distance.mean(): 1.1f} mm '
                f'and the max is {self._CentroidMatch.match_distance.max(): 1.1f}.'
                f'\nYou can continue by pressing any key, but you should visualize the data using the plot_3D_markers'
                f' method. If the match has failed, a simple remedy can be to trim the input data at the MarkerVolume'
                f'\n stage...')
            input("Press any key to continue...")

    def _handle_double_matched_markers(self):
        """
        Checks if any distorted markers have been matched to the same ground truth point
        """
        non_unique_ind = self._CentroidMatch.duplicated(subset=['x_gt', 'y_gt', 'z_gt'], keep=False)
        if any(non_unique_ind):
            logger.error(f'{np.count_nonzero(non_unique_ind)} distorted centroids were matched to the same ground truth'
                         f'marker. This may be indicative of a serious issue and you should proceed with caution. For '
                         f'now, handling this by removing these entries from the analysis pool and continuing...')
            self._CentroidMatch.drop(self._CentroidMatch.index[non_unique_ind])
            self._CentroidMatch.reset_index(drop=True)

        if self.distorted_centroidsRev is not None:
            non_unique_ind = self._CentroidMatchRev.duplicated(subset=['x_gt', 'y_gt', 'z_gt'], keep=False)
            if any(non_unique_ind):
                logger.error(
                    f'{np.count_nonzero(non_unique_ind)} distorted centroids were matched to the same ground truth'
                    f'marker. This may be indicative of a serious issue and you should proceed with caution. For '
                    f'now, handling this by removing these entries from the analysis pool and continuing...')
                self._CentroidMatchRev.drop(self._CentroidMatchRev.index[non_unique_ind])
                self._CentroidMatchRev.reset_index(drop=True)

    def _create_marker_data_prematched_markers(self):
        """
        If no autosorting is requested, this creates the MatchedCentroids DataFrame directly from the input data
        with no sorting etc.
        """
        self._CentroidMatch = self.ground_truth_centroids
        self._CentroidMatch = self._CentroidMatch.rename(columns={"x": "x_gt", "y": "y_gt", "z": "z_gt", "r": "r_gt"})
        self._CentroidMatch = pd.concat([self._CentroidMatch, self.distorted_centroids], axis=1)

        if self.distorted_centroidsRev is not None:
            self._CentroidMatchRev = self.ground_truth_centroids
            self._CentroidMatchRev = self._CentroidMatchRev.rename(
                columns={"x": "x_gt", "y": "y_gt", "z": "z_gt", "r": "r_gt"})
            self._CentroidMatchRev = pd.concat([self._CentroidMatchRev, self.distorted_centroidsRev], axis=1)

    def _generate_marker_position_data(self):
        """
        Convert the marker position data to the data we want to use for calculating gradient and B0 distortion
        If only one dataset has been entered, then this is trivial. However, if we have a forward and reverse dataset
        then we need to use these two seperate the gradient and B0 effects.
        """
        if self.distorted_centroidsRev is None:
            # then our life is easy
            self.MatchedCentroids = self._CentroidMatch
            self.MatchedCentroids = self.MatchedCentroids.rename(
                columns={'x': 'x_gnl', 'y': 'y_gnl', 'z': 'z_gnl', 'r': 'r_gnl'})
        else:
            # then we have to do some work to seperate G and B0

            # first we will want to correlate the marker positions by sorting by ground truth r
            self._CentroidMatch = self._CentroidMatch.sort_values('r_gt').reset_index(drop=True)
            self._CentroidMatchRev = self._CentroidMatchRev.sort_values('r_gt').reset_index(drop=True)
            # now find an index where the ground truth markers match; we will discard other markers
            rev_match_index = []
            forward_match_index = []
            tol = 1e-1  # just a small number not a sensitve parameter
            n_unmatched_markers = 0
            for index, row in self._CentroidMatch.iterrows():
                # find the closest ground truth point in the reversed data
                XA = self._CentroidMatchRev[["x_gt", "y_gt", "z_gt"]].to_numpy()
                XB = row[['x_gt', 'y_gt', 'z_gt']].values
                # Calculate distance from this marker to each item in search list
                Distances = cdist(XA, np.atleast_2d(XB), metric='euclidean')
                # calculate location of minimum
                min_dist_ind = np.argmin(Distances)
                # if it is less than tolerance, record the position
                if Distances[min_dist_ind] < tol:
                    rev_match_index.append(min_dist_ind)
                    forward_match_index.append(index)
                else:
                    n_unmatched_markers += 1  # don't actually need this as can infer from data sizes.

            if n_unmatched_markers > 0:
                logger.warning(
                    f'failed to match {n_unmatched_markers} of a total {self.ground_truth_centroids.shape[0]}'
                    f' input markers. These will be deleted from returned data')
            # now we have indexes linking our two datasets based on their ground truth match
            self._CentroidMatchRev = self._CentroidMatchRev.iloc[rev_match_index].reset_index()
            self._CentroidMatch = self._CentroidMatch.iloc[forward_match_index].reset_index()

            # marker position distortion from gradients is difference between forward and reversed:
            self.MatchedCentroids = (self._CentroidMatch[['x', 'y', 'z', 'r']] + self._CentroidMatchRev[
                ['x', 'y', 'z', 'r']]) / 2
            self.MatchedCentroids = self.MatchedCentroids.rename(
                columns={'x': 'x_gnl', 'y': 'y_gnl', 'z': 'z_gnl', 'r': 'r_gnl'})
            # add ground truth to data frame
            self.MatchedCentroids[['x_gt', 'y_gt', 'z_gt', 'r_gt']] = self._CentroidMatch[
                ['x_gt', 'y_gt', 'z_gt', 'r_gt']]
            # marker position distortion from B0 is difference between from total and gradient distortion:
            '''
            deta_z_B0=z_dis_AP-z_GNL; %% maximum distortion is 16mm
            deta_y_B0=y_dis_AP-y_GNL;  %% maximum distortio is 4mm  
            '''
            B0_dis = self.MatchedCentroids[['x_gnl', 'y_gnl', 'z_gnl']].to_numpy() - \
                     self._CentroidMatch[['x', 'y', 'z']].to_numpy()
            B0_pos = B0_dis + self.MatchedCentroids[['x_gt', 'y_gt', 'z_gt']].to_numpy()
            self.MatchedCentroids[['x_B0', 'y_B0', 'z_B0']] = B0_dis

    def _sort_distorted_centroids(self, data_to_sort):
        """
        sorts the data frame.
        sorting is important, because the search for a matching ground truth marker is carried out in order
        """

        if self.sorting_method == 'radial':
            # sort input data by r
            data_to_sort = data_to_sort.sort_values('r')
            data_to_sort.reset_index(inplace=True, drop=True)
        if self.sorting_method == 'closest':
            # find the distorted marker that has the closest match to a ground truth marker.
            min_distance = []
            XA = self.ground_truth_centroids[["x", "y", "z"]].to_numpy()
            for index, row in data_to_sort.iterrows():
                XB = row[0:3].values
                Distances = cdist(XA, np.atleast_2d(XB), metric='euclidean')
                min_distance.append(np.min(Distances))
            starting_marker_ind = np.argmin(min_distance)
            # now calculate the distances of every other distorted marker to this marker:
            XA = data_to_sort[["x", "y", "z"]].to_numpy()
            XB = np.atleast_2d(data_to_sort.loc[starting_marker_ind][["x", "y", "z"]].to_numpy())
            Distances = cdist(XA, XB)
            distance_to_start_ind = np.argsort(Distances, axis=0)
            assert distance_to_start_ind[0] == starting_marker_ind
            data_to_sort = data_to_sort.reindex(np.squeeze(distance_to_start_ind))
            data_to_sort.reset_index(inplace=True, drop=True)

        return data_to_sort

    # public methods

    def plot_3D_markers(self, add_arrows=True, title='3D marker positions'):  # pragma: no cover
        """
        Works very similarly to the MarkerVolume version, but plots both sets of markers and adds arrows

        :param add_arrows: if True, arrows drawn between matches
        :type add_arrows: bool, optional
        :param title: plot title
        :type title: str, optional
        """

        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')
        axs.scatter(self.MatchedCentroids.x_gt, self.MatchedCentroids.y_gt, self.MatchedCentroids.z_gt)
        axs.scatter(self.MatchedCentroids.x_gnl, self.MatchedCentroids.y_gnl, self.MatchedCentroids.z_gnl)
        if add_arrows:
            mag_color = np.sqrt((self.MatchedCentroids.x_gt - self.MatchedCentroids.x_gnl) ** 2
                                + (self.MatchedCentroids.y_gt - self.MatchedCentroids.y_gnl) ** 2
                                + (self.MatchedCentroids.z_gt - self.MatchedCentroids.z_gnl) ** 2)
            axs.quiver(self.MatchedCentroids.x_gnl, self.MatchedCentroids.y_gnl, self.MatchedCentroids.z_gnl,
                       self.MatchedCentroids.x_gt - self.MatchedCentroids.x_gnl,
                       self.MatchedCentroids.y_gt - self.MatchedCentroids.y_gnl,
                       self.MatchedCentroids.z_gt - self.MatchedCentroids.z_gnl,
                       edgecolors=("black",))

        axs.set_xlabel('X [mm]')
        axs.set_ylabel('Y [mm]')
        axs.set_zlabel('Z [mm]')
        axs.set_title(title)
        axs.set_box_aspect((np.ptp(self.MatchedCentroids.x_gt), np.ptp(self.MatchedCentroids.y_gt),
                            np.ptp(self.MatchedCentroids.z_gt)))
        plt.legend(['ground truth', 'distorted'])
        plt.show()

    def plot_compressed_markers(self, z_max=20, z_min=-20, add_arrows=True, title=None):  # pragma: no cover
        """
        compresses the 3D markers in the z plane, allowing a 2D visualisation.

        :param z_max: maximum z coordiante of markers to includes
        :type z_max: float, optional
        :param z_min: minimum z coordiante of markers to includes
        :type z_min: float, optional
        :param add_arrows: if True, arrows drawn between matches
        :type add_arrows: bool, optional
        :param title: plot title
        :type title: str, optional

        """

        def get_markers_as_function_of_z():
            """
            finds an index of markers within tolerence of z_pos
            """
            data_ind = np.logical_and(self.MatchedCentroids.z_gt > z_min, self.MatchedCentroids.z_gt < z_max)
            z_loc_data = self.MatchedCentroids[data_ind]
            return z_loc_data

        def plot_markers_inner(plot_data):
            """
            clears exising axes and plots plot_data
            """
            axs.clear()
            z_offset = (plot_data.z_gt - plot_data.z_gnl)
            axs.scatter(plot_data.x_gnl, plot_data.y_gnl, c='darkturquoise')  # c=z_offset,cmap="magma"
            axs.scatter(plot_data.x_gt, plot_data.y_gt, c='darkgoldenrod')
            # axs.legend(['Distorted', 'Ground Truth'])
            if add_arrows:
                axs.quiver(plot_data.x_gnl, plot_data.y_gnl,
                           plot_data.x_gt - plot_data.x_gnl,
                           plot_data.y_gt - plot_data.y_gnl,
                           angles="xy", units='xy', scale=1, width=1)
            axs.set_xlim([-150, 150])
            axs.set_ylim([-150, 150])
            axs.set_xlabel('x [mm]', fontsize=15)
            axs.set_ylabel('y [mm]', fontsize=15)
            axs.tick_params(axis='both', which='major', labelsize=12)
            axs.tick_params(axis='both', which='minor', labelsize=12)
            axs.axis("equal")
            axs.grid()
            plt.subplots_adjust(bottom=0.25)
            plt.legend(['distorted', 'ground truth'])
            if title:
                plt.title(title)

        fig, axs = plt.subplots(figsize=[8, 8], ncols=1, nrows=1)

        plot_data = get_markers_as_function_of_z()
        plot_markers_inner(plot_data)
        plt.show()

    def plot_reference_markers(self, title='Reference alignment markers'): # pragma: no cover

        try:
            fig = plt.figure()
            axs = fig.add_subplot(111, projection='3d')
            axs.scatter(self._gt_reference.x, self._gt_reference.y, self._gt_reference.z)
            axs.scatter(self._dist_reference.x, self._dist_reference.y, self._dist_reference.z)

        except AttributeError:
            logger.warning('Cannot plt reference  as it does not exist...')
            return

        axs.set_xlabel('X [mm]')
        axs.set_ylabel('Y [mm]')
        axs.set_zlabel('Z [mm]')
        axs.set_title(title)
        axs.set_box_aspect((np.ptp(self.MatchedCentroids.x_gt), np.ptp(self.MatchedCentroids.y_gt),
                            np.ptp(self.MatchedCentroids.z_gt)))
        plt.legend(['ground truth', 'distorted'])
        plt.show()

    def plot_3D_markers_with_color_scale(self, cmap=None, title=None, elevation=None, azimuth=None,
                                         vmin=None, vmax=None, return_axs=False): # pragma: no cover

        if cmap is None:
            cmap = plt.cm.viridis

        fig = plt.figure()
        axs = fig.add_subplot(111, projection='3d')

        x_dis = np.abs(self.MatchedCentroids.x_gt - self.MatchedCentroids.x_gnl)
        y_dis = np.abs(self.MatchedCentroids.y_gt - self.MatchedCentroids.y_gnl)
        z_dis = np.abs(self.MatchedCentroids.z_gt - self.MatchedCentroids.z_gnl)
        abs_dis = np.sqrt(x_dis ** 2 + y_dis ** 2 + z_dis ** 2)
        p = axs.scatter(self.MatchedCentroids.x_gt, self.MatchedCentroids.y_gt, self.MatchedCentroids.z_gt,
                        c=abs_dis, cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(p, ax=axs)
        if elevation:
            axs.elev = elevation
        if azimuth:
            axs.azim = azimuth
        axs.set_xlabel('X [mm]')
        axs.set_ylabel('Y [mm]')
        axs.set_zlabel('Z [mm]')
        if title:
            axs.set_title(title)
        axs.set_box_aspect((np.ptp(self.MatchedCentroids.x_gt), np.ptp(self.MatchedCentroids.y_gt),
                            np.ptp(self.MatchedCentroids.z_gt)))

        if return_axs:
            return axs
        else:
            plt.show()


    def report(self):
        x_dis = np.abs(self.MatchedCentroids.x_gt - self.MatchedCentroids.x_gnl)
        y_dis = np.abs(self.MatchedCentroids.y_gt - self.MatchedCentroids.y_gnl)
        z_dis = np.abs(self.MatchedCentroids.z_gt - self.MatchedCentroids.z_gnl)

        print(f'mean distortion: {self.MatchedCentroids.match_distance.mean(): 1.1f} mm, '
              f'std: {np.std(self.MatchedCentroids.match_distance): 1.1f} mm, '
              f'Max: {self.MatchedCentroids.match_distance.max(): 1.1f} mm'
              f'\nx: {np.mean(x_dis): 1.1f} +- {np.std(x_dis): 1.1f}. max: {np.max(x_dis): 1.1f}'
              f'\ny: {np.mean(y_dis): 1.1f} +- {np.std(y_dis): 1.1f}. max: {np.max(y_dis): 1.1f}'
              f'\nz: {np.mean(z_dis): 1.1f} +- {np.std(z_dis): 1.1f}. max: {np.max(z_dis): 1.1f}')
