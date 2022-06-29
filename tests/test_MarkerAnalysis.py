import os.path
import sys
from pathlib import Path
import numpy as np
import pandas as pd


this_dir = Path(__file__).parent
sys.path.insert(0, str(this_dir.parent))
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from MRI_DistortionQA.utilities import get_dicom_data


def test_mr_dicom_data_read_in():
    """
    - test that we read in dicom data and return a number at least in the vicinity of the right answer
    - tests that the dicom data exists and frequency encoding direction is correctly extracted
    - tests that we can export this data to slicer format, read it back in, and that those two datasets are correctly
        matched
    - tests that we can export and read back in the dicom_data
    """

    volume = MarkerVolume((this_dir / 'test_data' / 'MR_dicom').resolve())
    assert volume.MarkerCentroids.shape[0] > 50
    assert volume.MarkerCentroids.shape[0] < 65
    # there should be 58 I think so this is a very loose test but should catch gross errors

    assert volume.dicom_data is not None
    # test manual cutoff works:
    volume2 = MarkerVolume((this_dir / 'test_data' / 'MR_dicom').resolve(), cutoff_point=volume._cutoffpoint)
    assert volume.MarkerCentroids.shape == volume2.MarkerCentroids.shape
    assert volume.dicom_data['freq_encode_direction'] == 'x'

    volume.export_to_slicer(this_dir / 'test_data', filename='dicom_export_test')
    volume.save_dicom_data(this_dir / 'test_data', filename='dicom_data_export_test')
    new_dicom_data = get_dicom_data(this_dir / 'test_data' / 'dicom_data_export_test.json')
    assert new_dicom_data == volume.dicom_data
    volume2 = MarkerVolume(str(this_dir / 'test_data' / 'dicom_export_test.mrk.json'))
    assert np.allclose(volume.MarkerCentroids, volume2.MarkerCentroids, rtol=1e-01, atol=1e-01)


def test_ct_dicom_data_read_in():
    """
    - test that we read in dicom data and return a number at least in the vicinity of the right answer
    - tests that the dicom data exists and frequency encoding direction is correctly extracted
    - tests that we can export this data to slicer format, read it back in, and that those two datasets are correctly
        matched
    - tests that we can export and read back in the dicom_data
    """

    volume = MarkerVolume((this_dir / 'test_data' / 'CT_dicom').resolve())
    assert volume.MarkerCentroids.shape[0] > 1
    assert volume.MarkerCentroids.shape[0] < 4
    # there should be 3 I think so this is a very loose test but should catch gross errors
    assert volume.dicom_data is None


def test_pandas_read_in():
    """
    test we can read in numpy data and that the marker volume is the correct shape
    """
    r_outer = 150
    test_data = np.random.rand(100, 3) * r_outer
    test_data = pd.DataFrame(test_data, columns=['x', 'y', 'z'])
    volume = MarkerVolume(test_data)
    assert volume.MarkerCentroids.shape[0] == test_data.shape[0]
    volume.perturb_marker_positions(1)  # just so this code gets executed


def test_json_read_in():
    """
    test we can read in json data and that the marker volume is the correct shape
    """
    data_loc = this_dir / 'test_data' / 'CT.mrk.json'
    volume = MarkerVolume(data_loc)
    assert volume.MarkerCentroids.shape[0] == 336


def test_marker_matching_rigid():
    """
    this tests that two volumes offset from each other by a known rigid amount then randomised
    are correctly matched
    so far it always passes, but I think there is potential for it to fail if two markers are
    generated very close together, so we may need to add some extra condition on spacing
    """
    r_outer = 150  # just to scale random numbers with
    data_offset = 2
    expected_distance = np.sqrt(3*(data_offset**2))
    test_data = np.random.rand(100, 3) * r_outer
    test_data2 = test_data + data_offset
    np.random.shuffle(test_data2)  # randomise second dataset
    test_data = pd.DataFrame(test_data, columns=['x', 'y', 'z'])
    test_data2 = pd.DataFrame(test_data2, columns=['x', 'y', 'z'])
    volume1 = MarkerVolume(test_data)
    volume2 = MarkerVolume(test_data2)

    matched_volume = MatchedMarkerVolumes(volume1, volume2, sorting_method='radial')
    assert np.allclose(matched_volume.MatchedCentroids.match_distance, expected_distance)
    matched_volume = MatchedMarkerVolumes(volume1, volume2, sorting_method='nearest')
    assert np.allclose(matched_volume.MatchedCentroids.match_distance, expected_distance)

