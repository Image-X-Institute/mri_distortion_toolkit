import sys
from pathlib import Path
import numpy as np
import pandas as pd
this_dir = Path(__file__).parent
test_data_dir = this_dir / 'test_data'
sys.path.insert(0, str(this_dir.parent))
from MRI_DistortionQA.FieldAnalysis import SphericalHarmonicFit
from MRI_DistortionQA import calculate_harmonics
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.utilities import convert_spherical_harmonics
from MRI_DistortionQA.utilities import get_dicom_data


def test_spherical_harmonics_stability():
    """
    tests the creation of a SphericalHarmonicFit from an opera table file, and checks that
    the fitted harmonmics closely match a previously created set from the same data
    """

    data_loc = this_dir / 'test_data' / 'OperaTableTest.table'
    data = np.loadtxt(str(data_loc), skiprows=6)
    data_pd =  pd.DataFrame(data, columns = ['x','y','z','Bz'])
    harmonics = SphericalHarmonicFit(data_pd, n_order=10)

    test_harmonics = \
        pd.read_csv((this_dir / 'test_data'/ 'ground_truth_harmonics.csv').resolve(), index_col=0).squeeze("columns")
    assert np.allclose(harmonics.harmonics, test_harmonics, rtol=1e-05, atol=1e-06)


def test_calculate_harmonics():
    """
    test our wrapper function
    :return:
    """
    distorted_volume = MarkerVolume(test_data_dir / 'MR.mrk.json')
    ground_truth_volume = MarkerVolume(test_data_dir / 'CT.mrk.json')
    dicom_data_loc = test_data_dir / 'dicom_data.json'
    dicom_data = get_dicom_data(dicom_data_loc)
    gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3

    FieldData = pd.read_csv(test_data_dir / 'MagneticFieldData.csv', index_col=0).squeeze("columns")
    dicom_data_loc = Path(test_data_dir / 'dicom_data.json')  # previosly saved from a MarkerVolume
    dicom_data = get_dicom_data(dicom_data_loc)
    gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3
    normalisation_factor = [1 / gradient_strength[0], 1 / gradient_strength[1], 1 / gradient_strength[2],
                            1]  # this normalised gradient harmonics to 1mT/m

    G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(FieldData,
                                                                                    norm=normalisation_factor,
                                                                                    n_order=8)

    # can't be bothered saving in and reading all harmonics so just checking max:
    assert np.allclose(G_x_Harmonics.harmonics.max(), 0.40997688028675766)
    assert np.allclose(G_y_Harmonics.harmonics.max(), 0.7000613996928333)
    assert np.allclose(G_z_Harmonics.harmonics.max(), 0.5690328634698217)


def test_harmonic_conversion():
    """
    test that we can read in harmonics, convert them forward and back and the remain the same
    """
    g_x_harmonics = pd.read_csv(test_data_dir / 'G_x_harmonics.csv', index_col=0).squeeze("columns")
    g_x_harmonics_no_norm = convert_spherical_harmonics(g_x_harmonics, input_format='full', output_format='none')
    g_x_harmonics_norm = convert_spherical_harmonics(g_x_harmonics_no_norm, input_format='none', output_format='full')
    assert np.allclose(g_x_harmonics, g_x_harmonics_norm)
    g_x_harmonics_numpy = g_x_harmonics.to_numpy()
    g_x_harmonics_no_norm_numpy = convert_spherical_harmonics(g_x_harmonics_numpy, input_format='full', output_format='none')
    assert np.allclose(g_x_harmonics_no_norm, g_x_harmonics_no_norm_numpy)


def test_scale_harmonics():
    data_loc = this_dir / 'test_data' / 'OperaTableTest.table'
    data = np.loadtxt(str(data_loc), skiprows=6)
    data_pd =  pd.DataFrame(data, columns = ['x', 'y', 'z', 'Bz'])
    harmonics = SphericalHarmonicFit(data_pd, n_order=10)
    harmonics_scaled = SphericalHarmonicFit(data_pd, n_order=10, scale=10)
    scale_factor = np.mean(np.divide(harmonics_scaled.harmonics, harmonics.harmonics))
    assert np.isclose(10, scale_factor)