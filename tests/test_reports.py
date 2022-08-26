from pathlib import Path
import sys
this_dir = Path(__file__).parent
sys.path.insert(0, str(this_dir.parent))
import numpy as np
from MRI_DistortionQA.Reports import MRI_QA_Reporter
from MRI_DistortionQA.Harmonics import SphericalHarmonicFit
import pandas as pd

def test_linear_recon():
    """
    this tests that when only linear gradient harmoncis are used to reconstruct data,
     we have zero distortion.
     The perfect gradients are constructed by setting all other harmonics to zero.
    """
    n_order = 10
    # create dummy data to allow us to generate dummy harmonics:
    x = np.linspace(-150, 150, 10)
    y = np.linspace(-150, 150, 10)
    z = np.linspace(-150, 150, 10)
    [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
    Bz = np.random.rand(x_recon.flatten().shape[0])
    dummy_data = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten(), 'Bz': Bz})
    dummy_harmonics = SphericalHarmonicFit(dummy_data, n_order=n_order, AssessHarmonicPk_Pk=False, QuantifyFit=False)

    dicom_data = {'FOV': [330.0, 360.0, 330.0], 'bandwidth': 260.0, 'gama': 42.57747851, 'pixel_spacing': [2.578125, 4.0, 2.578125],
     'image_size': [128, 90, 128], 'magnetic_field_strength': '1.0', 'imaging_frequency': '41.980956',
     'acquisition_date': '28_April_2022', 'manufacturer': 'SIEMENS', 'chem_shift_magnitude': 0.8723592481911194,
     'InPlanePhaseEncodingDirection': 'ROW', 'phase_encode_direction': 'x', 'freq_encode_direction': 'z',
     'slice_direction': 'y', 'gradient_strength': [1e-3, 1e-3, 1e-3]}

    bandwidth = np.array(dicom_data['bandwidth'])
    image_size = np.array(dicom_data['image_size'])
    gama = np.array(dicom_data['gama'])
    FOV = np.array(dicom_data['FOV'])
    gradient_strength = bandwidth * image_size / (gama * 1e6 * FOV * 1e-3)  # unit(T / m)

    linear_Gx_harmonics = dummy_harmonics.harmonics.copy()
    linear_Gx_harmonics[:] = 0
    linear_Gx_harmonics[2] = -1*gradient_strength[0] / np.sqrt(0.75)
    # we are scaling each gradient strength by the legendre normalisation
    linear_Gy_harmonics = dummy_harmonics.harmonics.copy()
    linear_Gy_harmonics[:] = 0
    linear_Gy_harmonics[3] = -1*gradient_strength[1] / np.sqrt(0.75)
    linear_Gz_harmonics = dummy_harmonics.harmonics.copy()
    linear_Gz_harmonics[:] = 0
    linear_Gz_harmonics[1] = -1*gradient_strength[2] / np.sqrt(1+0.5)


    report = MRI_QA_Reporter(gradient_harmonics=[linear_Gx_harmonics, linear_Gy_harmonics, linear_Gz_harmonics],
                             dicom_data=dicom_data, r_outer=150)
    assert np.allclose(0, report._MatchedMarkerVolume.abs_dis)
    report.write_html_report()  # not totally sure if github will let me do this..


