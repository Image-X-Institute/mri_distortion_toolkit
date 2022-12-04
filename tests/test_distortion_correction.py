from pathlib import Path
import sys
this_dir = Path(__file__).parent
sys.path.insert(0, str(this_dir.parent))

from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from mri_distortion_toolkit.MarkerAnalysis import MatchedMarkerVolumes
from mri_distortion_toolkit.DistortionCorrection import KspaceDistortionCorrector
from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.utilities import generate_harmonic_names
import numpy as np

test_data_dir = (this_dir / 'test_data').resolve()

def test_k_space_corrector_runs():
    """
    execute it over some test data
    """
    assert (this_dir / 'test_data' / 'MR_dicom').resolve().is_dir()
    assert (this_dir / 'test_data' / 'G_x_Harmonics.csv').resolve().is_file()
    dis_volume = MarkerVolume(test_data_dir / 'MR_dicom')
    # correct input images
    GDC = KspaceDistortionCorrector(ImageDirectory=(test_data_dir / 'MR_dicom').resolve(),
                                    gradient_harmonics=[(test_data_dir / 'G_x_Harmonics.csv').resolve(),
                                                        (test_data_dir / 'G_y_Harmonics.csv').resolve(),
                                                        (test_data_dir / 'G_z_Harmonics.csv').resolve()],
                                    ImExtension='dcm',
                                    correct_through_plane=True)
    GDC.correct_all_images()
    GDC.save_all_images(save_loc=this_dir / 'test_data' / 'MR_dicom' / 'k_space_corrected_png')
    GDC.save_all_images_as_dicom(save_loc=this_dir / 'test_data' / 'MR_dicom' / 'k_space_corrected_dcm')


def test_image_domain_corrector_runs():
    """
    execute it over some test data
    """
    assert (this_dir / 'test_data' / 'MR_dicom').resolve().is_dir()
    assert (this_dir / 'test_data' / 'G_x_Harmonics.csv').resolve().is_file()
    dis_volume = MarkerVolume(test_data_dir / 'MR_dicom')
    # correct input images
    GDC = ImageDomainDistortionCorrector(ImageDirectory=(test_data_dir / 'MR_dicom').resolve(),
                                    gradient_harmonics=[(test_data_dir / 'G_x_Harmonics.csv').resolve(),
                                                        (test_data_dir / 'G_y_Harmonics.csv').resolve(),
                                                        (test_data_dir / 'G_z_Harmonics.csv').resolve()],
                                    ImExtension='dcm',
                                    correct_through_plane=True)
    GDC.correct_all_images()
    GDC.save_all_images_as_dicom(save_loc=this_dir / 'test_data' / 'MR_dicom' / 'im_domain_corrected_dcm')


def test_k_space_corrected_data():
    """
    test we can read in the corrected data, and test for stability of detected distortion
    """
    corrected_data = test_data_dir / 'MR_dicom' / 'k_space_corrected_dcm'
    if not corrected_data.is_dir():
        # shouldnt happen but to be safe
        test_k_space_corrector_runs()
    corrected_vol = MarkerVolume(corrected_data, threshold=39)
    distorted_vol = MarkerVolume(test_data_dir / 'MR_dicom')
    matched_vol = MatchedMarkerVolumes(corrected_vol, distorted_vol)
    assert np.mean(matched_vol.MatchedCentroids.match_distance) < 12


def test_image_domain_corrected_data():
    """
    test we can read in the corrected data, and test for stability of detected distortion
    """
    corrected_data = test_data_dir / 'MR_dicom' / 'im_domain_corrected_dcm'
    if not corrected_data.is_dir():
        # shouldnt happen but to be safe
        test_image_domain_corrector_runs()
    corrected_vol = MarkerVolume(corrected_data)
    distorted_vol = MarkerVolume(test_data_dir / 'MR_dicom')
    matched_vol = MatchedMarkerVolumes(corrected_vol, distorted_vol)
    print(np.mean(matched_vol.MatchedCentroids.match_distance))
    assert np.mean(matched_vol.MatchedCentroids.match_distance) < 10

