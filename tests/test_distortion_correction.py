from pathlib import Path
import numpy as np
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from mri_distortion_toolkit.MarkerAnalysis import MatchedMarkerVolumes
from mri_distortion_toolkit.FieldCalculation import ConvertMatchedMarkersToBz
from mri_distortion_toolkit import calculate_harmonics
from mri_distortion_toolkit.utilities import plot_distortion_xyz_hist
from mri_distortion_toolkit.K_SpaceCorrector import KspaceDistortionCorrector
from mri_distortion_toolkit.utilities import plot_matched_volume_hist, print_dict
from mri_distortion_toolkit.utilities import plot_MarkerVolume_overlay

def test_2D_execution():
    """
    execute it over some test data
    """
    dis_volume = MarkerVolume(Path(r'test_data/MR_dicom').resolve())
    # correct input images
    GDC = KspaceDistortionCorrector(ImageDirectory=Path(r'test_data/MR_dicom').resolve(),
                                    gradient_harmonics=[Path(r'test_data\G_x_harmonics.csv').resolve(),
                                                        Path(r'test_data\G_y_harmonics.csv').resolve(),
                                                        Path(r'test_data\G_z_harmonics.csv').resolve()],
                                    ImExtension='dcm',
                                    dicom_data=dis_volume.dicom_data,
                                    correct_through_plane=False)
    GDC.correct_all_images()
    save_loc = Path(r'2D_distortion_correction').resolve()
    GDC.save_all_images_as_dicom(save_loc)

def test_read_in_of_corrected_data():
    """
    once
    """
    save_loc = Path(r'2D_distortion_correction').resolve()
    if not save_loc.is_dir():
        # shouldnt happen but to be safe
        test_2D_execution()
    corrected_vol = MarkerVolume(save_loc)
    distorted_vol = MarkerVolume(Path(r'test_data/MR_dicom').resolve())
    matched_vol = MatchedMarkerVolumes(corrected_vol, distorted_vol)