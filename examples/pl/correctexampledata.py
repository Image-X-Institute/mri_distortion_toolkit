'''
this demonstrates using an image to generate harmonics, then using those harmonics
to correct the image
'''

from pathlib import Path
import numpy as np
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
from MRI_DistortionQA import calculate_harmonics
from MRI_DistortionQA.utilities import plot_distortion_xyz_hist
from MRI_DistortionQA.DistortionCorrector import KspaceDistortionCorrector, ImageDomainDistortionCorrector
from MRI_DistortionQA.utilities import plot_matched_volume_hist, print_dict
from MRI_DistortionQA.utilities import plot_MarkerVolume_overlay
from matplotlib import pyplot as plt

plt.rcParams["figure.dpi"] = 150  # for 4k screens

## Correct itself
data_loc = Path(
    r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220413 Full sphere foam Goam\20220413 MR Linac^Test')

# Ground truth
gt_volume = MarkerVolume(data_loc / 'slicer_centroids.mrk.json', r_max=300)
# gt_volume.plot_3D_markers(title='CT')

# Uncorrected markers
uncorrected_volume = MarkerVolume(data_loc / '05 gre_trans_PA_330' / 'Original' / 'slicer_centroids.mrk.json')
# uncorrected_volume.plot_3D_markers(title='Uncorrected')

uncorrected_match = matched_volume = MatchedMarkerVolumes(gt_volume, uncorrected_volume, n_refernce_markers=11)
# uncorrected_match.plot_3D_markers()
# plot_distortion_xyz_hist(uncorrected_match)


# # Correct images (kspace example)
# GDC = KspaceDistortionCorrector(ImageDirectory=data_loc / '05 gre_trans_PA_330' / 'Original',
#                                 Gx_Harmonics='G_x_Harmonics.csv',
#                                 Gy_Harmonics='G_y_Harmonics.csv',
#                                 Gz_Harmonics='G_z_Harmonics.csv',
#                                 ImExtension='dcm',
#                                 dicom_data=uncorrected_volume.dicom_data,
#                                 correct_through_plane=True)
# GDC.correct_all_images()
# # GDC.save_all_images_as_dicom()

# Correct images (image domain example)
IDC = ImageDomainDistortionCorrector(ImageDirectory=data_loc / '05 gre_trans_PA_330' / 'Original',
                                     Gx_Harmonics='G_x_Harmonics.csv',
                                     Gy_Harmonics='G_y_Harmonics.csv',
                                     Gz_Harmonics='G_z_Harmonics.csv',
                                     ImExtension='dcm',
                                     dicom_data=uncorrected_volume.dicom_data,
                                     correct_through_plane=False)
IDC.correct_all_images()
IDC.save_all_images()

# # Segment corrected
# corrected_volume = MarkerVolume(data_loc / '05 gre_trans_PA_330' / 'Original' / 'Corrected_dcm',
#                                 n_markers_expected=335, iterative_segmentation=True, r_max=165,
#                                 gaussian_image_filter_sd=0.7, verbose=True)
# corrected_volume.plot_3D_markers(title='Corrected')
#
# corrected_match = matched_volume = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=11)
# corrected_match.plot_3D_markers()
# plot_distortion_xyz_hist(corrected_match)
