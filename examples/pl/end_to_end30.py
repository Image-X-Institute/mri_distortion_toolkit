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
from MRI_DistortionQA.DistortionCorrector import KspaceDistortionCorrector
from MRI_DistortionQA.utilities import plot_matched_volume_hist, print_dict
from MRI_DistortionQA.utilities import plot_MarkerVolume_overlay
from matplotlib import pyplot as plt

plt.rcParams["figure.dpi"] = 150  # for 4k screens

# Data import
dis_centroid_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\MR\08 gre_cor_RL_330\Original\slicer_centroids1_cross.mrk.json')

# Ground truth
gt_centroid_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\CT\slicer_centroids_fixed_cross.mrk.json')

# extract markers:
gt_volume = MarkerVolume(gt_centroid_loc, r_max=300)
gt_volume.plot_3D_markers(title='CT')

dis_volume = MarkerVolume(dis_centroid_loc)
dis_volume.plot_3D_markers(title='MR')

# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, n_refernce_markers=6)
matched_volume.plot_3D_markers()
plot_distortion_xyz_hist(matched_volume)

# calculate fields
B_fields = ConvertMatchedMarkersToBz(matched_volume.MatchedCentroids, dis_volume.dicom_data)

# calculate harmonics
gradient_strength = np.array(dis_volume.dicom_data['gradient_strength'])
normalisation_factor = [1 / gradient_strength[0], 1 / gradient_strength[1], 1 / gradient_strength[2],
                        1]  # this normalised gradient harmonics to 1mT/m
# normalisation_factor = [1,1,1,1]
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(B_fields.MagneticFields,
                                                                                n_order=8,
                                                                                norm=normalisation_factor)

# save for downstream analysis:
G_x_Harmonics.harmonics.to_csv('G_x_Harmonics.csv')
G_y_Harmonics.harmonics.to_csv('G_y_Harmonics.csv')
G_z_Harmonics.harmonics.to_csv('G_z_Harmonics.csv')

# # Try to correct itself
# distorted_data_loc = r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\MR\08 gre_cor_RL_330\Original'
#
# # correct input images
# GDC = KspaceDistortionCorrector(ImageDirectory=distorted_data_loc,
#                                 Gx_Harmonics=G_x_Harmonics.harmonics,
#                                 Gy_Harmonics=G_y_Harmonics.harmonics,
#                                 Gz_Harmonics=G_z_Harmonics.harmonics,
#                                 ImExtension='dcm',
#                                 dicom_data=dis_volume.dicom_data,
#                                 correct_through_plane=True)
# GDC.correct_all_images()
# GDC.save_all_images()
# GDC.save_all_images_as_dicom()
#
# # Now we have the corrected images, we can compare the original and the corrected to the GT:
#
# corrected_volume = MarkerVolume(distorted_data_loc / 'Corrected_dcm',
#                                 n_markers_expected=336,
#                                 iterative_segmentation=True, r_max=160)
# remove_ind = np.logical_and(corrected_volume.MarkerCentroids.r >= 70, corrected_volume.MarkerCentroids.r <= 140)
# corrected_volume.MarkerCentroids = corrected_volume.MarkerCentroids.drop(
#     corrected_volume.MarkerCentroids.index[remove_ind])
# matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, ReferenceMarkers=11)
# plot_matched_volume_hist([matched_volume, matched_volume_corrected], ['original', 'corrected'])
# plot_distortion_xyz_hist(matched_volume_corrected)
