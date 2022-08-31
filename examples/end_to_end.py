'''
this demonstrates using an image to generate harmonics, then using those harmonics
to correct the image
'''

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



# Data import
dis_data_loc = Path(r'C:\Users\bwhe3635\cloudstor\Shared\Goam2^Mr\20220624 QA^QA')
dis_data = {'0': '01 localiser_gre',
            '1': '02 gre_trans_AP_330',
            '2': '03 gre_trans_PA_330',
            '3': '04 gre_sag_AP_330',
            '4': '05 gre_sag_PA_330',
            '5': '06 gre_cor_RL_330',
            '6': '07 gre_cor_LR_330',
            '7': 'k_space'}

distorted_data_loc = dis_data_loc / dis_data['1'] / 'Original'
distorted_data_loc_rev = dis_data_loc / dis_data['2'] / 'Original'
gt_data_loc = Path(r'C:\Users\bwhe3635\cloudstor\MRI_distortion_QA_sample_data\CT\slicer_centroids.mrk.json')

# extract markers:
gt_volume = MarkerVolume(gt_data_loc, r_max=300)
dis_volume = MarkerVolume(distorted_data_loc, n_markers_expected=336, iterative_segmentation=True)
dis_volume_rev = MarkerVolume(distorted_data_loc_rev, n_markers_expected=336, iterative_segmentation=True)
# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, reverse_gradient_data=None, n_refernce_markers=11)
# calculate fields
B_fields = ConvertMatchedMarkersToBz(matched_volume.MatchedCentroids, dis_volume.dicom_data)
# calculate harmonics
gradient_strength = np.array(dis_volume.dicom_data['gradient_strength'])
normalisation_factor = [1 / gradient_strength[0], 1 / gradient_strength[1], 1 / gradient_strength[2],
                        1]  # this normalised gradient harmonics to 1mT/m
# normalisation_factor = [1, 1, 1, 1]
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(B_fields.MagneticFields,
                                                                                n_order=8,
                                                                                norm=normalisation_factor)
# correct input images
GDC = KspaceDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                gradient_harmonics=[G_x_Harmonics.harmonics,
                                                    G_y_Harmonics.harmonics,
                                                    G_z_Harmonics.harmonics],
                                ImExtension='dcm',
                                dicom_data=dis_volume.dicom_data,
                                correct_through_plane=True)
GDC.correct_all_images()
GDC.save_all_images()
GDC.save_all_images_as_dicom()
# Now we have the corrected images, we can compare the original and the corrected to the GT:

corrected_volume = MarkerVolume(distorted_data_loc / 'Corrected_dcm',
                                n_markers_expected=336,
                                iterative_segmentation=True, r_max=160)
remove_ind = np.logical_and(corrected_volume.MarkerCentroids.r >= 70, corrected_volume.MarkerCentroids.r <= 140)
corrected_volume.MarkerCentroids = corrected_volume.MarkerCentroids.drop(
    corrected_volume.MarkerCentroids.index[remove_ind])
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=11)
plot_matched_volume_hist([matched_volume, matched_volume_corrected], ['original', 'corrected'])
plot_distortion_xyz_hist(matched_volume_corrected)