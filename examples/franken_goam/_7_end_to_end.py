'''
this rather long script demonstrates using an image to generate harmonics, using those harmonics
to correct the image, analysing the result, and generating a before/after report
'''

from pathlib import Path
import numpy as np
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from mri_distortion_toolkit.MarkerAnalysis import MatchedMarkerVolumes
from mri_distortion_toolkit.FieldCalculation import ConvertMatchedMarkersToBz
from mri_distortion_toolkit import calculate_harmonics
from mri_distortion_toolkit.utilities import plot_distortion_xyz_hist
from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.utilities import plot_matched_volume_hist
from mri_distortion_toolkit.utilities import plot_compressed_MarkerVolumes
from mri_distortion_toolkit.utilities import plot_MarkerVolume_overlay

data_loc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\20221107 MR Linac^Test')
scans = {'0': '01 localiser_gre',
         '1': '02 localiser_gre',
         '2': '03 localiser_gre',
         '3': '04 gre_trans_AP_330',
         '4': '05 gre_trans_PA_330',
         '5': '06 gre_sag_AP_330',
         '6': '07 gre_sag_PA_330',
         '7': '08 gre_cor_RL_330',
         '8': '09 gre_cor_LR_330',
         '9': '10 t1_tse_256_sag',
         '10': '11 t1_tse_256_sag_PA',
         '11': '12 t1_tse_256_tra_PA',
         '12': '13 t1_tse_256_sag_HF',
         '13': '14 t1_tse_256_sag_FH',  # for this one had to add gaussian_image_filter_sd=0.8
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR',
         '16': '17 localiser_gre',
         '17': '18 t1_tse_256_sag_HF_rot',
         '18': '19 t1_tse_256_sag_FH_rot',
         '19': '20 trufi_sag_128_torsocoil',
         '20': '21 trufi_sag_128_torsocoil',
         '21': '22 trufi_sag_128_torsocoil'}

# process TSE images
scans_to_segment = ['9', '11', '12', '13', '14']

# Data import
distorted_data_loc = data_loc / scans['17'] / 'Original' / 'slicer_centroids.mrk.json'
distorted_data_loc_rev = data_loc / scans['18'] / 'Original' / 'slicer_centroids.mrk.json'
gt_data_loc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\CT\slicer_centroids.mrk.json')

# extract markers:
gt_volume = MarkerVolume(gt_data_loc, r_max=300)

if False:
    gt_volume.rotate_markers(yaxis_angle=180)
    nmarkers = 618
else:
    gt_volume.rotate_markers(yaxis_angle=90)
    gt_volume.rotate_markers(zaxis_angle=180)
    gt_volume.MarkerCentroids = gt_volume.MarkerCentroids[gt_volume.MarkerCentroids.x > -150]
    gt_volume.translate_markers(x_shift=-14, y_shift=55, z_shift=-3)
    gt_volume.rotate_markers(xaxis_angle=180)
    nmarkers = 609
dis_volume = MarkerVolume(distorted_data_loc, n_markers_expected=nmarkers)

# plot_compressed_MarkerVolumes([dis_volume, gt_volume], projection_direction='x')
# plot_MarkerVolume_overlay([dis_volume, gt_volume], legend=['distorted', 'gt'])

dis_volume_rev = MarkerVolume(distorted_data_loc_rev, n_markers_expected=nmarkers)
# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, reverse_gradient_data=dis_volume_rev, n_refernce_markers=9)
# calculate fields
B_fields = ConvertMatchedMarkersToBz(matched_volume.MatchedCentroids, dis_volume.dicom_data)
# calculate harmonics
gradient_strength = np.array(dis_volume.dicom_data['gradient_strength'])
normalisation_factor = [1 / gradient_strength[0], 1 / gradient_strength[1], 1 / gradient_strength[2],
                        1]  # this normalised gradient harmonics to 1mT/m
# normalisation_factor = [1, 1, 1, 1]
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(B_fields.MagneticFields,
                                                                                n_order=5,
                                                                                scale=normalisation_factor)
G_x_Harmonics.harmonics.to_csv('_data/G_x_Harmonics_rot.csv')
G_y_Harmonics.harmonics.to_csv('_data/G_y_Harmonics_rot.csv')
G_z_Harmonics.harmonics.to_csv('_data/G_z_Harmonics_rot.csv')
if B0_Harmonics:  # None evaluates as False
    B0_Harmonics.harmonics.to_csv('_data/B0_Harmonics_rot.csv')
# correct input images
GDC = ImageDomainDistortionCorrector(ImageDirectory=distorted_data_loc.parent.resolve(),
                                gradient_harmonics=[Path('_data/Gx.csv'),
                                                    Path('_data/Gy.csv'),
                                                    Path('_data/Gz.csv')],
                                B0_harmonics=Path('_data/B0.csv'),
                                ImExtension='dcm',
                                dicom_data=dis_volume.dicom_data,
                                correct_through_plane=True,
                                correct_B0=True,
                                B0_direction='back')
GDC.correct_all_images()
# GDC.save_all_images()
GDC.save_all_images_as_dicom()
# Now we have the corrected images, we can compare the original and the corrected to the GT:

corrected_volume = MarkerVolume(distorted_data_loc.parent.resolve() / 'corrected_dcm',
                                n_markers_expected=nmarkers,
                                iterative_segmentation=True,
                                threshold=None,
                                gaussian_image_filter_sd=0.8)
corrected_volume.export_to_slicer()
corrected_volume.save_dicom_data()
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=9)
matched_volume_corrected.report()
plot_matched_volume_hist([matched_volume_corrected])
plot_distortion_xyz_hist(matched_volume_corrected)
