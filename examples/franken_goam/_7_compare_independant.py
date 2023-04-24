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
gt_data_loc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\CT\slicer_centroids.mrk.json')

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

distorted_data_loc_rev = data_loc / scans['18'] / 'Original' / 'slicer_centroids.mrk.json'


_orientations = ['tra', 'sag', 'cor']
orientation = 'tra'
_correct_B0 = True
_correct_through_plane = True
assert orientation in _orientations
if orientation == 'tra':
    distorted_data_loc = data_loc / scans['9'] / 'Original' / 'slicer_centroids.mrk.json'
    dis_volume = MarkerVolume(distorted_data_loc)
    _Gx = Path('_data/G_x_Harmonics_cor.csv')
    _Gy = Path('_data/G_y_Harmonics_sag.csv')
    _Gz = Path('_data/G_z_Harmonics_cor.csv')
    _B0 = Path('_data/B0_Harmonics_sag.csv')
    _B0_direction = 'forward'
if orientation == 'sag':
    distorted_data_loc = data_loc / scans['12'] / 'Original' / 'slicer_centroids.mrk.json'
    dis_volume = MarkerVolume(distorted_data_loc)
    _Gx = Path('_data/G_x_Harmonics_tra.csv')
    _Gy = Path('_data/G_y_Harmonics_tra.csv')
    _Gz = Path('_data/G_z_Harmonics_cor.csv')
    _B0 = Path('_data/B0_Harmonics_cor.csv')
    _B0_direction = 'forward'
if orientation == 'cor':
    distorted_data_loc = data_loc / scans['14'] / 'Original' / 'slicer_centroids.mrk.json'
    dis_volume = MarkerVolume(distorted_data_loc)
    _Gx = Path('_data/G_x_Harmonics_tra.csv')
    _Gy = Path('_data/G_y_Harmonics_tra.csv')
    _Gz = Path('_data/G_z_Harmonics_sag.csv')
    _B0 = Path('_data/B0_Harmonics_sag.csv')
    _B0_direction = 'forward'

gt_volume = MarkerVolume(gt_data_loc, r_max=300)
gt_volume.rotate_markers(yaxis_angle=180)

dis_volume = MarkerVolume(distorted_data_loc)
uncorrected_match = MatchedMarkerVolumes(dis_volume, gt_volume)
uncorrected_match.report()
plot_distortion_xyz_hist(uncorrected_match)

# correct input images
GDC = ImageDomainDistortionCorrector(ImageDirectory=distorted_data_loc.parent.resolve(),
                                gradient_harmonics=[_Gx, _Gy, _Gz],
                                B0_harmonics=_B0,
                                dicom_data=dis_volume.dicom_data,
                                correct_through_plane=_correct_through_plane,
                                correct_B0=_correct_B0,
                                B0_direction=_B0_direction)
print(dis_volume.dicom_data)
GDC.correct_all_images()
# GDC.save_all_images()
GDC.save_all_images_as_dicom()
# Now we have the corrected images, we can compare the original and the corrected to the GT:
gt_volume = MarkerVolume(gt_data_loc, r_max=300)
gt_volume.rotate_markers(yaxis_angle=180)
nmarkers = 618
corrected_volume = MarkerVolume(distorted_data_loc.parent.resolve() / 'corrected_dcm',
                                n_markers_expected=nmarkers,
                                iterative_segmentation=False,
                                threshold=323,
                                gaussian_image_filter_sd=0.8)
corrected_volume.export_to_slicer()
corrected_volume.save_dicom_data()
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=9)
matched_volume_corrected.report()
plot_matched_volume_hist([matched_volume_corrected])
plot_distortion_xyz_hist(matched_volume_corrected)