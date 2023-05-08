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
from mri_distortion_toolkit.DistortionCorrection import KspaceDistortionCorrector, ImageDomainDistortionCorrector
from mri_distortion_toolkit.utilities import plot_matched_volume_hist
from mri_distortion_toolkit.Reports import MRI_QA_Reporter
import matplotlib as mpi

mpi.rcParams['figure.dpi'] = 200
this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'

# Data import
distorted_data_loc = Path(
    r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\20221107 MR Linac^Test\18 t1_tse_256_sag_HF_rot\Original')
gt_data_loc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\CT\slicer_centroids.mrk.json')

# extract markers:
gt_volume = MarkerVolume(gt_data_loc, r_max=300)
gt_volume.rotate_markers(yaxis_angle=90)
gt_volume.rotate_markers(zaxis_angle=180)
gt_volume.rotate_markers(xaxis_angle=180)
dis_volume = MarkerVolume(distorted_data_loc, n_markers_expected=609,
                          iterative_segmentation=True,
                          gaussian_image_filter_sd=0.6)
# dis_volume_rev = MarkerVolume(distorted_data_loc_rev, n_markers_expected=336, iterative_segmentation=True)
# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, reverse_gradient_data=None, n_refernce_markers=9)
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
# correct input images
GDC = ImageDomainDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                     gradient_harmonics=[G_x_Harmonics.harmonics,
                                                         G_y_Harmonics.harmonics,
                                                         G_z_Harmonics.harmonics],
                                     ImExtension='dcm',
                                     dicom_data=dis_volume.dicom_data,
                                     correct_through_plane=True,
                                     correct_B0=True,
                                     B0_direction='back')
GDC.correct_all_images()
GDC.save_all_images(DSV_radius=150)
GDC.save_all_images_as_dicom()
# Now we have the corrected images, we can compare the original and the corrected to the GT:


corrected_volume = MarkerVolume(distorted_data_loc / 'corrected_dcm', n_markers_expected=609,
                                iterative_segmentation=True,
                                gaussian_image_filter_sd=0.6)
# remove_ind = np.logical_and(corrected_volume.MarkerCentroids.r >= 70, corrected_volume.MarkerCentroids.r <= 140)
# corrected_volume.MarkerCentroids = corrected_volume.MarkerCentroids.drop(
#     corrected_volume.MarkerCentroids.index[remove_ind])
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=9)
plot_matched_volume_hist([matched_volume, matched_volume_corrected], ['original', 'corrected'])
plot_distortion_xyz_hist(matched_volume_corrected)
