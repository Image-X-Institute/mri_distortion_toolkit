'''
this demonstrates using an image to generate harmonics, then using those harmonics
to correct the image
'''

from pathlib import Path
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
from MRI_DistortionQA import calculate_harmonics
import numpy as np
from MRI_DistortionQA.K_SpaceCorrector import KspaceDistortionCorrector
from MRI_DistortionQA.utilities import plot_matched_volume_hist

# the MR data can be downloaded from here
# https://cloudstor.aarnet.edu.au/plus/apps/files/?dir=/Shared/MRI-Linac%20Experimental%20Data/Goam2%5EMr/20220428%20MR%20Linac%5ETest/10%20gre_trans_AP_330&fileid=6603056421
distorted_data_loc = Path('/home/brendan/Documents/temp/MRILinac_DistortionCorrection/MrGoam Image Correction/MrGoam images')

# extract markers:
gt_volume = MarkerVolume('/home/brendan/Downloads/MRI_distortion_QA_sample_data/CT/slicer_centroids.mrk.json', r_max=300)
dis_volume = MarkerVolume(distorted_data_loc, n_markers_expected=336, iterative_segmentation=True)
# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, ReferenceMarkers=11)
# calculate fields
B_fields = ConvertMatchedMarkersToBz(matched_volume.MatchedCentroids, dis_volume.dicom_data)
# calculate harmonics
gradient_strength = np.array(dis_volume.dicom_data['gradient_strength'])
normalisation_factor = [1/gradient_strength[0], 1/gradient_strength[1], 1/gradient_strength[2], 1]  # this normalised gradient harmonics to 1mT/m
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(B_fields.MagneticFields,
                                                                                n_order=8,
                                                                                norm=normalisation_factor)
# correct input images
GDC = KspaceDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                Gx_Harmonics=G_x_Harmonics.harmonics,
                                Gy_Harmonics=G_y_Harmonics.harmonics,
                                Gz_Harmonics=G_z_Harmonics.harmonics,
                                ImExtension='dcm',
                                dicom_data=dis_volume.dicom_data)
GDC.correct_all_images()
GDC.save_all_images()
GDC.save_all_images_as_dicom()
# Now we have the corrected images, we can compare the original and the corrected to the GT:
corrected_volume = MarkerVolume(distorted_data_loc / 'Corrected_dcm', gaussian_image_filter_sd=1, n_markers_expected=336,
                                iterative_segmentation=True, r_max=160)
remove_ind = np.logical_and(corrected_volume.MarkerCentroids.r>=80,corrected_volume.MarkerCentroids.r <= 140)
corrected_volume.MarkerCentroids = corrected_volume.MarkerCentroids.drop(corrected_volume.MarkerCentroids.index[remove_ind])
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, ReferenceMarkers=11)
plot_matched_volume_hist([matched_volume, matched_volume_corrected], ['original', 'corrected'])
