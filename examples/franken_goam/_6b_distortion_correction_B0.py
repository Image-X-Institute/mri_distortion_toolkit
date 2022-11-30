from mri_distortion_toolkit.DistortionCorrection import KspaceDistortionCorrector
from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import plot_distortion_xyz_hist
from mri_distortion_toolkit.utilities import get_dicom_data
from pathlib import Path
import numpy as np

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''

distorted_data_loc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/20221107 MR Linac^Test/18 t1_tse_256_sag_HF_rot/Original')
dicom_data = get_dicom_data(distorted_data_loc / 'dicom_data.json')
GDC = ImageDomainDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                     gradient_harmonics=[Path('_data/Gx.csv').resolve(),
                                                         Path('_data/Gy.csv').resolve(),
                                                         Path('_data/Gz.csv').resolve()],
                                     B0_harmonics=Path('_data/B0.csv').resolve(),
                                     dicom_data=dicom_data,
                                     ImExtension='dcm',
                                     correct_through_plane=True,
                                     correct_B0=True,
                                     B0_direction='back')

GDC.correct_all_images()
GDC.save_all_images(DSV_radius=150)  # saves as png so you can quickly inspect results
GDC.save_all_images_as_dicom()  # saves as dicom which can be read into analysis packages.

# Assess correction:
gt_data_loc = Path(r'/home/brendan/HairyHome/Downloads/FrankenGoam^Mr/CT/slicer_centroids.mrk.json')
corrected_volume = MarkerVolume(distorted_data_loc / 'corrected_dcm', n_markers_expected=609,
                                iterative_segmentation=True,
                                gaussian_image_filter_sd=0.6)
gt_volume = MarkerVolume(gt_data_loc, r_max=300)
gt_volume.rotate_markers(yaxis_angle=90)
gt_volume.rotate_markers(zaxis_angle=180)
gt_volume.rotate_markers(xaxis_angle=180)
# gt_volume.MarkerCentroids = gt_volume.MarkerCentroids.drop(gt_volume.MarkerCentroids.index[remove_ind])
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=9)
# matched_volume_uncorrected = MatchedMarkerVolumes(gt_volume, uncorrected_volume, n_refernce_markers=9)
# matched_volume_corrected_versus_uncorrected = MatchedMarkerVolumes(uncorrected_volume, corrected_volume , n_refernce_markers=9)
matched_volume_corrected.report()
# matched_volume_uncorrected.report()
plot_distortion_xyz_hist(matched_volume_corrected)

