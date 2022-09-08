from mri_distortion_toolkit.DistortionCorrection import KspaceDistortionCorrector
from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import plot_distortion_xyz_hist
from pathlib import Path
import numpy as np

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''
distorted_data_loc = Path(r'C:\Users\bwhe3635\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data\MR\04 gre_trans_AP_330')
dis_volume = MarkerVolume(distorted_data_loc)

GDC = KspaceDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                gradient_harmonics=[Path('_example_data/G_x_Harmonics.csv').resolve(),
                                                    Path('_example_data/G_y_Harmonics.csv').resolve(),
                                                    Path('_example_data/G_z_Harmonics.csv').resolve()],
                                ImExtension='dcm',
                                dicom_data=dis_volume.dicom_data,
                                correct_through_plane=True)

GDC.correct_all_images()
GDC.save_all_images()  # saves as png so you can quickly inspect results
GDC.save_all_images_as_dicom()  # saves as dicom which can be read into analysis packages.

# Assess correction:
this_file_loc = Path(__file__).parent.resolve()
gt_data_loc = this_file_loc / '_example_data' / 'CT' / 'slicer_centroids.mrk.json'
gt_volume = MarkerVolume(gt_data_loc)

corrected_volume = MarkerVolume(distorted_data_loc / 'corrected_dcm',
                                n_markers_expected=336,
                                iterative_segmentation=True, r_max=160)
remove_ind = np.logical_and(corrected_volume.MarkerCentroids.r >= 70, corrected_volume.MarkerCentroids.r <= 140)
corrected_volume.MarkerCentroids = corrected_volume.MarkerCentroids.drop(
    corrected_volume.MarkerCentroids.index[remove_ind])
matched_volume_corrected = MatchedMarkerVolumes(gt_volume, corrected_volume, n_refernce_markers=11)
plot_distortion_xyz_hist(matched_volume_corrected)