from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA import calculate_harmonics
from pathlib import Path

'''
the fit_harmonics.py example shows you how to fit harmonics generally, but if you are happy to use default settings
in the marker matching step you can do it the easy way directly from two volumes:
'''

# download example data and unzip:
# https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data(2)\MRI_distortion_QA_sample_data')
# ^^ update to where you put the sample data!!

# distorted centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'MR.mrk.json', verbose=False)
# ground truth centroids
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'MR.mrk.json', verbose=False, r_max=300)
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume

B0_Harmonics, G_x_Harmonics, G_y_Harmonics, G_z_Harmonics = calculate_harmonics(ground_truth_volume, distorted_volume,
                                                                                dicom_data=dicom_data_loc)

# note that B0_harmonics is None as we did not provide distorted_volume_rev to calculate_harmonics