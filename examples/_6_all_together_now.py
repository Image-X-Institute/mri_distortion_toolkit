from pathlib import Path
import numpy as np
from MRI_DistortionQA.Reports import MRI_QA_Reporter
from MRI_DistortionQA.Reports import DefaultTestSuite
from MRI_DistortionQA import calculate_harmonics
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
from MRI_DistortionQA.utilities import get_dicom_data

"""
This example encapsulate all the previous examples; instead of doing things one file at a time we do them all 
in one script
"""

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'
dicom_data = get_dicom_data(dicom_data_loc)
# marker centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json', verbose=False)
distorted_volume_rev = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330' / 'slicer_centroids.mrk.json', verbose=False)
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)

# matched volumes including reversed gradient data
matched_volume_with_rev_data = MatchedMarkerVolumes(ground_truth_volume, distorted_volume,
                                                    ReverseGradientData=distorted_volume_rev, ReferenceMarkers=11)

# Field calculation:
Bz_field = ConvertMatchedMarkersToBz(matched_volume_with_rev_data.MatchedCentroids, dicom_data_loc)

# Harmonic calculation
gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3 # get gradient strength in mT/m
normalisation_factor = [1/gradient_strength[0], 1/gradient_strength[1], 1/gradient_strength[2], 1]
n_order = 8
# this normalised gradient harmonics to 1mT/m
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(Bz_field.MagneticFields,
                                                                                norm=normalisation_factor,
                                                                                n_order=n_order)
# Reporting
report = MRI_QA_Reporter(gradient_harmonics=[G_x_Harmonics.harmonics, G_y_Harmonics.harmonics, G_z_Harmonics.harmonics],
                         B0_harmonics=B0_Harmonics.harmonics, r_outer=150, dicom_data=dicom_data_loc,
                         tests_to_run=DefaultTestSuite)
report.write_html_report()