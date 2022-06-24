from MRI_DistortionQA.Reports import MRI_QA_Reporter
from MRI_DistortionQA.Reports import DefaultTestSuite
import pandas as pd
from pathlib import Path

# Direct data case: pass matched marker volume to MRI_QA_Reporter
# ---------------------------------------------------------------
data_loc = Path(r'_example_data').resolve()
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume using
Matched_Markers = pd.read_csv('_example_data/Matched_Markers.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(MatchedMarkerVolume=Matched_Markers, r_outer=150, dicom_data=dicom_data_loc)
report.write_html_report()

# Harmonic case: pass harmonics to MRI_QA_Reporter so that data can be recontructed
# ----------------------------------------------------------------------------------

G_x_harmonics = pd.read_csv('_example_data/G_x_harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv('_example_data/G_y_harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv('_example_data/G_z_harmonics.csv', index_col=0).squeeze("columns")
B0_harmonics  = pd.read_csv('_example_data/B0_harmonics.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics], B0_harmonics=B0_harmonics,
                         r_outer=150, dicom_data=dicom_data_loc, tests_to_run=DefaultTestSuite)
# report.write_html_report()








