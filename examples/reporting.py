from MR_DistortionQA.Reports import MRI_QA_Reporter
import pandas as pd
from pathlib import Path

# Direct data case: pass matched marker volume to MRI_QA_Reporter
# ---------------------------------------------------------------
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data')
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume using
Matched_Markers = pd.read_csv('Matched_Markers.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(MatchedMarkerVolume=Matched_Markers, r_outer=150, dicom_data=dicom_data_loc)
report.write_html_report()

# Harmonic case: pass harmonics to MRI_QA_Reporter so that data can be recontructed
# ----------------------------------------------------------------------------------
class CustomTestSuite:

    def test_case_1(self):
        # a test can return a bool:
        return True

    def test_case_2(self):
        # or a test can return a string:
        return "I am a string!"

    def test_case_3(self):
        # tests have access to the test data:
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=100)
        if test_data.abs_dis.max() < 2:
            return True
        else:
            return False


G_x_harmonics = pd.read_csv('G_x_harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv('G_y_harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv('G_z_harmonics.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics],
                         r_outer=150, dicom_data=dicom_data_loc, tests_to_run=CustomTestSuite)
report.write_html_report()








