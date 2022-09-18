from mri_distortion_toolkit.Reports import MRI_QA_Reporter
from mri_distortion_toolkit.Reports import DefaultTestSuite
from mri_distortion_toolkit.utilities import get_dicom_data
import pandas as pd
from pathlib import Path

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'
dicom_data_loc = Path(data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json')  # previosly saved from a MarkerVolume
dicom_data = get_dicom_data(dicom_data_loc)


G_x_harmonics = pd.read_csv(data_loc / 'G_x_Harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv(data_loc / 'G_y_Harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv(data_loc / 'G_z_Harmonics.csv', index_col=0).squeeze("columns")
B0_harmonics  = pd.read_csv(data_loc / 'B0_Harmonics.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics], B0_harmonics=B0_harmonics,
                         r_outer=150, dicom_data=dicom_data, tests_to_run=DefaultTestSuite)
report.write_html_report()






