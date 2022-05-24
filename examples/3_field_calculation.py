from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
import pandas as pd
from pathlib import Path


'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''
data_loc = Path(r'/home/brendan/Downloads/MRI_distortion_QA_sample_data')

# load the matched volume calculated in the previous step.
matched_volume = pd.read_csv('_example_data/Matched_Markers.csv', index_col=0).squeeze("columns")
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume

Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data_loc)
Bz_field.MagneticFields.to_csv('_example_data/Bfields.csv')

