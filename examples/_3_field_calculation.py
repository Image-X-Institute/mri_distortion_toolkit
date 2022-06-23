from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
import pandas as pd
from pathlib import Path

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'

# load the matched volume calculated in the previous step.
matched_volume = pd.read_csv(data_loc / 'Matched_Markers.csv', index_col=0).squeeze("columns")
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previously saved from a MarkerVolume

Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data_loc)
Bz_field.MagneticFields.to_csv(data_loc / 'Bfields.csv')

