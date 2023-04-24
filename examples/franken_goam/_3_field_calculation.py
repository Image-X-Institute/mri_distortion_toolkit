from mri_distortion_toolkit.FieldCalculation import ConvertMatchedMarkersToBz
import pandas as pd
from pathlib import Path

this_file_loc = Path(__file__).parent.resolve()
data_csv_loc = Path('_data').resolve()

data_loc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/20221107 MR Linac^Test')
scans = {'9': '10 t1_tse_256_sag',
         '11': '12 t1_tse_256_tra_PA',
         '12': '13 t1_tse_256_sag_HF',
         '13': '14 t1_tse_256_sag_FH',  # for this one had to add gaussian_image_filter_sd=0.8
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR'}

# load the matched volume calculated in the previous step.
matched_volume = pd.read_csv(data_csv_loc / 'tra_markers.csv', index_col=0).squeeze("columns")
dicom_data = data_loc / scans['9'] / 'Original' / 'dicom_data.json'  # previously saved from a MarkerVolume
Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data)
Bz_field.MagneticFields.to_csv(data_csv_loc / 'tra_Bfields.csv')

matched_volume = pd.read_csv(data_csv_loc / 'sag_markers.csv', index_col=0).squeeze("columns")
dicom_data = data_loc / scans['12'] / 'Original' / 'dicom_data.json'  # previously saved from a MarkerVolume
Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data)
Bz_field.MagneticFields.to_csv(data_csv_loc / 'sag_Bfields.csv')

matched_volume = pd.read_csv(data_csv_loc / 'cor_markers.csv', index_col=0).squeeze("columns")
dicom_data = data_loc / scans['14'] / 'Original' / 'dicom_data.json'  # previously saved from a MarkerVolume
Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data)
Bz_field.MagneticFields.to_csv(data_csv_loc / 'cor_Bfields.csv')


