from mri_distortion_toolkit import calculate_harmonics
from pathlib import Path
import pandas as pd
from mri_distortion_toolkit.utilities import get_dicom_data
import numpy as np

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_data'

scan_data_loc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/20221107 MR Linac^Test')
scans = {'9': '10 t1_tse_256_sag',
         '11': '12 t1_tse_256_tra_PA',
         '12': '13 t1_tse_256_sag_HF',
         '13': '14 t1_tse_256_sag_FH',  # for this one had to add gaussian_image_filter_sd=0.8
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR',
         '17': '18 t1_tse_256_sag_HF_rot',
         '18': '19 t1_tse_256_sag_FH_rot',
         }

FieldData = pd.read_csv(data_loc / 'rot_Bfields.csv', index_col=0).squeeze("columns")
dicom_data_loc = Path(scan_data_loc / scans['17'] / 'Original' / 'dicom_data.json')  # previosly saved from a MarkerVolume
dicom_data = get_dicom_data(dicom_data_loc)
gradient_strength = np.array(dicom_data['gradient_strength'])
normalisation_factor = [1 / gradient_strength[0], 1 / gradient_strength[1], 1 / gradient_strength[2],
                        1]  # this normalised gradient harmonics to 1mT/m
G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(FieldData,
                                                                                n_order=8,
                                                                                scale=normalisation_factor)
# save for downstream analysis:
G_x_Harmonics.harmonics.to_csv(data_loc / 'G_x_Harmonics_rot.csv')
G_y_Harmonics.harmonics.to_csv(data_loc / 'G_y_Harmonics_rot.csv')
G_z_Harmonics.harmonics.to_csv(data_loc / 'G_z_Harmonics_rot.csv')
if B0_Harmonics:  # None evaluates as False
    B0_Harmonics.harmonics.to_csv(data_loc / 'B0_Harmonics_rot.csv')
