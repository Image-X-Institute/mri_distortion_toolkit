from MRI_DistortionQA import calculate_harmonics
from pathlib import Path
import pandas as pd
from MRI_DistortionQA.utilities import get_dicom_data
import numpy as np

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'


FieldData = pd.read_csv(data_loc / 'Bfields.csv', index_col=0).squeeze("columns")
dicom_data_loc = Path(data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json')  # previosly saved from a MarkerVolume
dicom_data = get_dicom_data(dicom_data_loc)
gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3
normalisation_factor = [1/gradient_strength[0], 1/gradient_strength[1], 1/gradient_strength[2], 1]  # this normalised gradient harmonics to 1mT/m

G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics = calculate_harmonics(FieldData,
                                                                                norm=normalisation_factor,
                                                                                n_order=8)

# save for downstream analysis:
G_x_Harmonics.harmonics.to_csv(data_loc / 'G_x_Harmonics.csv')
G_y_Harmonics.harmonics.to_csv(data_loc / 'G_y_Harmonics.csv')
G_z_Harmonics.harmonics.to_csv(data_loc / 'G_z_Harmonics.csv')
B0_Harmonics.harmonics.to_csv(data_loc / 'B0_Harmonics.csv')