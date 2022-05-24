from MRI_DistortionQA.FieldAnalysis import SphericalHarmonicFit
import pandas as pd
from pathlib import Path
from MRI_DistortionQA.utilities import get_dicom_data
import numpy as np

FieldData = pd.read_csv('_example_data/Bfields.csv', index_col=0).squeeze("columns")
dicom_data_loc = Path('_example_data') / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume
dicom_data = get_dicom_data(dicom_data_loc)
gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3 # get gradient strength in mT/m
'''
This data contains columns ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz']
but the spherical harmonics code expects to receieve [x, y, z, Bz]
therefore, we will need to create a new dataframe with appropriately named columns
for each field we want to fit to:
'''

n_order = 8
# G_x Harmonics
GradXdata = FieldData[['x', 'y', 'z', 'B_Gx']]
GradXdata = GradXdata.rename(columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receive one field called Bz
G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[0])
G_x_Harmonics.harmonics.to_csv('_example_data/G_x_harmonics.csv')

# G_y Harmonics
GradYdata = FieldData[['x', 'y', 'z', 'B_Gy']]
GradYdata = GradYdata.rename(columns={"B_Gy": "Bz"})
G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[1])
G_y_Harmonics.harmonics.to_csv('_example_data/G_y_harmonics.csv')

# G_z Harmonics
GradZdata = FieldData[['x', 'y', 'z', 'B_Gz']]
GradZdata = GradZdata.rename(columns={"B_Gz": "Bz"})
G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[2])
G_z_Harmonics.harmonics.to_csv('_example_data/G_z_harmonics.csv')

# B0 Harmonics
GradZdata = FieldData[['x', 'y', 'z', 'B0']]
GradZdata = GradZdata.rename(columns={"B0": "Bz"})
G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[2])
G_z_Harmonics.harmonics.to_csv('_example_data/B0_harmonics.csv')

# some plotting examples
G_x_Harmonics.plot_cut_planes()
G_x_Harmonics.plot_harmonics_pk_pk(cut_off=.01)
G_x_Harmonics.print_key_harmonics(cut_off=.01)

