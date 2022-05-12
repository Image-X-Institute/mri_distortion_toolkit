from MR_DistortionQA.FieldAnalysis import SphericalHarmonicFit
import pandas as pd

FieldData = pd.read_csv('Bfields.csv', index_col=0).squeeze("columns")

'''
This data contains columns ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz']
but the spherical harmonics code expects to receieve [x, y, z, Bz]
therefore, we will need to create a new dataframe with appropriately named columns
for each field we want to fit to:
'''

n_order = 8
# G_x Harmonics
GradXdata = FieldData[['x', 'y', 'z', 'B_Gx']]
GradXdata = GradXdata.rename(columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order, rOuter=150)
G_x_Harmonics.harmonics.to_csv('G_x_harmonics.csv')

# G_y Harmonics
GradYdata = FieldData[['x', 'y', 'z', 'B_Gy']]
GradYdata = GradYdata.rename(columns={"B_Gy": "Bz"})
G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order, rOuter=150)
G_y_Harmonics.harmonics.to_csv('G_y_harmonics.csv')

# G_z Harmonics
GradZdata = FieldData[['x', 'y', 'z', 'B_Gz']]
GradZdata = GradZdata.rename(columns={"B_Gz": "Bz"})
G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, rOuter=150)
G_z_Harmonics.harmonics.to_csv('G_z_harmonics.csv')

# some plotting examples
G_x_Harmonics.plot_cut_planes()
G_x_Harmonics.plot_harmonics_pk_pk(cut_off=.01)
G_x_Harmonics.print_key_harmonics(cut_off=.01)