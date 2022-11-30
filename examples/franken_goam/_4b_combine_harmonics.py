import pandas as pd
from pathlib import Path
from mri_distortion_toolkit.Harmonics import SphericalHarmonicFit

data_loc = Path('_data')

# ok, for each harmonic we should have 2 estimates; 1 frequency and 1 phase

# tra:
# phase_encode_direction: y
# freq_encode_direction: x

# sag
# phase_encode_direction: z
# freq_encode_direction: y

# cor
# phase_encode_direction: x
# freq_encode_direction: z

# Gx_Harmonics:
G_x1 = pd.read_csv(data_loc / 'G_x_Harmonics_tra.csv', index_col=0).squeeze("columns")
G_x2 = pd.read_csv(data_loc / 'G_x_Harmonics_cor.csv', index_col=0).squeeze("columns")
G_x = (G_x1 + G_x2) / 2
G_x.to_csv('_data/Gx.csv')
# input_data = pd.DataFrame({'x': [0, 0, 0], 'y': [0, 0, 0], 'z': [0, 0, 0], 'Bz': [0, 1, 2]})
# dummy_harmonics1 = SphericalHarmonicFit(input_data)
# dummy_harmonics1.harmonics = G_x1
# dummy_harmonics1._assess_harmonic_pk_pk()
# dummy_harmonics1.plot_harmonics_pk_pk(cut_off=.005)
# dummy_harmonics2 = SphericalHarmonicFit(input_data)
# dummy_harmonics2.harmonics = G_x2
# dummy_harmonics2._assess_harmonic_pk_pk()
# dummy_harmonics2.plot_harmonics_pk_pk(cut_off=.005)

# Gy_Harmonics
G_y1 = pd.read_csv(data_loc / 'G_y_Harmonics_tra.csv', index_col=0).squeeze("columns")
G_y2 = pd.read_csv(data_loc / 'G_y_Harmonics_sag.csv', index_col=0).squeeze("columns")
G_y = (G_y1 + G_y2) / 2
G_y.to_csv('_data/Gy.csv')
input_data = pd.DataFrame({'x': [0, 0, 0], 'y': [0, 0, 0], 'z': [0, 0, 0], 'Bz': [0, 1, 2]})
# dummy_harmonics1 = SphericalHarmonicFit(input_data)
# dummy_harmonics1.harmonics = G_y1
# dummy_harmonics1._assess_harmonic_pk_pk()
# dummy_harmonics1.plot_harmonics_pk_pk(cut_off=.005)
# dummy_harmonics2 = SphericalHarmonicFit(input_data)
# dummy_harmonics2.harmonics = G_y2
# dummy_harmonics2._assess_harmonic_pk_pk()
# dummy_harmonics2.plot_harmonics_pk_pk(cut_off=.005)

# Gz_Harmonics
G_z1 = pd.read_csv(data_loc / 'G_z_Harmonics_sag.csv', index_col=0).squeeze("columns")
G_z2 = pd.read_csv(data_loc / 'G_z_Harmonics_cor.csv', index_col=0).squeeze("columns")
G_z = (G_z1 + G_z2) / 2
G_z.to_csv('_data/Gz.csv')
input_data = pd.DataFrame({'x': [0, 0, 0], 'y': [0, 0, 0], 'z': [0, 0, 0], 'Bz': [0, 1, 2]})
# dummy_harmonics1 = SphericalHarmonicFit(input_data)
# dummy_harmonics1.harmonics = G_z1
# dummy_harmonics1._assess_harmonic_pk_pk()
# dummy_harmonics1.print_key_harmonics()
# dummy_harmonics2 = SphericalHarmonicFit(input_data)
# dummy_harmonics2.harmonics = G_z2
# dummy_harmonics2._assess_harmonic_pk_pk()
# dummy_harmonics2.print_key_harmonics()

# B0
B01 = pd.read_csv(data_loc / 'B0_Harmonics_sag.csv', index_col=0).squeeze("columns")
B02 = pd.read_csv(data_loc / 'B0_Harmonics_cor.csv', index_col=0).squeeze("columns")
B03 = pd.read_csv(data_loc / 'B0_Harmonics_tra.csv', index_col=0).squeeze("columns")
G_x = (B01 + B02 + B03) / 3
B01.to_csv('_data/B0.csv')
# input_data = pd.DataFrame({'x': [0, 0, 0], 'y': [0, 0, 0], 'z': [0, 0, 0], 'Bz': [0, 1, 2]})
# dummy_harmonics1 = SphericalHarmonicFit(input_data)
# dummy_harmonics1.harmonics = B01
# dummy_harmonics1._assess_harmonic_pk_pk()
# dummy_harmonics1.print_key_harmonics()
# dummy_harmonics2 = SphericalHarmonicFit(input_data)
# dummy_harmonics2.harmonics = B02
# dummy_harmonics2._assess_harmonic_pk_pk()
# dummy_harmonics2.print_key_harmonics()
# dummy_harmonics3 = SphericalHarmonicFit(input_data)
# dummy_harmonics3.harmonics = B03
# dummy_harmonics3._assess_harmonic_pk_pk()
# dummy_harmonics3.print_key_harmonics()