from mri_distortion_toolkit.Harmonics import SphericalHarmonicFit
from pathlib import Path
import pandas as pd

harmonics = pd.read_csv('_data/Gx.csv', index_col=0).squeeze("columns")
Gx = SphericalHarmonicFit(harmonics)
Gx.plot_cut_planes(AddColorBar=True)

Gx.harmonics[2] = 0
Gx.plot_cut_planes(AddColorBar=True)