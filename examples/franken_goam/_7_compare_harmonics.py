import numpy as np

from mri_distortion_toolkit.Harmonics import SphericalHarmonicFit
from pathlib import Path
import pandas as pd
import matplotlib as mpi


class CompareHarmonics:
    """
    This is a class mostly intended for debugging/developing.
    It provides tools to compare the differences between two harmonic series.
    """

    def __init__(self, series1, series2):


        self._series1 = pd.read_csv(series1, index_col=0).squeeze("columns")
        self._series2 = pd.read_csv(series2, index_col=0).squeeze("columns")
        n_order1 = int(np.sqrt(self._series1.shape[0] - 1))
        n_order2 = int(np.sqrt(self._series2.shape[0] - 1))
        self.series1_fit = SphericalHarmonicFit(self._series1, n_order=n_order1)
        self.series2_fit = SphericalHarmonicFit(self._series2, n_order=n_order2)
        self._construct_data_frame()

    def _construct_data_frame(self):

        self.harmonics_comparison = pd.concat([self.series1_fit.HarmonicsPk_Pk,
                                               self.series2_fit.HarmonicsPk_Pk], axis=1)
        abs_error = np.abs(np.subtract(self.series1_fit.HarmonicsPk_Pk,
                                       self.series2_fit.HarmonicsPk_Pk))
        self.harmonics_comparison['absolute_error'] = abs_error
        percent_error = np.abs(np.divide(self.series1_fit.HarmonicsPk_Pk*100,
                                         self.series2_fit.HarmonicsPk_Pk)-100)
        self.harmonics_comparison['percent_error'] = percent_error

    def plot_cut_planes(self):
        self.series1_fit.plot_cut_planes()
        self.series2_fit.plot_cut_planes()

    def plot_harmonics(self):
        self.series1_fit.plot_harmonics_pk_pk(cut_off=.005)
        self.series2_fit.plot_harmonics_pk_pk(cut_off=.005)

if __name__ == '__main__':

    data_loc1 = Path('_data/G_z_Harmonics_tra.csv')
    data_loc2 = Path('_data/G_z_Harmonics_wtf.csv')
    Gz_compare = CompareHarmonics(data_loc1, data_loc2)