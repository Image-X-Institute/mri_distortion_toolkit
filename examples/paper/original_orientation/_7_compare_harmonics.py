import numpy as np

from mri_distortion_toolkit.Harmonics import SphericalHarmonicFit
from pathlib import Path
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


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
        percent_error = np.abs(np.divide(self.series1_fit.HarmonicsPk_Pk * 100,
                                         self.series2_fit.HarmonicsPk_Pk) - 100)
        self.harmonics_comparison['percent_error'] = percent_error

    def plot_cut_planes(self):
        self.series1_fit.plot_cut_planes()
        self.series2_fit.plot_cut_planes()

    def plot_harmonics(self, cut_off=.005, label_data=True, y_range=None):

        plt.figure(figsize=[10, 5])
        # ax = sns.barplot(self.HarmonicsPk_Pk, y='pk-pk [\u03BCT]', x=)
        CutOffInd1 = abs(self.series1_fit.HarmonicsPk_Pk) < cut_off * abs(self.series1_fit.HarmonicsPk_Pk).max()
        CutOffInd2 = abs(self.series2_fit.HarmonicsPk_Pk) < cut_off * abs(self.series2_fit.HarmonicsPk_Pk).max()
        CutOffInd = np.logical_not(np.logical_or(np.logical_not(CutOffInd1), np.logical_not(CutOffInd2)))
        HarmonicsToPlot_1 = self.series1_fit.HarmonicsPk_Pk.drop(self.series1_fit.HarmonicsPk_Pk[CutOffInd].index)
        HarmonicsToPlot_2 = self.series2_fit.HarmonicsPk_Pk.drop(self.series2_fit.HarmonicsPk_Pk[CutOffInd].index)

        HarmonicsToPlot_1 = HarmonicsToPlot_1.to_frame()
        HarmonicsToPlot_2 = HarmonicsToPlot_2.to_frame()
        HarmonicsToPlot_1['series'] = 'one'
        HarmonicsToPlot_2['series'] = 'two'
        HarmonicsToPlot = pd.concat([HarmonicsToPlot_1, HarmonicsToPlot_2])
        HarmonicsToPlot = HarmonicsToPlot.reset_index()
        HarmonicsToPlot = HarmonicsToPlot.rename(columns={0: "pk_pk"})

        axs = sns.barplot(data=HarmonicsToPlot, x="index", y="pk_pk", hue="series")
        if y_range:
            axs.set_ylim(y_range)
        axs.set_title(f'Principle Harmonics pk-pk (>{cut_off * 100: 1.1f}% of max)')

        axs.set_ylabel('pk-pk [\u03BCT]')
        for item in axs.get_xticklabels():
            item.set_rotation(45)
        if label_data:
            for i in axs.containers:
                axs.bar_label(i, )
        plt.show()

if __name__ == '__main__':
    data_loc1 = Path('_data/G_x_Harmonics_rot.csv')
    data_loc2 = Path('_data/Gx.csv')
    Gy_compare = CompareHarmonics(data_loc1, data_loc2)
    Gy_compare.plot_harmonics(cut_off=.005, label_data=False, y_range=[-10000, 10000])
