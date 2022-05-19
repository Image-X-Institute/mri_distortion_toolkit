import os.path
import sys
from pathlib import Path
import pydicom
import numpy as np
import pandas as pd


this_dir = Path(__file__).parent
sys.path.insert(0, str(this_dir.parent))
from MRI_DistortionQA.FieldAnalysis import SphericalHarmonicFit


def test_spherical_harmonics_stability():
    """
    tests the creation of a SphericalHarmonicFit from an opera table file, and checks that
    the fitted harmonmics closely match a previously created set from the same data
    """

    data_loc = this_dir / 'test_data' / 'OperaTableTest.table'
    data = np.loadtxt(str(data_loc), skiprows=6)
    data_pd =  pd.DataFrame(data, columns = ['x','y','z','Bz'])
    harmonics = SphericalHarmonicFit(data_pd, n_order=10)

    test_harmonics = \
        pd.read_csv((this_dir / 'test_data'/ 'ground_truth_harmonics.csv').resolve(), index_col=0).squeeze("columns")
    assert np.allclose(harmonics.harmonics, test_harmonics, rtol=1e-05, atol=1e-06)

def test_calculate_harmonics():
    """
    test our wrapper function
    :return:
    """