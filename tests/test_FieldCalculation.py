import sys
from pathlib import Path
import pandas as pd
import numpy as np
this_dir = Path(__file__).parent
test_data_dir = (this_dir / 'test_data').resolve()
sys.path.insert(0, str(this_dir.parent))
from mri_distortion_toolkit.FieldCalculation import ConvertMatchedMarkersToBz


def test_field_calculation():

    test_data_loc = test_data_dir / 'MatchedMarkerVolume.csv'
    test_data = pd.read_csv(test_data_loc, index_col=0).squeeze("columns")

    dicom_data = {'FOV': [330.0, 360.0, 330.0],
                  'bandwidth': 260.0,
                  'gama': 42.57747851,
                  'pixel_spacing': [2.578125, 4.0, 2.578125],
                  'image_size': [128, 90, 128],
                  'chem_shift_magnitude': 0.8723592481911194,
                  'InPlanePhaseEncodingDirection': 'ROW',
                  'phase_encode_direction': 'x',
                  'freq_encode_direction': 'z',
                  'slice_direction': 'y',
                  'gradient_strength': [0.00236859, 0.00152663, 0.00236859]}

    b_fields = ConvertMatchedMarkersToBz(test_data, dicom_data)
    expected_columns = ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz', 'B0']
    # check data looks like we expect
    assert all([el in b_fields.MagneticFields.columns for el in expected_columns])
    assert test_data.shape[0] == b_fields.MagneticFields.shape[0]

    # check the data looks like previous calculations:
    previous_data_loc = test_data_dir / 'MagneticFieldData.csv'
    previous_calc = pd.read_csv(previous_data_loc, index_col=0).squeeze("columns")
    assert np.allclose(previous_calc, b_fields.MagneticFields)


def test_bad_dicom_data():
    """
    check that putting in bad data results in an error
    """
    test_data_loc = test_data_dir / 'MatchedMarkerVolume.csv'
    test_data = pd.read_csv(test_data_loc, index_col=0).squeeze("columns")
    dicom_data = {'I_am': "missing_data"}
    try:
        b_fields = ConvertMatchedMarkersToBz(test_data, dicom_data)
    except Exception as e:
        assert isinstance(e, AttributeError)


