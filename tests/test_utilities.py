"""
like the name implies.
a lot of the functions in utilities relate to reading in data and converting it
it's hard to have ground truth tests for this, so most of these tests only check that the answer is the same
as a previous run.
Therefore, there may very well be legitimate cases where improvements are made to some function and the tests fail.
that's fine; just update the tests when that occurs.
"""
import os.path
import sys
from pathlib import Path
import pydicom
import numpy as np
import pandas as pd

this_dir = Path(__file__).parent
sys.path.insert(0, str(this_dir.parent))
import MRI_DistortionQA.utilities as ut


test_data_dir = (this_dir / 'test_data').resolve()
dicom_directory = (test_data_dir / 'MR_dicom').resolve()
test_coordinates = pd.read_csv((test_data_dir / 'test_coords.csv').resolve()).squeeze("columns")
test_coordinates_spherical = ut.convert_cartesian_to_spherical(test_coordinates)


def test_dicom_read_in():
    """
    tests all the functions related to read in of dicom files
    several of these tests are hardcoded for the dicom files in test data, and will fail if someone
    changes this data. It's not the most pure test but it should detect if someone breaks one of the functions!
    """

    dicom_files = ut.get_all_files(dicom_directory, '.dcm')
    assert dicom_files.__len__() == 11
    assert np.all([os.path.splitext(el)[1] == '.dcm' for el in dicom_files])
    CompletePathFiles = [str(Path(dicom_directory) / file) for file in dicom_files]
    dicom_slices = [pydicom.read_file(f) for f in CompletePathFiles]
    sorted_dicom_files = ut.sort_dicom_slices(dicom_slices)
    expected_order = ['MR000040.dcm', 'MR000041.dcm', 'MR000042.dcm', 'MR000043.dcm', 'MR000044.dcm', 'MR000045.dcm',
                      'MR000046.dcm', 'MR000047.dcm', 'MR000048.dcm', 'MR000049.dcm', 'MR000050.dcm']
    actual_order = [os.path.split(el.filename)[1] for el in sorted_dicom_files]
    assert np.array_equal(expected_order, actual_order)
    dicom_affine1 = ut.build_dicom_affine(sorted_dicom_files)
    expected_affine = np.array([[0., 2.578125, 0., -165.],
                                [2.578125, 0., 0., -165.],
                                [0., 0., 4., -22.],
                                [0., 0., 0., 0.]])
    assert np.allclose(dicom_affine1, expected_affine)
    DicomVolume, dicom_affine2, (X, Y, Z) = \
        ut.dicom_to_numpy(dicom_directory, FilesToReadIn=None, file_extension='dcm', return_XYZ=True)
    assert DicomVolume.shape == (128, 128, 11)
    assert np.allclose(dicom_affine1, dicom_affine2)
    assert DicomVolume.shape == X.shape == Y.shape == Z.shape


def test_marker_volume_zero_padding():

    data_loc = (this_dir / 'test_data' / 'MR_dicom').resolve()
    n_zero_pad = 2
    ImageArray, dicom_affine, (X, Y, Z) = ut.dicom_to_numpy(data_loc, file_extension='dcm',
                                                         return_XYZ=True,
                                                         zero_pad=0)
    ImageArray, dicom_affine_zp, (X_zp, Y_zp, Z_zp) = ut.dicom_to_numpy(data_loc, file_extension='dcm',
                                                         return_XYZ=True,
                                                         zero_pad=n_zero_pad)
    xyz_spacing = dicom_affine[0:3, 0:3].sum(axis=0)
    assert X_zp.min() == X.min() - xyz_spacing[0] * n_zero_pad
    assert X_zp.max() == X.max() + xyz_spacing[0] * n_zero_pad
    assert Y_zp.min() == Y.min() - xyz_spacing[1] * n_zero_pad
    assert Y_zp.max() == Y.max() + xyz_spacing[1] * n_zero_pad
    assert Z_zp.min() == Z.min() - xyz_spacing[2] * n_zero_pad
    assert Z_zp.max() == Z.max() + xyz_spacing[2] * n_zero_pad


def test_spherical_to_cartesian():
    """
    tests that the spherical coordinates are correctly converted back to cartesian
    """
    test_coordinates_2 = ut.convert_spherical_to_cartesian(test_coordinates_spherical)
    assert np.allclose(test_coordinates_2.x, test_coordinates.x)
    assert np.allclose(test_coordinates_2.y, test_coordinates.y)
    assert np.allclose(test_coordinates_2.z, test_coordinates.z)


def test_generate_legendre_basis():
    """
    test generation of a legendre basis
    - test the shape is correct
    - test that the generated basis closely matches a previously generated one
    """

    n_order = 10  # arbitrary order
    legendre = ut.generate_legendre_basis(test_coordinates_spherical, n_order)
    assert legendre.shape == (test_coordinates_spherical.shape[0], (n_order+1)**2)
    # read in old data
    legendre_test = pd.read_csv((test_data_dir / 'legendre_test.csv').resolve(), index_col=0).squeeze("columns")
    assert np.allclose(legendre, legendre_test)


def test_get_spherical_harmonics():
    """
    compare read in to ground truth.
    we need these harmonics for downstream testing so make them global variables
    """

    global Gx_Harmonics, Gy_Harmonics, Gz_Harmonics  # i need these for other tests
    # test csv
    Gx_Harmonics, Gy_Harmonics, Gz_Harmonics = ut.get_gradient_spherical_harmonics(test_data_dir / 'G_x_harmonics.csv',
                                                                          test_data_dir / 'G_y_harmonics.csv',
                                                                          test_data_dir / 'G_z_harmonics.csv')
    Gx_Harmonics2 = pd.read_csv((test_data_dir / 'G_x_harmonics.csv').resolve(), index_col=0).squeeze("columns")
    Gy_Harmonics2 = pd.read_csv((test_data_dir / 'G_y_harmonics.csv').resolve(), index_col=0).squeeze("columns")
    Gz_Harmonics2 = pd.read_csv((test_data_dir / 'G_z_harmonics.csv').resolve(), index_col=0).squeeze("columns")

    assert np.allclose(Gx_Harmonics, Gx_Harmonics2)
    assert np.allclose(Gy_Harmonics, Gy_Harmonics2)
    assert np.allclose(Gz_Harmonics, Gz_Harmonics2)


def test_reconstruct_Bz():
    """
    tests stability wrt to previous code
    """


    Bz_uT = ut.reconstruct_Bz(Gx_Harmonics, test_coordinates_spherical,
                      quantity='uT', r_outer=None)
    Bz_T = ut.reconstruct_Bz(Gy_Harmonics, test_coordinates_spherical,
                      quantity='T', r_outer=150)

    Bz_uT2 = pd.read_csv((test_data_dir / 'Bz_uT_test.csv').resolve(), index_col=0).squeeze("columns")
    Bz_T2 = pd.read_csv((test_data_dir / 'Bz_T_test.csv').resolve(), index_col=0).squeeze("columns")
    assert np.allclose(Bz_uT, Bz_uT2)
    assert np.allclose(Bz_T, Bz_T2)


