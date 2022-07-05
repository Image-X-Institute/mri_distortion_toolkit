from .MarkerAnalysis import MatchedMarkerVolumes
from .FieldCalculation import ConvertMatchedMarkersToBz
from .FieldAnalysis import SphericalHarmonicFit
from .utilities import get_dicom_data
import numpy as np

def calculate_harmonics(MagneticFields, n_order=8, norm=None):
    """
    This function is essentially a wrapper of convenience. Given a set of B fields and dicom_data,
    it will calculate and return the Gx, Gy, and Gz harmonics.
    Note that the gradient harmonics will be normalised to a gradient strengt of 1 mT/m

    :param MagneticFields: a pandas dataframe with columns ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz', 'B0']. B0 is optional.
        x, y, z are in mm and the fields are in any unit - the harmonics will reflect the units
    :type MagneticFields: pandas dataframe
    :param dicom_data: a dict containing required dicom data, such as that generated within the MarkerVolume class
    :type dicom_data: dict
    :param n_order: order of harmonic fit
    :type n_order: int, optional
    :param norm: Can pass a list of four valuues, the four returned harmonic lists will be normalised accordingly.
        e.g. norm = [2,2,2,2] will divide all harmonics by 2. This is useful so you can normalise the fields to some
        value
    :return: G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics
    """

    if norm is None:
        norm = [1, 1, 1, 1]
    if not np.shape(norm)[0] == 4:
        raise ValueError('norm must be a list such as [2,2,2,2], where the elements correspond to Gx, Gy, Gz, and B0')

    # calculate harmonics

    # Gx
    GradXdata = MagneticFields[['x', 'y', 'z', 'B_Gx']]
    GradXdata = GradXdata.rename(
        columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order, r_outer=150, scale=norm[0])

    # Gy
    GradYdata = MagneticFields[['x', 'y', 'z', 'B_Gy']]
    GradYdata = GradYdata.rename(
        columns={"B_Gy": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order, r_outer=150, scale=norm[1])

    # G_z
    GradZdata = MagneticFields[['x', 'y', 'z', 'B_Gz']]
    GradZdata = GradZdata.rename(
        columns={"B_Gz": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, r_outer=150, scale=norm[2])

    # B0
    try:
        B0_data = MagneticFields[['x', 'y', 'z', 'B0']]
        B0_data = B0_data.rename(
            columns={"B0": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
        B0_Harmonics = SphericalHarmonicFit(B0_data, n_order=n_order, r_outer=150, scale=norm[3])
    except KeyError:
        B0_Harmonics = None

    return G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics