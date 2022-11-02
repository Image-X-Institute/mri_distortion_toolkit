from .Harmonics import SphericalHarmonicFit
import numpy as np

def calculate_harmonics(MagneticFields, n_order=5, scale=None):
    """
    This function is essentially a wrapper of convenience. Given a set of B fields and dicom_data,
    it will calculate and return the Gx, Gy, and Gz harmonics.

    :param MagneticFields: a pandas dataframe with columns ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz', 'B0']. B0 is optional.
        x, y, z are in mm and the fields are in any unit - the harmonics will reflect the units
    :type MagneticFields: pandas dataframe
    :param n_order: order of harmonic fit
    :type n_order: int or list, optional
    :param scale: Can pass a list of four values, the four returned harmonic lists will be scalealised accordingly.
        e.g. scale = [2,2,2,2] will multiply all harmonics by 2. This is useful to normalise the fields to some
        value
    :return: G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics
    """

    if scale is None:
        scale = [1, 1, 1, 1]
    if not np.shape(scale)[0] == 4:
        raise ValueError('scale must be a list such as [2,2,2,2], where the elements correspond to Gx, Gy, Gz, and B0')

    if isinstance(n_order, int):
        n_order = [n_order, n_order, n_order, n_order]
    elif not len(n_order) == 4:
        raise ValueError('n_order must either be a single integer or a list of 4 integers corresponding to'
                         '[n_order_Gx,n_order_Gy,n_order_Gz,n_order_B0]')

    # calculate harmonics
    # Gx
    GradXdata = MagneticFields[['x', 'y', 'z', 'B_Gx']]
    GradXdata = GradXdata.rename(
        columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order[0], r_outer=150, scale=scale[0])

    # Gy
    GradYdata = MagneticFields[['x', 'y', 'z', 'B_Gy']]
    GradYdata = GradYdata.rename(
        columns={"B_Gy": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order[1], r_outer=150, scale=scale[1])

    # G_z
    GradZdata = MagneticFields[['x', 'y', 'z', 'B_Gz']]
    GradZdata = GradZdata.rename(
        columns={"B_Gz": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order[2], r_outer=150, scale=scale[2])

    # B0
    try:
        B0_data = MagneticFields[['x', 'y', 'z', 'B0']]
        B0_data = B0_data.rename(
            columns={"B0": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
        B0_Harmonics = SphericalHarmonicFit(B0_data, n_order=n_order[3], r_outer=150, scale=scale[3])
    except KeyError:
        B0_Harmonics = None

    return G_x_Harmonics, G_y_Harmonics, G_z_Harmonics, B0_Harmonics
