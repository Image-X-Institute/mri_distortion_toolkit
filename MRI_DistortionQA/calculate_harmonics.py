from .MarkerAnalysis import MatchedMarkerVolumes
from .FieldCalculation import ConvertMatchedMarkersToBz
from .FieldAnalysis import SphericalHarmonicFit
from .utilities import get_dicom_data
import numpy as np


def calculate_harmonics(ground_truth_volume, distorted_volume, distorted_volume_rev=None, dicom_data=None, n_order=8):
    """
    this function is essentially a wrapper of convenience. Given a ground truth MarkerVolume 
    and a forward (distorted) MarkerVolume, it will calculate and return the Gx, Gy, and Gz
    harmonics. Optionally, you can also pass a back MarkerVolume (reverse gradient) in which case
    it will also calcualte the B0 harmonics. Basically this makes your life easier, but you give up fine control over
    the marker matching, B calc, and harmonic fitting steps.
    Note that the gradient harmonics will be normalised to a gradient strengt of 1 mT/m
    
    :param ground_truth_volume: ground truth MarkerVolume
    :type ground_truth_volume: object
    :param distorted_volume: distorted MarkerVolume
    :type distorted_volume: object
    :param distorted_volume_rev: distorted MarkerVolume taken with reverse gradient
    :type distorted_volume_rev: object, optional
    :param n_order: order of harmonic fit
    :type n_order: int, optional
    :return: B0_Harmonics, G_x_Harmonics, G_y_Harmonics, G_z_Harmonics
    """
    if dicom_data is None:
        try:
            dicom_data = distorted_volume.dicom_data
        except AttributeError:
            raise AttributeError('please supply dicom data')
            return
    dicom_data = get_dicom_data(dicom_data)
    gradient_strength = np.array(dicom_data['gradient_strength']) * 1e3
    # match the markers
    matched_markers = MatchedMarkerVolumes(ground_truth_volume, distorted_volume, sorting_method='radial', ReferenceMarkers=11,
                                              WarpSearchData=True, ReverseGradientData=distorted_volume_rev)

    # calculate B fields
    B_fields = ConvertMatchedMarkersToBz(matched_markers.MatchedCentroids, dicom_data)

    # calculate harmonics
    # B0
    try:
        B0_data = B_fields.MagneticFields[['x', 'y', 'z', 'B0']]
        B0_data = B0_data.rename(
            columns={"B0": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
        B0_Harmonics = SphericalHarmonicFit(B0_data, n_order=n_order, r_outer=150)
    except KeyError:
        B0_Harmonics = None
    # Gx
    GradXdata = B_fields.MagneticFields[['x', 'y', 'z', 'B_Gx']]
    GradXdata = GradXdata.rename(
        columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[0])

    # Gy
    GradYdata = B_fields.MagneticFields[['x', 'y', 'z', 'B_Gy']]
    GradYdata = GradYdata.rename(
        columns={"B_Gy": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[1])

    # G_z
    GradZdata = B_fields.MagneticFields[['x', 'y', 'z', 'B_Gz']]
    GradZdata = GradZdata.rename(
        columns={"B_Gz": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
    G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, r_outer=150, scale=1/gradient_strength[2])

    return B0_Harmonics, G_x_Harmonics, G_y_Harmonics, G_z_Harmonics