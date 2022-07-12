import sys
from pathlib import Path

import MRI_DistortionQA.utilities
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA import utilities
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
import matplotlib.pyplot as plt

# Data was acquired on 28/04/2022
data_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                r'Staff)\Documents\MRI-Linac\Experiments\20220413 Full sphere foam Goam\20220428 MR '
                r'Linac^Test').resolve()

# Data was acquired on 02/06/2022
data_loc_distorted = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                          r'Staff)\Documents\MRI-Linac\Experiments\20220413 Full sphere foam Goam\20220428 MR '
                          r'Linac^Test\12 gre_sag_HF_corrected').resolve()

data_loc_CT = Path(r'..\_example_data').resolve()
ground_truth_volume = MarkerVolume(data_loc_CT / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)

# ####################
# # Data 1, Scan 12
# # Original (Otsu threshold)
# distorted_volume_12 = MarkerVolume(data_loc / '12 gre_sag_HF' / 'Original', verbose=True)
# #distorted_volume_12.plot_3D_markers(title='Otsu')
#
# Precise threshold
distorted_volume_12_precise = MarkerVolume(data_loc / '12 gre_sag_HF' / 'Original', verbose=True, iterative_segmentation=True, n_markers_expected=336)
# will still fail at other guassian values
distorted_volume_12_precise.plot_3D_markers(title='Original 12 gre_sag_HF')

####################
# Data 1, Scan 12
# Distortion corrected

# Precise threshold
distorted_volume_12_dc = MarkerVolume(data_loc_distorted, verbose=True, iterative_segmentation=True,
                                           n_markers_expected=336, gaussian_image_filter_sd=1.2)
distorted_volume_12_dc.plot_3D_markers(title='Distortion corrected 12 gre_sag_HF')

MRI_DistortionQA.utilities.plot_MarkerVolume_overlay([distorted_volume_12_precise, distorted_volume_12_dc])

plt.imshow(distorted_volume_12_dc.ThresholdVolume[:, :, 16])

# matched_volume_precise = MatchedMarkerVolumes(ground_truth_volume, distorted_volume_12_precise, ReferenceMarkers=11)
# matched_volume_precise.plot_3D_markers(title='Precise')
