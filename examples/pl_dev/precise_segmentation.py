import sys
from pathlib import Path
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes

# Data was acquired on 28/04/2022
data_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                r'Staff)\Documents\MRI-Linac\Experiments\20220413 Full sphere foam Goam\20220428 MR '
                r'Linac^Test').resolve()

# Data was acquired on 02/06/2022
data_loc2 = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                 r'Staff)\Documents\MRI-Linac\Experiments\20220413 Full sphere foam Goam\20220602 MR '
                 r'Linac^Test').resolve()

data_loc_CT = Path(r'..\_example_data').resolve()
ground_truth_volume = MarkerVolume(data_loc_CT / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)

# ####################
# # Data 1, Scan 9
# # Original (Otsu threshold)
# distorted_volume_12 = MarkerVolume(data_loc / '09 gre_trans_AP_330' / 'Original', verbose=True)
# #distorted_volume_12.plot_3D_markers(title='Otsu')
#
# # Precise threshold
# distorted_volume_12_precise = MarkerVolume(data_loc / '09 gre_trans_AP_330' / 'Original', verbose=True, precise_segmentation=True, n_markers_expected=336, gaussian_image_filter_sd=0.75)
# # will still fail at other guassian values
# distorted_volume_12_precise.plot_3D_markers(title='Precise')
#
# matched_volume_precise = MatchedMarkerVolumes(ground_truth_volume, distorted_volume_12_precise, ReferenceMarkers=11)
# matched_volume_precise.plot_3D_markers(title='Precise')


# ####################
# # Data 1, Scan 12
# # Original (Otsu threshold)
# distorted_volume_12 = MarkerVolume(data_loc / '12 gre_sag_HF' / 'Original', verbose=True)
# #distorted_volume_12.plot_3D_markers(title='Otsu')
#
# # Precise threshold
# distorted_volume_12_precise = MarkerVolume(data_loc / '12 gre_sag_HF' / 'Original', verbose=True, precise_segmentation=True, n_markers_expected=336, gaussian_image_filter_sd=0.75)
# # will still fail at other guassian values
# distorted_volume_12_precise.plot_3D_markers(title='Precise')
#
# matched_volume_precise = MatchedMarkerVolumes(ground_truth_volume, distorted_volume_12_precise, ReferenceMarkers=11)
# matched_volume_precise.plot_3D_markers(title='Precise')

# ####################
# # Data 1, Scan 14
# # Original (Otsu threshold)
# distorted_volume_14 = MarkerVolume(data_loc / '14 gre_cor_RL' / 'Original', verbose=True)
# distorted_volume_14.plot_3D_markers(title='Otsu')
#
# # Precise threshold
# distorted_volume_14_precise = MarkerVolume(data_loc / '14 gre_cor_RL' / 'Original', verbose=True, precise_segmentation=True, n_markers_expected=336, gaussian_image_filter_sd=0.7)
# # will still fail at other guassian values
# distorted_volume_14_precise.plot_3D_markers(title='Precise')
#
# matched_volume_precise = MatchedMarkerVolumes(ground_truth_volume, distorted_volume_14_precise, ReferenceMarkers=11)
# matched_volume_precise.plot_3D_markers(title='Precise match')

####################
# Data 2, Scan 2
# Original (Otsu threshold)
distorted_volume_2 = MarkerVolume(data_loc2 / '02 gre_tra_AP' / 'Original', verbose=True, gaussian_image_filter_sd=0.7)
distorted_volume_2.plot_3D_markers(title='Otsu')

# Precise threshold
distorted_volume_2_precise = MarkerVolume(data_loc2 / '02 gre_tra_AP' / 'Original', verbose=True, precise_segmentation=True, n_markers_expected=145, gaussian_image_filter_sd=0.75)
# will still fail at other guassian values
distorted_volume_2_precise.plot_3D_markers(title='Precise')
