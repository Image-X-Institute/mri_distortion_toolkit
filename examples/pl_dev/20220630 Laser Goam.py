import sys
from pathlib import Path

import MRI_DistortionQA.utilities
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA import utilities
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
import matplotlib.pyplot as plt

# Ground truth CT
data_loc_CT = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                   r'Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam').resolve()
ground_truth_volume = MarkerVolume(data_loc_CT / 'CT' / 'slicer_centroids_fixed_cross.mrk.json')
# ground_truth_volume = MarkerVolume(data_loc_CT / 'CT', verbose=True, n_markers_expected=528)
# ground_truth_volume.export_to_slicer()
ground_truth_volume.plot_3D_markers(title='CT')


# Data was acquired on 30/06/2022
# There should be 531 markers, but some of the ones near the vaseline are being missed as they wer too close
data_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney ('
                r'Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\MR').resolve()

# 04 gre_trans_AP_330
# distorted_volume = MarkerVolume(data_loc / '04 gre_trans_AP_330' / 'Original', verbose=True,
#                                 iterative_segmentation=True, n_markers_expected=529, gaussian_image_filter_sd=1)
# # distorted_volume = MarkerVolume(data_loc / '04 gre_trans_AP_330' / 'Original' / 'slicer_centroids.mrk.json')
# distorted_volume.plot_3D_markers(title='04 gre_trans_AP_330')
# # distorted_volume.export_to_slicer()

# 05 gre_trans_PA_330
# distorted_volume = MarkerVolume(data_loc / '05 gre_trans_PA_330' / 'Original', verbose=True,
#                                 iterative_segmentation=True, n_markers_expected=529, gaussian_image_filter_sd=1)
# # distorted_volume = MarkerVolume(data_loc / '04 gre_trans_AP_330' / 'Original' / 'slicer_centroids.mrk.json')
# # distorted_volume.plot_3D_markers(title='04 gre_trans_AP_330')
# distorted_volume.export_to_slicer()

# # 08 gre_cor_RL_330
# distorted_volume = MarkerVolume(data_loc / '08 gre_cor_RL_330' / 'Original', verbose=True,
#                                 iterative_segmentation=True, n_markers_expected=529, gaussian_image_filter_sd=1)
distorted_volume = MarkerVolume(data_loc / '08 gre_cor_RL_330' / 'Original' / 'slicer_centroids1_cross.mrk.json')
distorted_volume.plot_3D_markers(title='08 gre_cor_RL_330')
# distorted_volume.export_to_slicer()

# 09 gre_cor_LR_330
# distorted_volume_rev = MarkerVolume(data_loc / '09 gre_cor_LR_330' / 'Original', verbose=True,
#                                 iterative_segmentation=True, n_markers_expected=529, gaussian_image_filter_sd=1)
distorted_volume_rev = MarkerVolume(data_loc / '09 gre_cor_LR_330' / 'Original' / 'slicer_centroids1_cross.mrk.json')
distorted_volume_rev.plot_3D_markers(title='09 gre_cor_LR_330')
# distorted_volume_rev.export_to_slicer()


#Match
matched_volume = MatchedMarkerVolumes(ground_truth_volume, distorted_volume_rev,
                                      ReferenceMarkers=6)
matched_volume.plot_3D_markers(add_arrows=True)
matched_volume.export_to_csv(save_path=data_loc / '08 gre_cor_RL_330' / 'Original')
