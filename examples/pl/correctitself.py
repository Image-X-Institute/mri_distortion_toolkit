'''
this demonstrates using an image to generate harmonics, then using those harmonics
to correct the image
'''

from pathlib import Path
import numpy as np
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes
from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
from MRI_DistortionQA import calculate_harmonics
from MRI_DistortionQA.utilities import plot_distortion_xyz_hist
from MRI_DistortionQA.DistortionCorrector import KspaceDistortionCorrector
from MRI_DistortionQA.utilities import plot_matched_volume_hist, print_dict
from MRI_DistortionQA.utilities import plot_MarkerVolume_overlay
from matplotlib import pyplot as plt

plt.rcParams["figure.dpi"] = 150  # for 4k screens

## Correct itself

# Data import
corrected_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\MR\08 gre_cor_RL_330\Original\Corrected_dcm')

# Ground truth
gt_centroid_loc = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\CT\slicer_centroids_fixed_cross.mrk.json')

# extract markers:
gt_volume = MarkerVolume(gt_centroid_loc, r_max=300)
gt_volume.plot_3D_markers(title='CT')

dis_volume = MarkerVolume(corrected_loc.resolve(), n_markers_expected=530, iterative_segmentation=True, r_max=155, gaussian_image_filter_sd=0.4)
dis_volume.plot_3D_markers(title='MR')

# match markers:
matched_volume = MatchedMarkerVolumes(gt_volume, dis_volume, n_refernce_markers=6)
matched_volume.plot_3D_markers()


# Graph
# plot_matched_volume_hist([matched_volume, matched_volume_corrected], ['original', 'corrected'])
plot_distortion_xyz_hist(matched_volume)


# ## Brendan's coefficients from foam prototype
#
# # Data import
# corrected_loc_bw = Path(r'C:\Users\pliu4890\OneDrive - The University of Sydney (Staff)\Documents\MRI-Linac\Experiments\20220630 Laser Goam\MR\04 gre_trans_AP_330\Original\Corrected_dcm')
#
# corrected_volume_bw = MarkerVolume(corrected_loc_bw.resolve(), n_markers_expected=531, iterative_segmentation=True, r_max=155, gaussian_image_filter_sd=0.5, verbose=True)
# corrected_volume_bw.plot_3D_markers(title='MR')
#
# matched_volume_bw = MatchedMarkerVolumes(gt_volume, corrected_volume_bw, ReferenceMarkers=6,Rotation=False)
# matched_volume_bw.plot_3D_markers()
# plot_distortion_xyz_hist(matched_volume_bw)