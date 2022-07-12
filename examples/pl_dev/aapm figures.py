import sys
from pathlib import Path
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes


# Data was acquired on 28/04/2022
data_loc_CT = Path(r'..\_example_data').resolve()

distorted_volume = MarkerVolume(data_loc_CT / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json', verbose=False)
# distorted_volume.plot_3D_markers(title='Otsu')

data_loc_CT = Path(r'..\_example_data').resolve()
ground_truth_volume = MarkerVolume(data_loc_CT / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)
# ground_truth_volume.plot_3D_markers(colour='#1f77b4')

# matched volumes
matched_volume = MatchedMarkerVolumes(ground_truth_volume, distorted_volume, ReverseGradientData=distorted_volume,
                                      ReferenceMarkers=11)
matched_volume.plot_3D_markers()



