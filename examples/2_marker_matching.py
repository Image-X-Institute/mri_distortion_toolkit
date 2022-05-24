import sys
from pathlib import Path
sys.path.insert(0, 'C:/Users/Brendan/Documents/python/MRI_DistCorrectionPhantom')
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes


data_loc = Path(r'_example_data').resolve()

# distorted centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json', verbose=False)
distorted_volume_rev = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330' / 'slicer_centroids.mrk.json', verbose=False)
# ^ data taken with reverse phase and frequency encoding

# ground truth centroids
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)

# matched volumes
matched_volume = MatchedMarkerVolumes(ground_truth_volume, distorted_volume, ReverseGradientData=distorted_volume_rev,
                                      ReferenceMarkers=11)
matched_volume.MatchedCentroids.to_csv('_example_data/Matched_Markers.csv')  # for use in later examples

# plot the match
matched_volume.plot_3D_markers()

# matched volumes no reference
matched_volume_no_ref = MatchedMarkerVolumes(ground_truth_volume, distorted_volume)

# plot the match
matched_volume_no_ref.plot_3D_markers(title='matched markers, no ref (failure!)')

# distorted centroids, reversed gradient
distorted_volume_rev = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330' / 'slicer_centroids.mrk.json', verbose=False)

# matched volumes including reversed gradient data
matched_volume_with_rev_data = MatchedMarkerVolumes(ground_truth_volume, distorted_volume,
                                                    ReverseGradientData=distorted_volume_rev, ReferenceMarkers=11)

