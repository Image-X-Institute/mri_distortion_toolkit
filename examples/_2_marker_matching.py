import sys
from pathlib import Path
sys.path.insert(0, 'C:/Users/Brendan/Documents/python/MRI_DistCorrectionPhantom')
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes

this_file_loc = Path(__file__).parent.resolve()
data_loc = this_file_loc / '_example_data'

# distorted centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json', verbose=False)
distorted_volume_rev = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330' / 'slicer_centroids.mrk.json', verbose=False)

# ground truth centroids
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)

# matched volumes
matched_volume = MatchedMarkerVolumes(ground_truth_volume, distorted_volume, ReverseGradientData=distorted_volume_rev,
                                      ReferenceMarkers=11)
matched_volume.MatchedCentroids.to_csv(data_loc / 'Matched_Markers.csv')  # for use in later examples
