import sys
from pathlib import Path
sys.path.insert(0, 'C:/Users/Brendan/Documents/python/MRI_DistCorrectionPhantom')
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA.MarkerAnalysis import MatchedMarkerVolumes

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data')
# ^^ update to where you put the sample data!!

# distorted centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'MR.mrk.json', verbose=False)

# ground truth centroids
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'MR.mrk.json', verbose=False, r_max=300)

# matched volumes
matched_volume = MatchedMarkerVolumes(ground_truth_volume, distorted_volume, ReferenceMarkers=11)
matched_volume.MatchedCentroids.to_csv('Matched_Markers.csv')  # for use in later examples

# plot the match
matched_volume.plot_3D_markers()

# matched volumes no reference
matched_volume_no_ref = MatchedMarkerVolumes(ground_truth_volume, distorted_volume)

# plot the match
matched_volume_no_ref.plot_3D_markers(title='matched markers, no ref (failure!)')

# distorted centroids, reversed gradient
distorted_volume_rev = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330' / 'MR.mrk.json', verbose=False)

# matched volumes including reversed gradient data
matched_volume_with_rev_data = MatchedMarkerVolumes(ground_truth_volume, distorted_volume,
                                                    ReverseGradientData=distorted_volume_rev, ReferenceMarkers=11)

