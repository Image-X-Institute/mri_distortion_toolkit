

from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from pathlib import Path
import numpy as np
import pandas as pd

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''

data_loc = Path('C:\\Users\\Brendan\\Downloads\\MRI_distortion_QA_sample_data(2)\\MRI_distortion_QA_sample_data\\')
# ^^ update to where you put the sample data!!
marker_volume = MarkerVolume(data_loc / 'CT', verbose=True)
marker_volume.plot_3D_markers()  # produce a quick plot of marker positions

# # pandas data frame read in
# r_outer = 150
# test_data = np.random.rand(100, 3) * r_outer  # create some random data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

# test_data = pd.DataFrame(test_data, columns=['x', 'y', 'z'])  # convert to data frame
# pandas_volume = MarkerVolume(test_data)  # create MarkerVolume
#
# # json read in
# json_file = data_loc / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json'
# json_volume = MarkerVolume(json_file)  # create MarkerVolume