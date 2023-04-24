from pathlib import Path
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import enumerate_subfolders
import matplotlib as mpi
mpi.rcParams['figure.dpi'] = 200

dataloc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\20221107 MR Linac^Test')
ct_loc = Path(r'C:\Users\bwhe3635\Downloads\FrankenGoam^Mr\FrankenGoam^Mr\CT\slicer_centroids.mrk.json')

scans = {'9': '10 t1_tse_256_sag',
         '11': '12 t1_tse_256_tra_PA',
         '12': '13 t1_tse_256_sag_HF',
         '13': '14 t1_tse_256_sag_FH',  # for this one had to add gaussian_image_filter_sd=0.8
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR',
         '17': '18 t1_tse_256_sag_HF_rot',
         '18': '19 t1_tse_256_sag_FH_rot',
         }

gt_vol = MarkerVolume(ct_loc)
gt_vol.rotate_markers(yaxis_angle=90)
gt_vol.rotate_markers(zaxis_angle=180)
gt_vol.rotate_markers(xaxis_angle=180)
tra_forward_vol = MarkerVolume(dataloc / scans['17'] / 'Original' / 'slicer_centroids.mrk.json')
tra_back_vol = MarkerVolume(dataloc / scans['18'] / 'Original' / 'slicer_centroids.mrk.json')
tra_match = MatchedMarkerVolumes(gt_vol, tra_forward_vol, reverse_gradient_data=tra_back_vol, n_refernce_markers=9)
tra_match.MatchedCentroids.to_csv('_data/rot_markers.csv')

vol_temp = MarkerVolume(tra_match._gt_reference)

vol_temp.rotate_markers(xaxis_angle=-180)
vol_temp.rotate_markers(zaxis_angle=-180)
vol_temp.rotate_markers(yaxis_angle=-90)
vol_temp.export_to_slicer('_data/', 'gt_ref')
# vol_temp = MarkerVolume(tra_match._dist_reference)
# vol_temp.export_to_slicer('_data/', 'dis_ref')