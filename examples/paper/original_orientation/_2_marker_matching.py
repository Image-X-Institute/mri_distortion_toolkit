from pathlib import Path
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import enumerate_subfolders

dataloc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/20221107 MR Linac^Test')
gt_data_loc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/CT/slicer_centroids.mrk.json')
scans = {'0': '01 localiser_gre',
         '1': '02 localiser_gre',
         '2': '03 localiser_gre',
         '3': '04 gre_trans_AP_330',
         '4': '05 gre_trans_PA_330',
         '5': '06 gre_sag_AP_330',
         '6': '07 gre_sag_PA_330',
         '7': '08 gre_cor_RL_330',
         '8': '09 gre_cor_LR_330',
         '9': '10 t1_tse_256_sag',
         '10': '11 t1_tse_256_sag_PA',
         '11': '12 t1_tse_256_tra_PA',
         '12': '13 t1_tse_256_sag_HF',
         '13': '14 t1_tse_256_sag_FH',  # for this one had to add gaussian_image_filter_sd=0.8
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR',
         '16': '17 localiser_gre',
         '17': '18 t1_tse_256_sag_HF_rot',
         '18': '19 t1_tse_256_sag_FH_rot',
         '19': '20 trufi_sag_128_torsocoil',
         '20': '21 trufi_sag_128_torsocoil',
         '21': '22 trufi_sag_128_torsocoil'}



gt_vol = MarkerVolume(gt_data_loc)
gt_vol.rotate_markers(yaxis_angle=180)
tra_forward_vol = MarkerVolume(dataloc / scans['9'] / 'Original' / 'slicer_centroids.mrk.json')
tra_back_vol = MarkerVolume(dataloc / scans['11'] / 'Original' / 'slicer_centroids.mrk.json')
tra_match = MatchedMarkerVolumes(gt_vol, tra_forward_vol, reverse_gradient_data=tra_back_vol, n_refernce_markers=9)
tra_match.MatchedCentroids.to_csv('_data/tra_markers.csv')

sag_forward_vol = MarkerVolume(dataloc / scans['12'] / 'Original' / 'slicer_centroids.mrk.json')
sag_back_vol = MarkerVolume(dataloc / scans['13'] / 'Original' / 'slicer_centroids.mrk.json')
sag_match = MatchedMarkerVolumes(gt_vol, sag_forward_vol, reverse_gradient_data=sag_back_vol, n_refernce_markers=9)
sag_match.MatchedCentroids.to_csv('_data/sag_markers.csv')

cor_forward_vol = MarkerVolume(dataloc / scans['14'] / 'Original' / 'slicer_centroids.mrk.json')
cor_back_vol = MarkerVolume(dataloc / scans['15'] / 'Original' / 'slicer_centroids.mrk.json')
cor_match = MatchedMarkerVolumes(gt_vol, cor_forward_vol, reverse_gradient_data=cor_back_vol, n_refernce_markers=9)
cor_match.MatchedCentroids.to_csv('_data/cor_markers.csv')