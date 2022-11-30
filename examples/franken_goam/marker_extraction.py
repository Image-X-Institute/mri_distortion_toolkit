from pathlib import Path
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import enumerate_subfolders

dataloc = Path(r'/home/brendan/Downloads/FrankenGoam^Mr/20221107 MR Linac^Test')

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
         '13': '14 t1_tse_256_sag_FH',
         '14': '15 t1_tse_256_cor_RL',
         '15': '16 t1_tse_256_cor_LR',
         '16': '17 localiser_gre',
         '17': '18 t1_tse_256_sag_HF_rot',
         '18': '19 t1_tse_256_sag_FH_rot',
         '19': '20 trufi_sag_128_torsocoil',
         '20': '21 trufi_sag_128_torsocoil',
         '21': '22 trufi_sag_128_torsocoil'}

# process TSE images
scans_to_segment = ['9', '11', '12', '13', '14', '15', '18']
for scan in scans_to_segment:
    volume = MarkerVolume(dataloc / scans[scan] / 'Original', n_markers_expected=618, iterative_segmentation=True)
    print(f'for {scans[scan]}, {volume.MarkerCentroids.shape[0]} markers found')
    volume.export_to_slicer()
    volume.save_dicom_data()

# test marker matching
# volume_forward = MarkerVolume(dataloc / scans['9'] / 'Original' / 'slicer_centroids.mrk.json')
# volume_back = MarkerVolume(dataloc / scans['11'] / 'Original' / 'slicer_centroids.mrk.json')
# matched_vol = MatchedMarkerVolumes(volume_back, volume_forward, n_refernce_markers=9)