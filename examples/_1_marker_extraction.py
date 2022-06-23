from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from pathlib import Path

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''

data_loc = Path('/home/brendan/Downloads/MRI_distortion_QA_sample_data')
# ^^ update to where you put the sample data!!
this_file_loc = Path(__file__).parent.resolve()

gt_data_loc = data_loc / 'CT'
mr_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330'
mr_data_loc_reverse_gradient = data_loc / 'MR' / '05 gre_trans_PA_330'

gt_volume = MarkerVolume(gt_data_loc, r_max=300)
gt_volume.plot_3D_markers()  # produce a quick plot of marker positions
gt_volume.export_to_slicer(save_path=this_file_loc / '_example_data' / 'CT')

mr_volume = MarkerVolume(mr_data_loc, correct_fat_water_shift=True, fat_shift_direction=-1)
mr_volume.save_dicom_data(save_path=this_file_loc / '_example_data' / 'MR' / '04 gre_trans_AP_330')
# save necessary acquisition data as json for easy use later (only works for MR data)
mr_volume.export_to_slicer(save_path=this_file_loc / '_example_data' / 'MR' / '04 gre_trans_AP_330')

mr_volume_rev = MarkerVolume(mr_data_loc_reverse_gradient, correct_fat_water_shift=True, fat_shift_direction=1)
mr_volume_rev.save_dicom_data(save_path=this_file_loc / '_example_data' / 'MR' / '05 gre_trans_PA_330')
mr_volume_rev.export_to_slicer(save_path=this_file_loc / '_example_data' / 'MR' / '05 gre_trans_PA_330')
