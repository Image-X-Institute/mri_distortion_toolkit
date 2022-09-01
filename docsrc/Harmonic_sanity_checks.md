# Harmonic sanity checks

## Intro

This write up concerns [6 pairs of images](https://cloudstor.aarnet.edu.au/plus/apps/files/?dir=/Shared/MRI-Linac%20Experimental%20Data/Goam2%5EMr/20220428%20MR%20Linac%5ETest&fileid=6603039901) taken on the Aus MRI-Linac magnet, summarized below:

| Slice direction | Phase encoding direction |
| --------------- | ------------------------ |
| transverse      | Posterior-Anterior       |
| transverse      | Anterior-Posterior       |
| sagittal        | Head-Feet                |
| sagittal        | Feet-Head                |
| coronal         | Right-Left               |
| coronal         | Left-Right               |

With this data, we can get three independent estimates of B0/ gradient harmonics. The purpose of this analysis is to carry out two sanity checks:

1. the order of the forward/reverse images shouldn't effect the harmonics. If we switch between e.g. the AP and PA sets in the marker matching step, we should still get very similar answers
2. the estimates of the harmonics for all the different slice directions should also be similar.

## Sanity check 1: the order of the forward/reverse images shouldn't effect the harmonics

for this data we have

 'phase_encode_direction': 'y', 'freq_encode_direction': 'x', 'slice_direction': 'z'

### Gx harmonics

|      | AP_PA | PA_AP |
| ---- | ----- | ----- |
| A11  | -721  | -721  |
| A31  | -24   | -24   |
| A51  | 24    | 24    |

### Gy harmonics

|      | AP_PA | PA_AP |
| ---- | ----- | ----- |
| A10  | -9    | -14   |
| B11  | -720  | -720  |
| B31  | -24   | -21   |
| B51  | 23    | 23    |

### Gz harmonics

|      | AP_PA | PA_AP |
| ---- | ----- | ----- |
| A10  | -444  | -444  |
| A30  | -27   | -27   |
| A50  | 10    | 20    |

### B0 harmonics

there are a lot of these so these are the ones with 30% of the dominant harmonic

|      | AP_PA | PA_AP |
| ---- | ----- | ----- |
| B11  | 6     | 6     |
| A20  | 8     | 8     |
| B31  | 9     | 9     |
| B51  | 6     | 6     |
| A60  | -5    | -5    |
| A80  | 7     | 7     |

## Sanity check 2: Different imaging plane directions should give the same harmonics

### B0 harmonics

The AP/PA and RL/LR look quite similar, except for the reversal in sign - which I think is more likely to be an error in the field calculation than in the harmonics.

The HF/HF look completely different - the same harmonics aren't even being excited. Further more, the fit data indicates completely non physical data:

```
Initial pk-pk:        50.709 μT
Reconstructed pk-pk:  134.032 μT
Residual pk-pk:       104.138 μT
```

Ok and a quick look at the volumes tells us why!

![](__resources/HFFH_Fail.png)

This was with warping of data turned off. When I turned it on, I got a much  more sensible looking result and that is reported below.. 

|      | AP_PA | HF_FH | RL_LR |
| ---- | ----- | ----- | ----- |
| B11  | 6     | -6    | -8    |
| A20  | 8     | -8    | -8    |
| B31  | 9     | -9    | -10   |
| B51  | 6     | -7    | -7    |
| A60  | -5    | 3     | 5     |
| A80  | 7     | -10   | -10   |

OK!! not perfect, but reasonably stable. All estimates suggest pk-pk perturbation of ~20 uT

### Gx harmonics

|      | AP_PA       | HF_FH       | RL_LR       |
| ---- | ----------- | ----------- | ----------- |
| A11  | **-721 (100%)** | -465 (100%) | -721 (100%) |
| A31  | **-24 (3.3%)**  | -14 (3.0%)  | -22 (3.1%)  |
| A51  | **24 (3.3%)**   | 15 (3.2%)   | 25 (3.5%)   |

### Gy harmonics

|      | AP_PA       | HF_FH       | RL_LR       |
| ---- | ----------- | ----------- | ----------- |
| A10  | -9 (1.3%)   | **-9 (1.3%)**   | -4 (1.0%)   |
| B11  | -720 (100%) | **-720 (100%)** | -465 (100%) |
| B31  | -24 (3.3%)  | **-24 (3.3%)**  | -16 (3.4%)  |
| B51  | 23 (3.2%)   | **22 (3.1%)**   | 13 (2.8%)   |

### Gz harmonics

|      | AP_PA       | HF_FH       | RL_LR       |
| ---- | ----------- | ----------- | ----------- |
| A10  | -444 (100%) | -689 (100%) | **-689 (100%)** |
| A30  | -27 (6.1%)  | -43 (6.2%)  | **-40 (5.8%)**  |
| A50  | 10 (2.3%)   | 16 (-2.3%)  | **17 (2.5%)**   |



## Conclusion

overall, this worked much better than I was expecting - I would go so far as to say the harmonics are exhibiting remarkable stability!

The main concern is the reversal in sign of the harmonics for B0, so we have to get to the bottom of that. 





## python script

in case I ever want this again:

```python
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from pathlib import Path
from mri_distortion_toolkit.FieldCalculation import ConvertMatchedMarkersToBz
from mri_distortion_toolkit.Harmonics import SphericalHarmonicFit

'''
data is here
https://cloudstor.aarnet.edu.au/plus/apps/files/?dir=/Shared/MRI-Linac%20Experimental%20Data/Goam2%5EMr/20220428%20MR%20Linac%5ETest&fileid=6603039901
'''
data_loc = Path(r'X:\PRJ-RPL\2RESEARCH\2_ProjectData\MRI-Linac\20220428 MR Linac^Test')
all_scans = {'1': '01 localiser_gre',
             '2': '02 gre_trans_AP_330',
             '3': '03 gre_trans_PA_330',
             '4': '04 gre_sag_HF_330',
             '5': '05 gre_sag_FH_330',
             '6': '06 gre_cor_RL_330',
             '7': '07 gre_cor_RL_330',
             '8': '08 gre_trans_AP_330_F_reset',
             '9': '09 gre_trans_AP_330',
             '10': '10 gre_trans_AP_330',
             '11': '11 gre_trans_PA',
             '12': '12 gre_sag_HF',
             '13': '13 gre_sag_FH',
             '14': '14 gre_cor_RL',
             '15': '15 gre_cor_LR',
             '16': '16 gre_tran_AP_large_BW',
             '17': '17 gre_tran_PA_large_BW',
             '18': '18 gre_trans_PA_reshim_refreq',
             '19': '19 gre_trans_AP_reshim_refreq'}

correct_FW = True
ct_volume = MarkerVolume('CT.mrk.json')

# AP/PA Harmonics
forward_volume = MarkerVolume(data_loc / all_scans['14'] / 'Original', gaussian_image_filter_sd=1,
                              n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                              correct_fat_water_shift=correct_FW, fat_shift_direction=-1)
back_volume = MarkerVolume(data_loc / all_scans['15'] / 'Original', gaussian_image_filter_sd=1,
                           n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                           correct_fat_water_shift=correct_FW, fat_shift_direction=1)
B0_Harmonics_AP, G_x_Harmonics_AP, G_y_Harmonics_AP, G_z_Harmonics_AP = calculate_harmonics(forward_volume, back_volume)

# RL/LR Harmonics
forward_volume = MarkerVolume(data_loc / all_scans['14'] / 'Original', gaussian_image_filter_sd=1,
                              n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                              correct_fat_water_shift=correct_FW, fat_shift_direction=-1)
back_volume = MarkerVolume(data_loc / all_scans['15'] / 'Original', gaussian_image_filter_sd=1,
                           n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                           correct_fat_water_shift=correct_FW, fat_shift_direction=1)
B0_Harmonics_RL, G_x_Harmonics_RL, G_y_Harmonics_RL, G_z_Harmonics_RL = calculate_harmonics(forward_volume, back_volume)

# HF/FH Harmonics
forward_volume = MarkerVolume(data_loc / all_scans['12'] / 'Original', gaussian_image_filter_sd=1,
                              n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                              correct_fat_water_shift=correct_FW, fat_shift_direction=-1)
forward_volume.save_dicom_data(save_path=Path(__file__).parent)
forward_volume.export_to_slicer()
back_volume = MarkerVolume(data_loc / all_scans['13'] / 'Original', gaussian_image_filter_sd=1,
                           n_markers_expected=336, cutoff_point=50, verbose=False, r_max=165,
                           correct_fat_water_shift=correct_FW, fat_shift_direction=1)

B0_Harmonics_HF, G_x_Harmonics_HF, G_y_Harmonics_HF, G_z_Harmonics_HF =
    calculate_harmonics(back_volume, forward_volume)
```

