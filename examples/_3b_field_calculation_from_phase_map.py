from mri_distortion_toolkit.FieldCalculation import ConvertPhaseMapsToBz
from pathlib import Path
from mri_distortion_toolkit.utilities import dicom_to_numpy
from mri_distortion_toolkit.utilities import get_all_files
import pydicom
import numpy as np
from mri_distortion_toolkit.utilities import sort_dicoms_by_echo_time

'''
Data location
https://cloudstor.aarnet.edu.au/plus/apps/files/?dir=/MRI-Linac%20Experimental%20Data/UQ%20PHANTOM%5ESHIMMING%20TEST/20210520%20QA%5EQA&fileid=5564960679
'''
echo1_loc = Path(r'/home/brendan/Documents/temp/B0 field map data/GRE_FIELD_MAPPING_POS4_XY_0005')
# sort echo times into differeent folders (you may be able to skip this depending on data format)
data_locs = sort_dicoms_by_echo_time(echo1_loc, file_extension='IMA')
# data_locs = [Path('/home/brendan/Documents/temp/B0 field map data/GRE_FIELD_MAPPING_POS4_XY_0005/echo_10.0'),
             # Path('/home/brendan/Documents/temp/B0 field map data/GRE_FIELD_MAPPING_POS4_XY_0005/echo_17.22')]

FieldData = ConvertPhaseMapsToBz(data_locs[0], data_locs[1], file_extension='IMA')

