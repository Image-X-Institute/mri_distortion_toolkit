from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit.utilities import plot_distortion_xyz_hist
from mri_distortion_toolkit.utilities import get_dicom_data
from pathlib import Path

'''
download example data and unzip:
C:\Users\Brendan\cloudstor\MRI-Linac Experimental Data\UQ PHANTOM^SHIMMING TEST
'''

distorted_data_loc = Path(r'C:\Users\Brendan\Documents\temp\mr_testing\03 t2_tse_tra_fullphantom_z axis long')
try:
    dicom_data = get_dicom_data(distorted_data_loc / 'dicom_data.json')
except FileNotFoundError:
    volume = MarkerVolume(distorted_data_loc, ImExtension='dcm')
    volume.save_dicom_data()
    dicom_data = get_dicom_data(distorted_data_loc / 'dicom_data.json')
GDC = ImageDomainDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                     gradient_harmonics=[Path('../_data/Gx.csv').resolve(),
                                                         Path('../_data/Gy.csv').resolve(),
                                                         Path('../_data/Gz.csv').resolve()],
                                     B0_harmonics=Path('../_data/B0.csv').resolve(),
                                     dicom_data=dicom_data,
                                     ImExtension='dcm',
                                     correct_through_plane=True,
                                     correct_B0=True,
                                     B0_direction='back')

GDC.correct_all_images()
GDC.save_all_images(DSV_radius=150)  # saves as png so you can quickly inspect results
GDC.save_all_images_as_dicom()  # saves as dicom which can be read into analysis packages.
