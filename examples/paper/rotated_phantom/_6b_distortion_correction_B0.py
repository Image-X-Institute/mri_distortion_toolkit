from _data_loc import dataloc, scans
from mri_distortion_toolkit.utilities import get_dicom_data
from mri_distortion_toolkit.DistortionCorrection import ImageDomainDistortionCorrector
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from pathlib import Path

distorted_data_loc = dataloc / scans['17'] / 'Original'  # this is a 'rotated' acquisition
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
GDC.save_all_images()  # saves as png so you can quickly inspect results
GDC.save_all_images_as_dicom()  # saves as dicom which can be read into analysis packages.