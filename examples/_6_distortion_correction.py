from mri_distortion_toolkit.DistortionCorrection import KspaceDistortionCorrector
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from pathlib import Path

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''
distorted_data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data\MR\04 gre_trans_AP_330')
dis_volume = MarkerVolume(distorted_data_loc)


GDC = KspaceDistortionCorrector(ImageDirectory=distorted_data_loc.resolve(),
                                gradient_harmonics=[Path('_example_data/G_x_Harmonics.csv').resolve(),
                                                    Path('_example_data/G_y_Harmonics.csv').resolve(),
                                                    Path('_example_data/G_z_Harmonics.csv').resolve()],
                                ImExtension='dcm',
                                dicom_data=dis_volume.dicom_data,
                                correct_through_plane=False)

GDC.correct_all_images()
GDC.save_all_images()  # saves as png so you can quickly inspect results
GDC.save_all_images_as_dicom()  # saves as dicom which can be read into analysis packages.