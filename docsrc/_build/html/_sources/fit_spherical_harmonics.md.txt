# Spherical Harmonics

## Simple example

Create a new file 'fit_harmonics.py' and copy the below code into it.

```python
from MRI_DistortionQA.FieldAnalysis import SphericalHarmonicFit
import pandas as pd

FieldData = pd.read_csv('Bfields.csv', index_col=0).squeeze("columns")
# load previously saved data

'''
This data contains columns ['x', 'y', 'z', 'B_Gx', 'B_Gy', 'B_Gz']
but the spherical harmonics code expects to receieve [x, y, z, Bz]
therefore, we will need to create a new dataframe with appropriately named columns
for each field we want to fit to:
'''

n_order = 8
# G_x Harmonics
GradXdata = FieldData[['x', 'y', 'z', 'B_Gx']]
GradXdata = GradXdata.rename(columns={"B_Gx": "Bz"})  # spherical harmonics code expects to receieve one field called Bz
G_x_Harmonics = SphericalHarmonicFit(GradXdata, n_order=n_order, r_outer=150)
G_x_Harmonics.harmonics.to_csv('G_x_harmonics.csv')

# some plotting examples
G_x_Harmonics.plot_cut_planes()
G_x_Harmonics.plot_harmonics_pk_pk(cut_off=.01)
G_x_Harmonics.print_key_harmonics(cut_off=.01)
```

![](__resources/x_gradient_cut_planes.png)

**Reconstructed fields in each cardinal plane for the X gradient coil. Note that there is strong variantion in X (as expected) and the field is close to 0 in the ZY plane (as expected)**

![](__resources/x_gradient_harmonics_bar.png)

**This figure shows the dominant harmonics for the X gradient. If you are a harmonics nerd, you will know that the A11 harmonic corresponds to a perfect X gradient field; therefore it is gratifying to see that this is by far the most strongly expressed harmonic for the X gradient!**

## Explaining the code output:

The code probably printed the following to the screen:

```
[FieldAnalysis.py: line 86  WARNING] input sample points do not appear to cover a full sphere
Initial pk-pk:        693.401 μT
Reconstructed pk-pk:  693.795 μT
Residual pk-pk:       2.023 μT
```

The warning here is telling us that the sample points do not appear to cover a full sphere. We can ignore this in situations where we are confident that we have sufficient sampling of points for the order of harmonics we are fitting. 

- [ ] ToDo: automate this check!!

The second part is telling us the peak-to-peak perturbation over the surface of r_outer (150 mm in this case). We would like to see that the reconstructed pk-pk closely matches the input, and that the residual pk-pk is low relative to the total. In this case, the reconstructed pk-pk is within 0.4 μT and the residual is < 1%, so the fit is pretty good!

## Remaining harmonics

Now we have the X-harmonics; we need to do the same thing for the other two gradient coils:

````python
# G_y Harmonics
GradYdata = FieldData[['x', 'y', 'z', 'B_Gy']]
GradYdata = GradYdata.rename(columns={"B_Gy": "Bz"})
G_y_Harmonics = SphericalHarmonicFit(GradYdata, n_order=n_order, r_outer=150)
G_y_Harmonics.harmonics.to_csv('G_y_harmonics.csv')

# G_z Harmonics
GradZdata = FieldData[['x', 'y', 'z', 'B_Gz']]
GradZdata = GradZdata.rename(columns={"B_Gz": "Bz"})
G_z_Harmonics = SphericalHarmonicFit(GradZdata, n_order=n_order, r_outer=150)
G_z_Harmonics.harmonics.to_csv('G_z_harmonics.csv')
````

## The easy way...

Calculating harmonics from marker volumes involves three steps:

1. Matching the volumes
2. Calculating fields from the markers
3. Calculating harmonics from the fields

The way we just showed you gives you a lot of fine control over every step. However, if you are willing to give up this control, we have written a wrapper function that allows you to do this all in one step:

```python
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from MRI_DistortionQA import calculate_harmonics
from pathlib import Path

# download example data and unzip:
# https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data(2)\MRI_distortion_QA_sample_data')
# ^^ update to where you put the sample data!!

# distorted centroids
distorted_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330' / 'slicer_centroids.mrk.json', verbose=False)
# ground truth centroids
ground_truth_volume = MarkerVolume(data_loc / 'CT' / 'slicer_centroids.mrk.json', verbose=False, r_max=300)
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume

B0_Harmonics, G_x_Harmonics, G_y_Harmonics, G_z_Harmonics = 	calculate_harmonics(ground_truth_volume, distorted_volume, dicom_data=dicom_data_loc)
# note that B0_harmonics is None as we did not provide distorted_volume_rev to calculate_harmonics
```

## Next steps

You are ready to move onto [Reporting](https://acrf-image-x-institute.github.io/MRI_DistortionQA/reporting.html)!
