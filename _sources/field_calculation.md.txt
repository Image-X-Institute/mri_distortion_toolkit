# Field Calculation



Once we have the marker positions created from each field, we can convert these positions into magnetic fields for each of the gradient coils (and, if you used the reverse gradient technique, and estimate of B0). This example shows how to do this.

 coverage-badge -f -o docsrc/__resources/coverage.svgbash

```python
from MRI_DistortionQA.FieldCalculation import ConvertMatchedMarkersToBz
import pandas as pd
from pathlib import Path

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data')

# load the matched volume calculated in the previous step.
matched_volume = pd.read_csv('Matched_Markers.csv', index_col=0).squeeze("columns")
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume using save_dicom_data()

Bz_field = ConvertMatchedMarkersToBz(matched_volume, dicom_data_loc)
Bz_field.MagneticFields.to_csv('Bfields.csv')  # save for later
```

This is really only an intermidiate step; this data isn't particularly useful by itself, but it does allow us to fit spherical harmonics, which is the [next step](https://acrf-image-x-institute.github.io/MRI_DistortionQA/fit_spherical_harmonics.html)!!

