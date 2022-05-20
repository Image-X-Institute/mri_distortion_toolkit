# Reporting

OK, the final step of our journey!

This section demonstrates how you can automatically generate interactive html reports from the data we have generated. There are two ways for you to do this, both of which will be demonstrated:
An example of the types of reports you can generate is [here](_static/MR_QA_report_12_05_2022.html).

1. Pass the reporting code a data frame containg ground truth and distorted marker positions
2. Pass the reporting code spherical harmonics

## Case 1: passing data directly 

Create a new file called 'reporting.py'. Copy the below code into it.

```python
from MRI_DistortionQA.Reports import MRI_QA_Reporter
import pandas as pd
from pathlib import Path

# Direct data case: pass matched marker volume to MRI_QA_Reporter
# ---------------------------------------------------------------
data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data')
dicom_data_loc = data_loc / 'MR' / '04 gre_trans_AP_330' / 'dicom_data.json'  # previosly saved from a MarkerVolume using
Matched_Markers = pd.read_csv('Matched_Markers.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(MatchedMarkerVolume=Matched_Markers, r_outer=150, dicom_data=dicom_data_loc)
report.write_html_report()
```

This code will generate a report at {your_home_directory} / 'Documents' / 'MRI_QA_reports'.

To be honest, this report looks pretty bad - this is because with this phantom, there are very limited data points, so most of the plotting routines don't really work. If you have a more conventional phantom with lots of datapoints, this should work  a lot better.

This phantom was actually designed to get a good measurement of data on the surface of a sphere for the purpose of fitting spherical harmonics; therefore, let's move on and use the data we have more appropriately!! 

> :warning: I have occasionally seen firefox fail to display these reports. I haven't figured out what, but if there is an error "The address wasn’t understood" try a different browser.

## Case 2: harmonic reconstruction

the code for harmonic reconstruction is below:

```python
from MRI_DistortionQA.Reports import MRI_QA_Reporter
import pandas as pd
from pathlib import Path

# Harmonic case: pass harmonics to MRI_QA_Reporter so that data can be recontructed
# ----------------------------------------------------------------------------------
G_x_harmonics = pd.read_csv('G_x_harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv('G_y_harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv('G_z_harmonics.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics],
                         r_outer=150, dicom_data=dicom_data_loc)
report.write_html_report()

```

You will now have a new report sitting {your_home_directory} / 'Documents' / 'MRI_QA_reports'.  This one should [look a lot better](https://acrf-image-x-institute.github.io/MRI_DistortionQA/_static/MR_QA_report_12_05_2022.html)!! 

If you complete the B0 estimate parts of the previous tutorials, and have a 'B0_harmonics.csv' file sitting in your working directory, you can also add this to the call to include a plot of B0 homogeneity:

```python
report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics],
                         r_outer=150, dicom_data=dicom_data_loc, B0_harmonics='B0_harmonics.csv')
report.write_html_report()
```



## Adding custom tests

You will notice that some tests have been run (and failed) from 'DefaultTestSuite'. What is that and how do you add your own tests?

Code demonstration the creation of a custom test suite is below:

```python
from MRI_DistortionQA.Reports import MRI_QA_Reporter
import pandas as pd
from pathlib import Path


class CustomTestSuite:

    def test_case_1(self):
        # a test can return a bool:
        return True

    def test_case_2(self):
        # or a test can return a string:
        return "I am a string!"

    def test_case_3(self):
        # tests have access to the test data:
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=100)
        if test_data.abs_dis.max() < 2:
            return True
        else:
            return False


# Harmonic case: pass harmonics to MRI_QA_Reporter so that data can be recontructed
# ----------------------------------------------------------------------------------
G_x_harmonics = pd.read_csv('G_x_harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv('G_y_harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv('G_z_harmonics.csv', index_col=0).squeeze("columns")

report = MRI_QA_Reporter(gradient_harmonics=[G_x_harmonics, G_y_harmonics, G_z_harmonics],
                         r_outer=150, dicom_data=dicom_data_loc,
                         tests_to_run=CustomTestSuite)  # custom test class passed to tests_to_run
report.write_html_report()

G_x_harmonics = pd.read_csv('G_x_harmonics.csv', index_col=0).squeeze("columns")
G_y_harmonics = pd.read_csv('G_y_harmonics.csv', index_col=0).squeeze("columns")
G_z_harmonics = pd.read_csv('G_z_harmonics.csv', index_col=0).squeeze("columns")
```



