# Marker Extraction

## Basic examples

Say that you have [built](https://acrf-image-x-institute.github.io/MRI_DistortionPhantom/phantom_construction.html) and [imaged](https://acrf-image-x-institute.github.io/MRI_DistortionPhantom/phantom_imaging.html) a marker-based distortion phantom. To use this data within this library, you first have to extract the position of the markers and create a 'MarkerVolume'. This example shows you how do that.


> :warning: For this part you will need some data. **Example data is provided [here](https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU)**. Download and unzip this data somewhere and take note of the path.

First, create a directory called 'MRI_QA_tutorial' or something like that. Within that directory, create a new python file called 'MarkerExtractionExample'. Copy the below code into it, and update

```python
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from pathlib import Path
import numpy as np
import pandas as pd

'''
download example data and unzip:
https://cloudstor.aarnet.edu.au/plus/s/Wm9vndV47u941JU
'''

data_loc = Path(r'C:\Users\Brendan\Downloads\MRI_distortion_QA_sample_data\MRI_distortion_QA_sample_data')
# ^^ update to where you put the sample data!!
marker_volume = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330', verbose=False)
marker_volume.export_to_slicer()  # save this data as json for easy use later
marker_volume.save_dicom_data()  # save this data as json for easy use later
marker_volume.plot_3D_markers()  # produce a quick plot of marker positions
```

This is the simplest interface to creating an instance of 'MarkerVolume' from a series of dicom slices. If you want to get at the underlying data, it is stored in ```marker_volume.MarkerCentroids``` as a pandas dataframe. The MarkerVolume object is one of the building blocks of this code. As such, there are multiple ways you can instantiate it:

```python
# pandas data frame read in
r_outer = 150
test_data = np.random.rand(100, 3) * r_outer  # create some random data
test_data = pd.DataFrame(test_data, columns=['x', 'y', 'z'])  # convert to data frame
pandas_volume = MarkerVolume(test_data)  # create MarkerVolume

# json read in
json_file = data_loc / 'MR' / '04 gre_trans_AP_330' / 'MR.mrk.json'
# this file was created with the export_to_slicer step above; this is also the format used by slicer
json_volume = MarkerVolume(json_file)  # create MarkerVolume
```

At this point, you are free to move on to the next step: automatic matching of markers, or you can read on to find out a little bit more about how you can handle marker extraction in challenging datasets...

## Do we guarantee that we will find your markers?

Short answer: No! 

Although the automatic extraction works quite well in most cases, because there are so many variables in MR, we have no knowledge of the signal-to-noise, contrast-to-noise, contrast type, voxel size, etc. that you may be using. This means that it is very difficult to automatically know what settings to use for marker extraction. In some low SNR cases, no matter what settings you use automatic extraction is difficult, but in most cases you should be able to find a reliable combination of settings for a given scan and scanner.

## Is that a major issue?

Also no!

We provide an easy interface to [slicer](https://www.slicer.org/) via the ```export_to_slicer``` method; we also read these slicer .json files back in as demonstrated in the example above. This means that in situations where the automatic marker processing fails, you are free to move, delete and add markers through the excellent slicer GUI. Once you are satisfied, you can go file>>save data and save the *.mrk.json file for reading back into this workflow. A screenshot of the process of editing marker locations in slicer is below:

![](__resources/Slicer_Markers_screengrab.PNG)

## Handling Fat-water shift 

If you are using a phantom with oil filled markers, your images may be subject to [fat-water shift](https://acrf-image-x-institute.github.io/MRI_DistortionPhantom/phantom_imaging.html#fat-water-chemical-shift).

The best way to check this is to take a forward/reverse gradient pair of images, and compare the markers in the middle of the DSV. Since B0 homogeneity is very good here, if the markers are offset from each other it is due to fat/water shift.

You can read about the different options for handling this effect [here](https://acrf-image-x-institute.github.io/MRI_DistortionPhantom/phantom_imaging.html#fat-water-chemical-shift); but one option is to correct for this in software. If you want to do that, the code would look this:

```python
from MRI_DistortionQA.MarkerAnalysis import MarkerVolume
from pathlib import Path
from MRI_DistortionQA.utilities import plot_MarkerVolume_overlay

marker_volume_forward = MarkerVolume(data_loc / 'MR' / '04 gre_trans_AP_330', verbose=False,
                             correct_fat_water_shift=True, fat_shift_direction=-1)
marker_volume_back = MarkerVolume(data_loc / 'MR' / '05 gre_trans_PA_330', verbose=False,
                             correct_fat_water_shift=True, fat_shift_direction=1)
plot_MarkerVolume_overlay([marker_volume_forward, marker_volume_back])
```

- This will apply a shift to the marker positions based on the estimate of fat/water chemical shift
- We know which axis the shift will occur in (the frequency encode direction) but we are not yet confident we can predict the direction (forward/back). ```fat_shift_direction``` controls this. Basically you need to compare the markers in the center of the phantom for the forward/ reverse gradient images. If it worked, you should see that markers in the center of the phantom are closely aligned. If it moved them further apart, change the sign. If they still aren't, then log an issue!
- If you use this feature, please let us know because as you can tell it is still under development a bit!



