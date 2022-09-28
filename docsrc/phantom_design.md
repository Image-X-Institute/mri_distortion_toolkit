# Phantom Design

You can use the phantom design module to create phantom designs, which will be automatically build in [FreeCAD](https://www.freecadweb.org/). Please see [here](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/FreeCADsetup.html#setting-up-freecad) for instructions on installing and setting up FreeCAD. Detailed notes on phantom [construction](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/phantom_construction.html) and [imaging](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/phantom_imaging.html) are also provided; the purpose of this section is to provide some examples of generating phantom designs.

The Phantom design module is based around the concept of a `Slice`. You can stack multiple slices to build full phantom. The following script demonstrates the creation of a simple Slice with the [default parameters](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#module-mri_distortion_toolkit.PhantomDesigner):

```python
from mri_distortion_toolkit import PhantomDesigner
try:
    import FreeCAD
    running_in_free_cad = True
except ImportError:
    running_in_free_cad = False

Slice = PhantomDesigner.PhantomSlice()
if running_in_free_cad:
    Slice.draw_slice()
```

Note that this to create any CAD, this script has to be [executed as a FreeCAD macro](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/FreeCADsetup.html). Otherwise, it simply create a python object that holds the geometric parameters.

You have almost complete freedom to alter the slice shape and size, and also change where the marker positions are. Obviously to build a useful 3D phantom, you will need to stack multiple slices on top of each other. A simple example of building a multi slice phantom (again with mostly default [parameters](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#module-mri_distortion_toolkit.PhantomDesigner)) is below:

```python
from pathlib import Path
import importlib
import numpy as np
import sys
path_to_source = Path(__file__).parent.parent
sys.path.insert(0, str(path_to_source))
from mri_distortion_toolkit import PhantomDesigner
try:
    import FreeCAD
    run_in_free_cad = True
except ImportError:
    run_in_free_cad = False
importlib.reload(PhantomDesigner)

Nslices = 11 # make this an odd number to make sure you have a slize at z=0
SliceZPositions = np.linspace(-150, 150, Nslices)
SliceThickness = np.mean(np.abs(np.diff(SliceZPositions)))

for i, z_pos in enumerate(SliceZPositions):
    Slice = PhantomDesigner.PhantomSlice(slice_thickness=SliceThickness, z_pos=z_pos)
    if run_in_free_cad:
        Slice.draw_slice()
        Slice.add_full_scale_drawing()
```

## Phantom Customization

The baseline model of a distortion phantom is highly customizable. You can change any of the parameters in the FoamProtoypeExample.py file. One of main reasons you may wish to do this is that different scanner have different field of views, so you may wish to make your phantom larger or smaller.
All options for the AddPhantomSlice class are described within the [code docs](), but we provide some additional notes on some of the more common things you may wish to change below:

### Specifying marker locations 

The marker locations are specified on each slice object. We provide two methods to automatically generate marker locations: ```HoleCentroids=cartesian``` will generate a cartesian grid of markers, while ```ROI_polar```  will generate concentric rings of markers. Both wil the ```hole_spacing``` parameter to space out markers. If you specify a ```DSV```,  the ```ROI_polar``` option will ensure good marker coverage over the surface of this sphere, and will provide an etch of the intersection of the DSV on each slice surface so you can tell where the DSV is on each slice.

You can specify a crosshair of markers using the ```ReferenceCrosshairRadius``` option. This will add a crosshair of markers within ```ReferenceCrosshairRadius```. This is a good idea to add to the central slice, as it makes alignment with CT/Ground truth much easier.

Finally, you may not wish to use any of the existing methods for defining marker positions. In that case, you are free to simply specify them as a list: ```HoleCentroids = [[x1,x2,x3],[y1,y2,y3]]```

### Specifying a load region

This phantom consists of a base material that does not give MRI signal, and is then packed with oil capsules, which also don't generate much signal. This can result in the RF coil of the scanner not being properly loaded. To avoid this, it is a good idea to add some load to your phantom. You can specify a region to be cut from the center of each slice using e.g. ```LoadRegion={'shape': 'rectangle', 'width': 100, 'height': 200}``` (see code docs for other options).

In our experience, not much load is required: during development we simple put a container of oil capsules into a zip lock bag. The exact location of the load also shouldn't be especially sensitive, just put it somewhere near the middle. 

### Specifying a DSV

Specifying a Diameter of Spherical Volume (DSV) has two effects

1. the intersection of the DSV with each slice will be etched on the surface of the slice
2. If you specify ```HoleCentroids=ROI_polar``` then the code will ensure good marker coverage over the surface of the DSV sphere. This can be important if you wish to fit spherical harmonics using this data.

### Specifying Guide Rods

This phantom is based on the concept of individual slices which are stacked on top of each other. A number of methods can be envisaged to hold all of these slices together, but internally we have been using nylon guide rods with great success. 

To specify guide rods, simply use ```GuideRods={'radius': 5, 'position': 30, 'height': 370}```. This will add four holes to the corner of your slice. Each hole will have a radius of` radius` and be  `position` mm from the edge of the slice.

## Complete phantom design

Below is the complete script for the phantom we [built]() which incorporates all these elements

```python
from pathlib import Path
import importlib
import numpy as np
import sys
path_to_source = Path(__file__).parent.parent
sys.path.insert(0, str(path_to_source))
from mri_distortion_toolkit import PhantomDesigner
importlib.invalidate_caches()  # maybe is not needed
importlib.reload(PhantomDesigner)

'''
This is the script which I used to generate a design which was sent to Evolution Gear
'''
Nslices = 11 # make this an odd number to make sure you have a slize at z=0
SliceZPositions = np.linspace(-150, 150, Nslices)
SliceThickness = np.mean(np.abs(np.diff(SliceZPositions)))

try:
    import FreeCAD
    run_in_free_cad = True
except ImportError:
    run_in_free_cad = False

for i, z_pos in enumerate(SliceZPositions):
    # setup load:
    if not int(z_pos) == 0 and ((z_pos) < 120) and ((z_pos) > 20):
        load = {'shape': 'rectangle', 'width': 110, 'height': 510}
    else:
        load = None
    # set up crosshair
    if int(z_pos) == 0:
        referenceRadius = 70
    else:
        referenceRadius = None
    # set up end slices:
    if abs(int(z_pos)) == 150:
        HoleStart = 0
    else:
        HoleStart = None

    Slice = PhantomDesigner.PhantomSlice(slice_shape='rectangle',
                                           slice_thickness=SliceThickness, HVL_x=390 / 2, HVL_Y=390 / 2,
                                           hole_depth=17, hole_spacing=25,
                                           hole_radius=8.5 / 2,
                                           DSV=150, z_pos=z_pos,
                                           LoadRegion=load,
                                           GuideRods={'radius': 5, 'position': 20, 'height': 370},
                                           HoleCentroids='ROI_polar',
                                           ReferenceCrosshairRadius=referenceRadius,
                                           bottom_cut=30,
                                           hole_start=HoleStart)
    if run_in_free_cad:
        Slice.draw_slice()
        Slice.add_full_scale_drawing()

    z_array = np.ones(np.shape(Slice.HoleCentroids)[1]) * z_pos
    marker_positions_temp = np.vstack([np.array(Slice.HoleCentroids), z_array])
    try:
        marker_positions = np.hstack([marker_positions, marker_positions_temp])
    except NameError:
        marker_positions = marker_positions_temp

if run_in_free_cad:
    Slice.draw_DSV()
    Slice.draw_Guide()
else:
    marker_positions = np.array(marker_positions)
    np.savetxt(r'evolution_phantom_marker_positions.txt', marker_positions)

```
