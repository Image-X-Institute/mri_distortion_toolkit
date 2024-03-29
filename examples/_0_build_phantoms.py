from pathlib import Path
import importlib
import numpy as np
import sys

path_to_source = Path(__file__).parent.parent
sys.path.insert(0, str(path_to_source))
importlib.invalidate_caches()  # maybe is not needed
from mri_distortion_toolkit import PhantomDesigner
importlib.reload(PhantomDesigner)
try:
    import FreeCAD
    running_in_free_cad = True
except ImportError:
    running_in_free_cad = False

"""
This script was set up to demonstrate some of the different phantom designs that can be generated with this code.
It is intended that you execute it as a macro from inside FreeCAD
"""
different_designs =['elipsoid_cartesian', 'rectangle_polar', 'circle_spiral']
design_to_build = different_designs[0]  # change this to change slice design
SliceZPositions = 0

if design_to_build == 'elipsoid_cartesian':
    Slice = PhantomDesigner.PhantomSlice(slice_shape='ellipse',
                                           slice_thickness=35, HVL_x=200, HVL_Y=160,
                                           hole_depth=17, hole_spacing=20,
                                           hole_radius=8.7 / 2,
                                           z_pos=-SliceZPositions,
                                           HoleCentroids='cartesian',
                                           DSV=None,
                                           bottom_cut=50)


elif design_to_build == 'rectangle_polar':
    Slice = PhantomDesigner.PhantomSlice(slice_shape='rectangle',
                                           slice_thickness=20, HVL_x=100, HVL_Y=100,
                                           hole_depth=17, hole_spacing=20,
                                           hole_radius=5,
                                           z_pos=-SliceZPositions,
                                           HoleCentroids='ROI_polar',
                                           DSV=50,
                                           bottom_cut=0)

elif design_to_build == 'circle_spiral':
    # lets generate a spiral of markers just for fun
    theta = np.radians(np.linspace(0, 360*5, 100))
    r = theta ** 2
    r = (r / np.max(r)) * 180 # scale to slice dimension
    del_ind = r < 30  # remove overlapping centroids
    r = np.delete(r, del_ind)
    theta = np.delete(theta, del_ind)
    r = np.insert(r, 0, 0)  # put back in zero point
    theta = np.insert(theta, 0, 0)
    x_2 = r * np.cos(theta)
    y_2 = r * np.sin(theta)
    hole_centroids = [x_2, y_2]

    Slice = PhantomDesigner.PhantomSlice(slice_shape='ellipse',
                                           slice_thickness=35, HVL_x=200, HVL_Y=200,
                                           hole_depth=17, hole_spacing=20,
                                           hole_radius=8.7 / 2,
                                           z_pos=-SliceZPositions,
                                           HoleCentroids=hole_centroids,
                                           DSV=None,
                                           bottom_cut=0)

else:
    raise AttributeError(f'no method defined for {design_to_build}')

if running_in_free_cad:
    Slice.draw_slice()
    Slice.add_full_scale_drawing()
    Slice.draw_Guide()

