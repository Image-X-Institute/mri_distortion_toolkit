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
                                           hole_depth=19, hole_spacing=25,
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