"""
This package is for the characterisation and reporting of geometric distortion in MRI
Please see here for detailed documentation
https://acrf-image-x-institute.github.io/MRI_DistortionQA/index.html
"""

__version__ = '0.14.6'
try:
    import FreeCAD
except ImportError:
    # this line crashes free CAD
    from .calculate_harmonics import calculate_harmonics
