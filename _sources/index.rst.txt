.. Template documentation master file, created by
   sphinx-quickstart on Tue May 12 15:57:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MRI Distortion QA
=================
`This code <https://github.com/ACRF-Image-X-Institute/mri_distortion_toolkit>`_ enables characterization and reporting of geometric distortion in Magnetic Resonance Imaging. The basic end-to-end workflow is below, but all steps have well defined input/output so you can use any part of this code independently from the other parts. For a tutorial on each step, click on the diagram below.

.. mermaid::

   flowchart LR

   AA[Phantom Design]

   A[Marker <br>Extraction]--->B[Marker <br>Matching]
   B[Marker <br>Matching]--->C[Field <br> Calculation] & E[Automated <br>reporting]
   C[Field <br> Calculation]-->D[Spherical Harmonic <br>Analysis]
   D[Spherical Harmonic <br>Analysis]-->E[Automated <br>reporting];
   D[Spherical Harmonic <br>Analysis]-->F[Distortion Correction]

   click AA "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/phantom_notes.html"
   click A "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/marker_extraction.html"
   click B "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/marker_matching.html"
   click C "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/field_calculation.html"
   click D "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/fit_spherical_harmonics.html"
   click E "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/reporting.html"



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   examples
   phantom_notes
   code_docs

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
