# mri_distortion_toolkit  
[![codecov](https://codecov.io/gh/ACRF-Image-X-Institute/mri_distortion_toolkit/branch/main/graph/badge.svg?token=3MCT7S6KVK)](https://codecov.io/gh/ACRF-Image-X-Institute/mri_distortion_toolkit) ![](docsrc/__resources/interrogate.svg)  ![tests](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/run_tests.yml/badge.svg) ![docs](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/build_docs.yml/badge.svg)[![PyPI version](https://badge.fury.io/py/mri_distortion_toolkit.svg)](https://badge.fury.io/py/mri_distortion_toolkit)

This code enables characterization, reporting, and correction of geometric distortion in Magnetic Resonance Imaging.

The workflow steps are below. All steps have well defined input/output so you can use any part of this code independently from the other parts. For an example of our automated reporting template see [here](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/_static/MR_QA_report_20_05_2022.html)

```mermaid
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
```



## Setup/Build/Install

```bash
pip install mri_distortion_toolkit
```

## Usage

Detailed documentation is [here](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/). 

## Directory Structure

- *docsrc* markdown/rst source documentation
- *tests* test cases
- *mri_distortion_toolkit* source code 
- *examples* source code for the [worked examples](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/examples.html)
