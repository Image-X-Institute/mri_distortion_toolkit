# MRI_DistortionQA  
![](docsrc/__resources/interrogate.svg)

**Authors:** Brendan Whelan, Paul Liu, Shanshan Shan

This code enables characterization and reporting of geometric distortion in Magnetic Resonance Imaging. The basic end-to-end workflow is below, but all steps have well defined input/output so you can use any part of this code independently from the other parts. For a tutorial on each step, click on the diagram below.

```mermaid
flowchart LR
    A[Marker <br>Extraction]--->B[Marker <br>Matching]
    B[Marker <br>Matching]--->C[Field <br> Calculation] & E[Automated <br>reporting]
    C[Field <br> Calculation]-->D[Spherical Harmonic <br>Analysis]
    D[Spherical Harmonic <br>Analysis]-->E[Automated <br>reporting];
    click A "https://acrf-image-x-institute.github.io/MRI_DistortionQA/marker_extraction.html"
    click B "https://acrf-image-x-institute.github.io/MRI_DistortionQA/marker_matching.html"
    click C "https://acrf-image-x-institute.github.io/MRI_DistortionQA/field_calculation.html"
    click D "https://acrf-image-x-institute.github.io/MRI_DistortionQA/fit_spherical_harmonics.html"
    click E "https://acrf-image-x-institute.github.io/MRI_DistortionQA/reporting.html"
```

## Setup/Build/Install

```bash
pip install MRI_DistortionQA  # soon!!
```


## Usage

Detailed documentation is AT SPHINX DOCS

## Directory Structure

- *docs* contains html documentation
- *docsrc* markdown/rst source documentation
- *tests* test cases
- *MRI_DistortionQA* source code
