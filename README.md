# mri_distortion_toolkit  
[![codecov](https://codecov.io/gh/ACRF-Image-X-Institute/mri_distortion_toolkit/branch/main/graph/badge.svg?token=3MCT7S6KVK)](https://codecov.io/gh/ACRF-Image-X-Institute/mri_distortion_toolkit) ![](docsrc/__resources/interrogate.svg)  ![tests](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/run_tests.yml/badge.svg) ![docs](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/build_docs.yml/badge.svg)

This code enables characterization, reporting, and correction of geometric distortion in Magnetic Resonance Imaging.

For the measurement of such distortions, see [here](https://github.com/ACRF-Image-X-Institute/MRI_DistortionPhantom). 

The workflow steps are below, but all steps have well defined input/output so you can use any part of this code independently from the other parts. For a tutorial on each step, click on the diagram below. For an example of our automated reporting template see [here](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/_static/MR_QA_report_20_05_2022.html)

```mermaid
flowchart LR
    A[Marker <br>Extraction]--->B[Marker <br>Matching]
    B[Marker <br>Matching]--->C[Field <br> Calculation] & E[Automated <br>reporting]
    C[Field <br> Calculation]-->D[Spherical Harmonic <br>Analysis]
    D[Spherical Harmonic <br>Analysis]-->E[Automated <br>reporting];
    D[Spherical Harmonic <br>Analysis]-->F[Distortion Correction]

    click A "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MarkerVolume"
    click B "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MatchedMarkerVolumes"
    click C "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#MRI_DistortionQA.FieldCalculation.ConvertMatchedMarkersToBz"
    click D "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#MRI_DistortionQA.FieldAnalysis.SphericalHarmonicFit"
    click E "https://acrf-image-x-institute.github.io/mri_distortion_toolkit/code_docs.html#MRI_DistortionQA.Reports.MRI_QA_Reporter"

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
- *MRI_DistortionQA* source code
- *examples* source code for the [worked examples](https://acrf-image-x-institute.github.io/mri_distortion_toolkit/examples.html)
