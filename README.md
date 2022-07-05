# MRI_DistortionQA  
![](docsrc/__resources/coverage.svg) ![](docsrc/__resources/interrogate.svg) ![](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/blob/main/.github/workflows/run_tests.yml/coverage.svg) ![tests](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/run_tests.yml/badge.svg) ![docs](https://github.com/ACRF-Image-X-Institute/MRI_DistortionQA/actions/workflows/build_docs.yml/badge.svg)

**Authors:** Brendan Whelan, Paul Liu, Shanshan Shan

This code enables characterization and reporting of geometric distortion in Magnetic Resonance Imaging. For the measurement of such distortions, see [here](https://github.com/ACRF-Image-X-Institute/MRI_DistortionPhantom). The basic end-to-end workflow is below, but all steps have well defined input/output so you can use any part of this code independently from the other parts. For a tutorial on each step, click on the diagram below. For an example of our automated reporting template see [here](https://acrf-image-x-institute.github.io/MRI_DistortionQA/_static/MR_QA_report_20_05_2022.html)

```mermaid
flowchart LR
    A[Marker <br>Extraction]--->B[Marker <br>Matching]
    B[Marker <br>Matching]--->C[Field <br> Calculation] & E[Automated <br>reporting]
    C[Field <br> Calculation]-->D[Spherical Harmonic <br>Analysis]
    D[Spherical Harmonic <br>Analysis]-->E[Automated <br>reporting];
    click A "https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MarkerVolume"
    click B "https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.MarkerAnalysis.MatchedMarkerVolumes"
    click C "https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.FieldCalculation.ConvertMatchedMarkersToBz"
    click D "https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.FieldAnalysis.SphericalHarmonicFit"
    click E "https://acrf-image-x-institute.github.io/MRI_DistortionQA/code_docs.html#MRI_DistortionQA.Reports.MRI_QA_Reporter"
```

## Setup/Build/Install

```bash
pip install MRI_DistortionQA
```


## Usage

Detailed documentation is [here](https://acrf-image-x-institute.github.io/MRI_DistortionQA/).

## Directory Structure

- *docsrc* markdown/rst source documentation
- *tests* test cases
- *MRI_DistortionQA* source code
- *examples* source code for the [worked examples](https://acrf-image-x-institute.github.io/MRI_DistortionQA/examples.html)
