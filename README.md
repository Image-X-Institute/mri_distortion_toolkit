# MRI_DistortionQA  
![](docsrc/__resources/coverage.svg) ![](docsrc/__resources/interrogate.svg)

**Authors:** Brendan Whelan, Paul Liu, Shanshan Shan

This code enables characterization and reporting of geometric distortion in Magnetic Resonance Imaging. The basic end-to-end workflow is below, but all steps have well defined input/output so you can use any part of this code independently from the other parts. For a tutorial on each step, click on the diagram below.


```mermaid
%% Example of sequence diagram
graph LR;
    A[Marker <br>Extraction]-->B[Marker <br>Matching]-->C[Field <br> Calculation]-->D[Spherical Harmonic <br>Analysis]-->E[Automated <br>reporting];
    click A "http://www.github.com"
    click B "http://www.github.com"
    click C "http://www.github.com"
    click D "http://www.github.com"
    click E "http://www.github.com"
    
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
