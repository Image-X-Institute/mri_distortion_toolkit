# Developer notes

## check docstring coverage with interrogate

from root directory:

```bash
 interrogate -c setup.cfg -v MRI_DistortionQA/ --generate-badge  docsrc/__resources/interrogate.svg
```

to run tests, from root directory:
```bash
pytest
```

to assess test coverage, from root directory:
```bash
coverage run -m pytest
```

following that, to update the test coverage badge:

```bash
coverage-badge -f -o docsrc/__resources/coverage.svg
```
