# Developer notes

## check docstring coverage with interrogate

from MR_DistortionQA run

```bash
 interrogate -vv --ignore-init-method --generate-badge  ../docsrc/__resources/interrogate.svg
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
