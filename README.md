# TMC Breast Cancer Risk — Auxiliary Code

Supporting code and analysis scripts developed for the TMC breast cancer risk
prediction project. This repository holds the **auxiliary material** that builds
toward the model — it is a companion to the model deployment application, not the
deployment itself.

> **No patient data is included.** All scripts read from data files that are
> kept outside this repository. See [Data](#data) below.

## Repository layout

```
code/
  analysis/                       Main analysis workflow
    breast_cancer_risk_model.Rmd    Current end-to-end risk-model analysis (PRIMARY file)
    R/
      reference_imputation.R        Hot-deck imputation for the reference population
  docs/
    icare_example.R                 iCARE package reference example
```

## Data

The analysis scripts expect data files (e.g. the TMC case-control extracts,
reference/control datasets, and incidence/mortality rate tables) that are
**not distributed here** for privacy reasons. In the original working tree the
scripts reference a sibling `Data/` directory. To run the code, place the
required files where the scripts expect them (paths are defined near the top of
each script) and run from the project root.

## Running

The code is written in R (R Markdown for the main analyses). Open
`code/analysis/breast_cancer_risk_model.Rmd` in RStudio and knit, or run the
`.R` scripts directly. The analysis sources
`code/analysis/R/reference_imputation.R`, so keep that path intact. Key package:
`iCARE`; install dependencies as prompted.

## Notes

`breast_cancer_risk_model.Rmd` is the current, primary workflow; it fits the
stratified premenopausal/postmenopausal risk model used by the project.
