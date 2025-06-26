# ipwcch

This is version 0.1.0.

Provides helper functions for evaluating prognostic biomarkers in prospective case-cohort study designs. 
Designed for survival outcomes and dynamic prediction in epidemiologic or clinical research settings.
Supports time-dependent ROC curve analysis and AUC with CI estimation using inverse probability weighting.
This package includes also plotting helper functions to vizualize estimates of absolute and relative risks over time assessed with other packages.

## Requirements

R version minimally 4.4.1 with tidyverse installed.

Following packages need to be installed in order to reproduce the /inst/simulations/benchmark.Rmd:

* pROC
* timeROC
* survival
* mvna
* bshazard
* etm
* mstate
* ggpubr
* kableExtra
* patchwork

## Installation

```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("eutops/ipwcch")
```

## Citation

Please cite:

Charlotte D. Vavourakis, Chiara Herzog, Elisa Redl, Christian Munk, Kirsten Frederiksen, Nora Pashayan, Susanne K. Kjaer, Martin Widschwendter (under review).
The WID-BC: a cervicovaginal DNA methylation-based biomarker for long-term breast cancer risk prediction in young women.
