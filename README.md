# lepto-blca-colombia

Analysis code for the study:

> **Evaluation of three rapid diagnostic tests and IgM ELISA for human leptospirosis in Colombia: a Bayesian Latent Class Analysis**
> Parra Barrera EL, Marshal G, Bello S, Salas D, Duarte C, Moreno J, Galloway R, Walke H, Undurraga E, Schafer I.
> *PLOS Neglected Tropical Diseases* (under review, 2026).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> **TODO before release:** replace `XXXXXXX` (Zenodo DOI), `<user>` (GitHub username) and `XXXXX` (OSF project id) throughout this file, in `CITATION.cff` and in `.zenodo.json`.

---

## Overview

This repository contains the complete analysis code used to evaluate the diagnostic
accuracy of three commercial rapid diagnostic tests (LeptoCheck-WB, LifeAssay
Test-It™ IgM, SD Bioline *Leptospira* IgG/IgM) and the Panbio® IgM ELISA for human
leptospirosis, using 479 paired acute- and convalescent-phase serum samples collected
through Colombia's national surveillance program (2015–2018).

The analysis combines:

- **Classical (frequentist) accuracy** against the combined MAT result as the
  operational reference standard (sensitivity, specificity, PPV, NPV with exact 95% CIs).
- **Bayesian Latent Class Analysis (BLCA)** estimating true sensitivity, specificity and
  prevalence without assuming any single test (including MAT) is a perfect reference,
  under two model structures (conditional independence and conditional dependence),
  compared by DIC.
- **Secondary analysis** using phase-specific MAT references (MAT-acute, MAT-convalescent).
- **Prior sensitivity analysis** (skeptical / neutral / optimistic scenarios).
- **ROC / AUC analysis** of the Panbio® IgM ELISA (Youden-optimal cut-off; DeLong
  comparison between phases; AUC by days from symptom onset).
- A **Stan replication** of the BLCA models to assess robustness to the MCMC software.

Reporting follows the STARD-BLCM guideline.

---

## Repository structure

```
lepto-blca-colombia/
├── R/                         # Analysis scripts (run in numbered order)
│   ├── 00_setup.R             # Packages, paths, global seed
│   ├── 01_data_preparation.R  # Read data, derive MAT definitions & test coding
│   ├── 02_classical_analysis.R# 2×2 tables vs MAT; Se/Sp/PPV/NPV (Tables 3, 5)
│   ├── 03_blca_jags_independence.R  # Primary BLCA model (Table 3)
│   ├── 04_blca_jags_dependence.R    # Dependence model + DIC comparison (Table 3)
│   ├── 05_secondary_phase_specific_mat.R  # MAT-acute / MAT-convalescent (Table 4)
│   ├── 06_prior_sensitivity.R # Alternative prior scenarios (Table S2)
│   ├── 07_blca_stan.R         # Stan replication (Table S6)
│   ├── 08_roc_analysis.R      # pROC: AUC, Youden, DeLong, AUC by days (Figure 1)
│   └── 09_figures_tables.R    # Build figures and Word tables
├── models/                    # Model definitions
│   ├── blca_independence.jags
│   ├── blca_dependence.jags
│   └── blca_stan.stan
├── data/                      # NO raw data deposited (see data/README.md)
├── outputs/
│   ├── figures/               # Generated figures
│   └── tables/                # Generated tables
├── docs/
│   └── data_dictionary.md     # Variable definitions (for reviewers)
├── CITATION.cff
├── LICENSE
├── .zenodo.json
└── .gitignore
```

---

## Requirements

| Software | Version used |
|----------|--------------|
| R | 4.5.1 |
| JAGS | 4.3.1 |
| R2jags | 0.5-7 |
| rstan | 2.32.7 (Stan 2.32.2) |
| pROC | (as installed) |
| tidyverse, readxl | (as installed) |

JAGS must be installed system-wide before `R2jags` (https://mcmc-jags.sourceforge.io).
A full list of loaded packages and versions is written to `outputs/sessionInfo.txt`
when `00_setup.R` is run.

---

## How to reproduce

1. Install R 4.5.1 and JAGS 4.3.1 (and, for the Stan replication, a working C++
   toolchain for `rstan`).
2. Place the source data file where `01_data_preparation.R` expects it
   (see `data/README.md`). The data are **not** distributed with this repository.
3. Run the scripts in `R/` in numerical order, starting with `00_setup.R`.
   Each script is self-contained and writes its outputs to `outputs/`.

```r
# from the repository root, in R:
source("R/00_setup.R")
source("R/01_data_preparation.R")
# ... through 09_figures_tables.R
```

A fixed random seed is set in `00_setup.R` for reproducibility of the MCMC runs.

---

## Data availability

Consistent with the manuscript's data availability statement, **individual-level data
are not deposited in this repository.** The samples derive from Colombia's national
laboratory-based leptospirosis surveillance system and are held by the Microbiology
Group at the Instituto Nacional de Salud (INS), Bogotá. Aggregate results are reported
in the article and its supplementary material. Requests for access to the underlying
de-identified data should be directed to the corresponding author and the INS, subject
to applicable Colombian regulations (Resolution 8430 of 1993) and the approving ethics
committees.

`docs/data_dictionary.md` documents the variables and coding the scripts expect, so the
pipeline can be inspected and re-run by anyone with authorized access to the data.

---

## Ethics

The study was approved by the U.S. Centers for Disease Control and Prevention
(NCEZID Tracking Number 051817IS) and the Ethics Committee of the Pontificia
Universidad Católica de Chile (Protocol ID 251104016). The INS authorized the use of
archived, de-identified surveillance samples; individual informed consent was waived as
risk-free research under Resolution 8430 of 1993. The study followed the Declaration of
Helsinki.

---

## How to cite

Please cite both the article and this software (see `CITATION.cff`). Once the Zenodo
archive is created, cite the version DOI shown on the Zenodo record.

---

## License

Code is released under the MIT License (see `LICENSE`).

---

## Contact

Eliana L. Parra Barrera — Pontificia Universidad Católica de Chile / Instituto Nacional
de Salud, Colombia. *(add corresponding-author email before release.)*
