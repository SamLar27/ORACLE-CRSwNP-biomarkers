# Type-2 Biomarkers in Asthma With Nasal Polyps: Discriminative and Prognostic Value

## Authors

Samuel Mailhot-Larouche\* [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--3028--2315-brightgreen)](https://orcid.org/0000-0002-3028-2315),
Philip Wickham\* [![ORCID](https://img.shields.io/badge/ORCID-0009--0009--0736--0387-brightgreen)](https://orcid.org/0009-0009-0736-0387),
Fleur L. Meulmeester [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--8910--9795-brightgreen)](https://orcid.org/0000-0002-8910-9795),
Morgane Gronnier [![ORCID](https://img.shields.io/badge/ORCID-0009--0002--7138--4518-brightgreen)](https://orcid.org/0009-0002-7138-4518),
Michael E. Wechsler [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3505--2946-brightgreen)](https://orcid.org/0000-0003-3505-2946),
Guy Brusselle [![ORCID](https://img.shields.io/badge/ORCID-0000--0001--7021--8505-brightgreen)](https://orcid.org/0000-0001-7021-8505),
Christopher E. Brightling [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--9345--4903-brightgreen)](https://orcid.org/0000-0002-9345-4903),
Mario Castro [![ORCID](https://img.shields.io/badge/ORCID-0000--0001--6328--8994-brightgreen)](https://orcid.org/0000-0001-6328-8994),
Nicola A. Hanania [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--3087--953X-brightgreen)](https://orcid.org/0000-0003-3087-953X),
David J. Jackson [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--2299--868X-brightgreen)](https://orcid.org/0000-0002-2299-868X),
Neil Martin,
Alison Moore,
Peter Howarth,
Megan E. Hardin,
Cindy Burg [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--8983--6669-brightgreen)](https://orcid.org/0000-0002-8983-6669),
Vennila Dharman [![ORCID](https://img.shields.io/badge/ORCID-0009--0007--1993--0522-brightgreen)](https://orcid.org/0009-0007-1993-0522),
Rebecca Gall,
Paul Klekotka [![ORCID](https://img.shields.io/badge/ORCID-0009--0003--8529--6399-brightgreen)](https://orcid.org/0009-0003-8529-6399),
Mark Holliday [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--2344--6091-brightgreen)](https://orcid.org/0000-0003-2344-6091),
Krystelle Godbout,
Richard W. Beasley [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--0337--406X-brightgreen)](https://orcid.org/0000-0003-0337-406X),
Jacob K. Sont [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--5840--0651-brightgreen)](https://orcid.org/0000-0002-5840-0651),
Ewout W. Steyerberg [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7787--0122-brightgreen)](https://orcid.org/0000-0002-7787-0122),
Ian D. Pavord [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--4288--5973-brightgreen)](https://orcid.org/0000-0002-4288-5973),
Simon Couillard [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--4057--6886-brightgreen)](https://orcid.org/0000-0002-4057-6886)

\* These authors contributed equally.

---

## Overview

This repository contains the complete R analysis code for the paper submitted to *Allergy*. The analysis uses individual participant data (IPD) from multiple randomised controlled trials (the ORACLE consortium) to evaluate the discriminative and prognostic value of three type-2 inflammatory biomarkers — fractional exhaled nitric oxide (FeNO), blood eosinophil count (BEC), and total immunoglobulin E (IgE) — in asthma patients with and without comorbid chronic rhinosinusitis with nasal polyps (CRSwNP).

---

## Repository Structure

```
.
├── Analysis_CRSwNP_Part_A.R   # Data preparation & multiple imputation
├── Analysis_CRSwNP_Part_B.R   # Discriminative value (ROC / AUC analysis)
├── Analysis_CRSwNP_Part_C.R   # Prognostic value (IPD two-stage meta-analysis)
└── README.md
```

### Part A — Data preparation
- Loads and harmonises raw trial data across ORACLE trials
- Applies inclusion/exclusion criteria
- Performs multiple imputation (mice) for missing biomarker values
- Saves `DATA_Raw.RData` and `DATA_imputed.RData` for use in Parts B and C

### Part B — Discriminative value
- Meta-analytic ROC curves (FeNO, BEC, total IgE) for detecting CRSwNP
- Pooled AUC estimation using Rubin's rules on the logit scale (DeLong + random-effects)
- Sensitivity / specificity at selected clinical thresholds
- Meta-analytic Youden index
- Comparison of AUCs across biomarkers (DeLong-based)

### Part C — Prognostic value
- Two-stage IPD meta-analysis: negative binomial models with offset for follow-up duration
- Interaction test: biomarker × CRSwNP status on exacerbation rate
- Subgroup forest plots (With CRSwNP / Without CRSwNP)
- Per-IQR rescaled estimates for FeNO (per 29 ppb), BEC (per 0.27 ×10⁹/L), and IgE (per 361 kU/L)
- Restricted cubic spline analysis of IgE × CRSwNP interaction

---

## Requirements

**R version:** 4.4.2 or later

Install all required packages with:

```r
install.packages(c(
  # Core data manipulation
  "tidyverse", "magrittr", "readxl", "writexl", "lubridate", "glue",

  # Statistics & modelling
  "MASS", "mice", "mitools", "rms", "MuMIn", "boot", "AICcmodavg",
  "caret", "gam", "rstatix", "metafor", "lme4", "glmmTMB", "leaps",
  "pROC", "broom",

  # Survival
  "survival", "survivalAnalysis", "ggsurvfit", "tidycmprsk",

  # Visualization
  "ggpubr", "ggvenn", "cowplot", "pheatmap", "patchwork",
  "gridExtra", "grid", "DiagrammeR", "rsvg", "forplo",

  # Tables & reporting
  "gtsummary", "tableone", "table1", "flextable", "knitr",

  # Spirometry reference equations
  "rspiro"
))
```

> **Note:** The `Gmisc` package is also used; install via `install.packages("Gmisc")`.

---

## How to Run

The scripts must be run **in order** (A → B → C), as each part depends on outputs saved by the previous one.

1. **Set your username** — each script reads `user_name <- Sys.info()[["user"]]` to build file paths automatically. No manual path editing is needed as long as the folder structure mirrors the original.

2. **Run Part A** (`Analysis_CRSwNP_Part_A.R`) — produces `DATA_Raw.RData` and `DATA_imputed.RData`.

3. **Run Part B** (`Analysis_CRSwNP_Part_B.R`) — reads those `.RData` files and saves figures/tables to `Figures_and_Tables/Part_B/`.

4. **Run Part C** (`Analysis_CRSwNP_Part_C.R`) — reads those `.RData` files and saves figures/tables to `Figures_and_Tables/Part_C/`.

> Raw data files are not included in this repository. Output directories (`Figures_and_Tables/`) are also excluded via `.gitignore`.

---

## Citation

> Mailhot-Larouche S, Wickham P, Meulmeester FL, et al. Type-2 Biomarkers in Asthma With Nasal Polyps: Discriminative and Prognostic Value. *Allergy* (submitted).

---

## License

This code is released under the [MIT License](https://opensource.org/licenses/MIT).
