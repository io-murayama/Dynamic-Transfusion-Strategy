# Dynamic-Transfusion-Strategy
This repository contains the full set of R scripts and environment configuration files used to conduct the statistical analyses for a study evaluating dynamic red blood cell transfusion strategies in ICU patients. Using data from eight ICUs in Japan (2013–2025), the scripts estimate the per-protocol effects of clinically relevant dynamic transfusion strategies, defined by hemoglobin thresholds of 9, 8, or 7 g/dL, on seven-day and 28-day all-cause mortality.

---
## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

---

## Overview

This project evaluates the effects of clinically relevant dynamic transfusion strategies on short-term mortality in ICU patients using multicenter observational data from Japan.

All statistical analyses were conducted in R. This repository provides:

- All R scripts used for data processing, modeling, and visualization
- Files required to reproduce the R computational environment (`renv`, Docker)

---

## Repository Structure
```
Dynamic-Transfusion-Strategy
├── README.md
├── LICENSE
├── renv.lock
├── repository.Rproj
├── renv/
│   ├── activate.R
│   └── settings.json
├── rocker/
│   ├── Dockerfile
│   ├── docker-compose.yml
│   └── up.sh
└── scripts/
    ├── 1_concat_csv_and_save_as_rdata.R
    ├── 2_parametric_gformula_simulations.R
    ├── 3_tableone.R
    ├── 4_hb_tf_rate_figures.R
    ├── 5_survival_curves.R
    ├── 6_forest_plots.R
    ├── 7_forest_plots_combined.R
    └── main.sh
```
- renv
  - Contains files required to reproduce the exact R package environment used in the analysis, including project activation and dependency settings.
- rocker
  - Contains a Dockerfile for building a containerized R environment consistent with the analysis setup.
- scripts
  - Contains all R scripts used for data preparation, parametric g-formula simulations, and visualization.

---
## Requirements

### R Environment

- R version: 4.5.1 (R Foundation for Statistical Computing); the exact package versions are recorded in `renv.lock`
- R packages: fully specified in `renv.lock`

---
## Usage

This section describes how to run the parametric g-formula simulations used in the analysis.

### Run parametric g-formula simulations

The main analysis is performed using the script:

```
scripts/2_parametric_gformula_simulations.R
```

This script implements a parametric g-formula framework to estimate per-protocol effects of dynamic transfusion strategies on seven-day and 28-day all-cause mortality. All results are written to the `./output/` directory.

Before execution, analysis settings can be modified in the Configurations section of the script. In particular, the `total_ram` parameter should be set according to the available memory of the machine used for computation.

#### Basic usage
From the project root directory, run:
```
Rscript scripts/2_parametric_gformula_simulations.R --sg {subgroup_name}
```

#### Command-line arguments
- `--sg`
  - Specifies the subgroup to be analyzed. The value must correspond to an entry defined in `subgroup_filters` in the Configurations section of the script.
  - Examples:
    ```
    Rscript scripts/2_parametric_gformula_simulations.R --sg all
    Rscript scripts/2_parametric_gformula_simulations.R --sg age_70_or_older
    ```

- `--no-sim-main`
  - Skips the main simulation with bootstrap resampling. By default, the main simulation (with bootstrap) is executed.

- `--sim-sub`
  - Runs an additional sub-simulation without bootstrap resampling. This option is typically used to generate secondary outputs such as hemoglobin trajectories for figure construction. By default, sub-simulations are not executed.

- `--cov-inv`
  - Performs a sensitivity analysis by reversing the order of time-varying covariates modeled in the g-formula.

- `--time-window-width`
  - Specifies the width of the time window (in hours) used in the longitudinal data structure. The default value is 6. This argument is used for a sensitivity analysis.

---
## Contact
For questions or collaboration inquiries, please reach out to us by email:
 - [MeDiCU, Inc.](mailto:info@medicu.co.jp)

---
## License
This project is licensed under the GNU General Public License (GPL) - see the [LICENSE](https://github.com/io-murayama/Dynamic-Transfusion-Strategy/blob/main/LICENSE) file for details.

---
**Disclaimer:**
The code in this repository is provided for academic research and educational purposes. Individual patient data are not provided.

