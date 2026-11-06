# vitek2-n442-ast-validation

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24015045-blue)](https://doi.org/10.5281/zenodo.24015045)

R analysis pipeline for the manuscript:

> Bouland G, Isnard C, Daurel C, Fines-Guyon M, Le Hello S, Malandain D, Gravey F. **Performance of the VITEK AST-N442 Card for Susceptibility Testing of Hospital-Acquired Enterobacterales: High Rates of Very Major Errors for Ertapenem and Nitrofurantoin.** *Microbiology Spectrum* (Spectrum01157-26).

This repository contains the complete R code used to:

- clean and harmonize VITEK 2 and reference MIC values (out-of-range imputation per ISO 20776-2:2021);
- compute Essential Agreement (EA), bias, Categorical Agreement (CA) and categorical error rates (VME / ME / miE) per ISO 20776-2:2021, ISO 20776-1:2019 and CLSI M52;
- run McNemar tests comparing VITEK 2 vs. reference interpretations (Table 1 of the manuscript);
- compute subgroup error distributions (Enterobacterales groups and resistance phenotypes, Supplementary Tables S6 and S7) using two-sided Fisher's exact tests on 2×2 contingency tables, with odds ratios and exact 95% confidence intervals (**unadjusted p-values; exploratory analyses**);
- compare VITEK 2 Advanced Expert System (AES) calls to WGS-derived resistomes (Supplementary Table S8);
- generate the MIC scattergrams and stacked bar figures (Figure 1).

## Repository structure

```text
.
├── R/
│   ├── nouvelle_norme.R   # EA/bias/CA/VME-ME-miE pipeline, Fisher subgroup analyses (Tables S6/S7), McNemar, scattergrams
│   └── article_vitek.R    # Per-antibiotic error classification, figures, supplementary exports
├── data/                  # Input: put "Tableaux_final_analyses_pour_R.xlsx" here (Supplementary Table S1 of the manuscript)
├── tables/                # Generated (created at runtime)
├── figures/               # Generated (created at runtime)
└── results/               # Generated (created at runtime)
```

## Requirements

- R ≥ 4.3.0
- Packages: `dplyr`, `ggplot2`, `openxlsx`, `readxl`, `stringr`, `tidyr`, `readr`, `rstatix`, `tibble`, `lubridate`, `RColorBrewer`, `xlsx`, `viridis`, `tidyverse`, `gtsummary`, `writexl`, `scales`

## Usage

1. Place the input file `Tableaux_final_analyses_pour_R.xlsx` (i.e., Supplementary Table S1: raw VITEK 2 and reference MICs, interpretations, species, phenotypes, resistomes) in `./data/`.
2. Open `R/nouvelle_norme.R` and `R/article_vitek.R` and run the sections in order. All input/output paths are **relative to the project root** (`./data/`, `./tables/`, `./figures/`, `./results/`), which are created automatically if missing.
3. Raw genomic reads are available in the European Nucleotide Archive under BioProject **PRJEB105945** (study accession ERP187096).

## Statistical notes

- Fisher's exact tests are **two-sided, with unadjusted p-values** (significance threshold p < 0.05). Given the large number of subgroup comparisons, these analyses are exploratory and hypothesis-generating, as stated in the manuscript.
- Odds ratios and their 95% confidence intervals are the conditional maximum-likelihood estimates from `fisher.test()` (exact CIs). ORs are reported as non-estimable (NE) when a contingency cell contains zero events.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).

## Citation

Gravey F, Bouland G. *vitek2-n442-ast-validation: R pipeline for evaluation of the VITEK 2 AST-N442 card against reference MIC methods*. Version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.24015045

Please cite both the manuscript and this release DOI.
