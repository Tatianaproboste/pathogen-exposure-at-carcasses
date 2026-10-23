# Pathogen exposure at carcasses

Clean reviewer-ready repository for reproducing the analysis and figures from the paper.

## Repository structure

- `Scripts/` — analysis scripts (methods unchanged).
- `Raw_data/` — place required input datasets here locally (not tracked by git).
- `Outputs/` — generated outputs (tables, RDS files, figures; not tracked by git).

## Reproducibility workflow

Run scripts in this order:

1. `Scripts/00_Normalize species names.r`
2. `Scripts/01_Network_construction.R`
3. `Scripts/02_Network_measures.R`
4. `Scripts/03_Disease_Exposure_Risk_Models.R`
5. `Scripts/04_Species_centrality.r`
6. `Scripts/05_Order_Visualization.r`
7. `Scripts/06_Pairwise_precedence.R`

## Data placement

- Keep the original expected file/folder paths used by the scripts under `Raw_data/`.
- This repository intentionally ignores `Raw_data/` and `Outputs/` contents to avoid publishing raw/private data and generated files.

## Notes

- Script logic, functions, and methods were not changed.
- Only repository structure and presentation were cleaned for reviewer use.
