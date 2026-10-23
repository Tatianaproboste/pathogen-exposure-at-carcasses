# Scavenger community assembly and host susceptibility shapes pathogen exposure at carcasses

## Repository structure

- `Scripts/` — analysis scripts
- `Raw_data/` — place required input datasets here locally
- `Outputs/` — generated outputs
- 
## Reproducibility workflow

Run scripts in this order:

1. `Scripts/01_Network_construction.R`
2. `Scripts/02_Network_measures.R`
3. `Scripts/03_Disease_Exposure_Risk_Models.R`
4. `Scripts/04_Species_centrality.r`
5. `Scripts/05_Pairwise_precedence.R`

## Data placement

- Keep the original expected file/folder paths used by the scripts under `Raw_data/`.
- This repository intentionally ignores `Raw_data/` and `Outputs/` contents to avoid publishing raw/private data and generated files.

# Artificial Intelligence (AI) statement
GitHub Copilot (agent mode), Claude Sonnet 4.6 (Anthropic), and GPT-5.3 Codex (OpenAI) were used to assist with writing, debugging and polishing and organising R code used for data processing, network construction, disease exposure risk modelling and visualizations. These tools were not used to design the analytical approach, select model parameters, or interpret results. All code was reviewed, tested, and validated by the authors, who take full responsibility for its accuracy and for the results presented in this manuscript.
 
