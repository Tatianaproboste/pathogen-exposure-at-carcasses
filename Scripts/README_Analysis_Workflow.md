# Scavenger Network Analysis - Workflow Guide

## Overview
This project analyzes disease transmission risk in scavenger communities using camera trap data and network analysis.

---

## Analysis Pipeline

### **Script Flow:**
```
00_Normalize_species_names.r
         ↓
    (cleaned species data)
         ↓
01_Network_construction.R
         ↓
    (all_networks.rds)
         ↓
02_Network_measures.R
         ↓
    (network metrics)
         ↓
03_Random_network.r (VALIDATION - test if structure is real)
         ↓
    (validated networks)
         ↓
    ┌────┴────┬──────────┐
    ↓         ↓          ↓
   04        05      Extra_KNP.r
Carcass_  Species_
Risk      Centrality
```

**Dependencies:**
- **Script 00:** Standalone - standardizes species names
- **Script 01:** Uses cleaned data from Script 00 → creates `all_networks.rds`
- **Script 02:** Uses `all_networks.rds` → calculates descriptive metrics
- **Script 03:** Uses `all_networks.rds` → validates network structure (null models)
- **Script 04:** Uses validated networks from Script 03 → carcass-level risk analysis
- **Script 05:** Uses validated networks from Script 03 → species-level centrality

**Execution Order:**
1. **Must run sequentially:** 00 → 01 → 02 → 03
2. **Can run in parallel:** 04 & 05 (after Script 03 validation)
3. **Optional:** Extra_KNP.r (location-specific analysis)

---

## Scripts & Outputs

### **00_Normalize_species_names.r**
**What it does:**
- Standardizes species names across datasets
- Resolves taxonomic inconsistencies
- Creates clean species lookup table

**Key Outputs:**
- Cleaned species names in main dataset
- Species taxonomy reference file

**Main Results:**
- Consistent species nomenclature for network analysis

---

### **01_Network_construction.R**
**What it does:**
- Constructs co-occurrence networks from camera trap data
- Creates temporal windows (24h) for species interactions
- Builds networks for each context (season × location × year)

**Key Outputs:**
- `all_networks.rds` - List of igraph network objects (16 networks)
- Network edge lists by context

**Main Results:**
- 16 scavenger networks across different environmental contexts
- Network sizes: 10-40 species, 20-150 edges per network

---

### **02_Network_measures.R**
**What it does:**
- Calculates structural properties for all networks
- Computes global metrics (density, transitivity, modularity)
- Provides descriptive network statistics

**Key Outputs:**
- `network_structural_metrics.csv` - Global network properties
- `network_comparison_plots.png` - Visual comparisons across contexts

**Main Results:**
- Network density varies by season (higher in winter)
- Transitivity appears elevated (suggesting clustering)
- Modularity low/non-significant (weak species grouping)
- **Note:** Statistical validation in Script 03 shows clustering is due to degree heterogeneity, not true community structure

---

### **03_Random_network.r**
**What it does:**
- Tests if observed networks differ from random expectations
- Compares against two null models:
  - **Erdős-Rényi:** Simple random (preserves only nodes + edges)
  - **Configuration Model:** Degree-preserving (accounts for species activity levels)
- Validates non-random network structure

**Key Outputs:**
- `null_model_comparison.csv` - Statistical test results for all networks
- `null_model_comparison.png` - Observed vs. random transitivity plots

**Main Results:**
- Networks show MORE clustering than simple random (Erdős-Rényi: 6/16 networks significant, p < 0.05)
- Networks do NOT differ from degree-preserving random (Configuration: 0/16 significant)
- **Interpretation:** Network structure is explained by species activity heterogeneity (some species are generalists visiting many carcasses, others are specialists), but specific co-occurrence patterns are stochastic (random)

**Biological Meaning:**
- Scavenger communities assemble through neutral, opportunity-driven processes
- No evidence for preferential associations or active avoidance between species
- Disease transmission pathways are unpredictable beyond species activity levels

---

### **04_Carcass_transmissio_Risk.r**
**What it does:**
- Analyzes individual carcasses as disease transmission nodes
- Calculates pathogen-specific risk scores (anthrax, rabies, avian flu, etc.)
- Builds temporal co-occurrence networks within 24h windows

**Key Outputs:**
- `carcass_transmission_metrics.csv` - Risk scores per carcass
- `species_superspreader_risk.csv` - Network-based species rankings
- `network_season_*.png` - Seasonal network plots
- `species_centrality_heatmap.png` - Species roles by centrality metrics

**Main Results:**
- Which carcasses are high-risk for disease spread
- Pathogen-specific transmission potential varies by scavenger composition
- Seasonal variation in carcass-level network structure

**Note:** Interpretation informed by Script 03 validation (stochastic assembly with degree heterogeneity)

---

### **05_Disease_transmission_Risk.r** (formerly Species_centrality.r)
**What it does:**
- Identifies superspreader species across all contexts
- Compares species importance by functional group (carnivores, corvids, raptors, etc.)
- Analyzes centrality metrics (degree, betweenness, closeness, eigenvector)
- Tracks how species importance changes across seasons/locations

**Key Outputs:**
- `species_centrality_all.csv` - Full dataset with functional groups
- `top_species_by_functional_group.csv` - Top 3 per group
- `heatmap_improved_composite.png` - Important species by location×season
- `ranking_dynamics_improved.png` - How top 5 species change across contexts
- `multi_metric_heatmap_important.png` - Centrality metrics for consistently important species (≥50% contexts)
- `multi_metric_heatmap_all.png` - All species appearing in top 10

**Main Results:**
- ~20 species consistently important for transmission (≥50% of contexts)
- Corvids (ravens, magpies): high degree centrality → contact hubs
- Carnivores (dingo, fox): high betweenness → bridge different groups
- Raptors (eagles): high closeness → rapid disease spread potential
- Species roles vary by season/location (e.g., corvids more important in summer)

**Note:** Centrality metrics reflect species activity levels (validated by Script 03) rather than true preferential associations

---

### **Extra_KNP.r**
**What it does:**
- Supplementary analysis specific to Kruger National Park (KNP) data
- Location-specific deep dive or methodological tests

**Key Outputs:**
- Location-specific results (varies)

**Main Results:**
- Supplementary findings for KNP context

---

## Key Findings Summary

### **Disease Transmission Risk:**

#### **1. Top superspreader species** (from Script 05):
| Functional Group | Species | Role | Centrality |
|-----------------|---------|------|------------|
| Corvids | Raven, Magpie, Currawong | Contact hubs | High degree |
| Carnivores | Dingo, Fox, Cat | Network bridges | High betweenness |
| Raptors | Wedge-Tailed Eagle, Falcon | Fast spreaders | High closeness |
| Herbivores | Red Kangaroo, Deer | Occasional visitors | Low centrality |

#### **2. Carcass-level risk** (from Script 04):
- Mean 5-8 species per carcass
- High-risk carcasses: >10 species + rapid succession (<6h gaps between visitors)
- Pathogen-specific risk varies by functional group composition:
  - **Anthrax risk:** Higher with large herbivores + carnivores
  - **Rabies risk:** Higher with carnivore-dominated assemblages
  - **Avian flu risk:** Higher with corvid + raptor co-occurrence

#### **3. Network structure** (from Script 03):

**Observed patterns (Script 02):**
- Networks show clustering (high transitivity)
- Appear more structured than random

**Statistical validation (Script 03):**
- Clustering is **real** vs. simple random (Erdős-Rényi: 6/16 significant, p < 0.05)
- BUT clustering is **explained by degree distribution** (Configuration: 0/16 significant)

**Biological interpretation:**
- **Apparent structure:** Some species (corvids, carnivores) visit many carcasses → create clustering
- **True randomness:** Which specific species meet at each carcass is random
- **Conclusion:** Stochastic assembly with degree heterogeneity (no preferential associations)

**What this means:**
- Networks LOOK structured (high transitivity)
- But this is just because some species are active generalists
- Specific co-occurrence patterns are random → unpredictable transmission

---

## Interpretation for Disease Ecology

### **What the results mean:**
**Random assembly + degree heterogeneity = unpredictable transmission**

1. **Some species are active generalists** (visit many carcasses) → natural superspreaders
2. **Some species are rare specialists** (visit few carcasses) → low transmission risk
3. **BUT specific pairings are random** → disease can spread anywhere in network once introduced
4. **No "safe" species pairs** or predictable transmission chains

### **Management Implications:**
✅ **Surveillance:** Must monitor entire scavenger community (not just carnivores or single taxa)  
✅ **Outbreak prediction:** Cannot rely on species pair associations; need community-wide tracking  
✅ **Seasonal risk:** Winter networks denser → higher transmission potential  
✅ **Carcass management:** 24h removal window critical for most pathogens  
✅ **Targeted control:** Focus on generalist species (corvids, common carnivores) but expect spillover to any species

---

## File Dependencies

### **Key Data Files:**
| File | Created by | Used by | Description |
|------|-----------|---------|-------------|
| `cleaned_species_data` | Script 00 | Script 01 | Standardized species names |
| `all_networks.rds` | Script 01 | Scripts 02, 03, 04, 05 | Network objects for all contexts |
| `network_structural_metrics.csv` | Script 02 | Script 04 | Global network properties |
| `null_model_comparison.csv` | Script 03 | Scripts 04, 05 | Validation results (inform interpretation) |

### **Outputs Location:**
- `./Outputs/Networks/Structure_Analysis/` (Scripts 01, 02)
- `./Outputs/Networks/Random_Tests/` (Script 03 - validation)
- `./Outputs/Disease_Transmission/` (Script 04 - carcass risk)
- `./Outputs/Networks/Disease_Risk/` (Script 05 - species centrality)

---

## Quick Start

### **To run full analysis:**
```r
# 1. Standardize species names
source("00_Normalize_species_names.r")

# 2. Build networks (constructs co-occurrence networks)
source("01_Network_construction.R")

# 3. Calculate descriptive network metrics
source("02_Network_measures.R")

# 4. Validate network structure (null model tests) - CRITICAL STEP
source("03_Random_network.r")

# 5. Analyze carcass-level transmission risk
source("04_Carcass_transmissio_Risk.r")

# 6. Analyze species-level centrality (superspreaders)
source("05_Disease_transmission_Risk.r")

# Optional: KNP-specific analysis
source("Extra_KNP.r")
```

### **To run only specific analyses:**
```r
# Just species centrality (requires all_networks.rds from Script 01):
source("05_Disease_transmission_Risk.r")

# Just null model tests (requires all_networks.rds from Script 01):
source("03_Random_network.r")

# Just carcass risk (requires network metrics from Script 02):
source("04_Carcass_transmissio_Risk.r")
```

**Note:** Scripts 04 & 05 should be run AFTER Script 03 validation to properly interpret results in light of network structure findings.

---

## Research Questions Answered

| Script | Research Question | Key Finding |
|--------|------------------|-------------|
| **00** | Are species names consistent? | Standardized nomenclature across datasets |
| **01** | How do scavenger networks vary? | 16 networks across contexts; 10-40 species each |
| **02** | What are network-level patterns? | Higher density in winter; apparent clustering |
| **03** | Are networks structured? | Stochastic assembly with degree heterogeneity |
| **04** | Which carcasses are high-risk? | >10 species + <6h succession = high risk |
| **05** | Which species are superspreaders? | Corvids (hubs), carnivores (bridges), raptors (speed) |

---

## Software Requirements

```r
# Required packages
install.packages(c(
  "igraph",         # Network analysis
  "data.table",     # Fast data manipulation
  "tidyverse",      # Data wrangling & visualization
  "viridis",        # Color palettes
  "pheatmap",       # Heatmaps
  "patchwork"       # Multi-panel plots
))
```

**Minimum versions:**
- R ≥ 4.0.0
- igraph ≥ 1.3.0
- tidyverse ≥ 1.3.0

---

## Troubleshooting

### **Common Errors:**

**1. "all_networks.rds not found"**
- **Solution:** Run Script 01 first to create the file
- **Check:** File should be in `./Outputs/Networks/Structure_Analysis/`

**2. "Network too small" warnings**
- **Cause:** Some contexts have <3 species
- **Solution:** These are automatically skipped; expected behavior
- **Impact:** Reduces total networks analyzed (some contexts excluded)

**3. "Cannot calculate closeness centrality"**
- **Cause:** Disconnected network components
- **Solution:** Scripts use normalized metrics; handled automatically
- **Interpretation:** Species in separate components cannot transmit to each other

**4. "Negative betweenness values"**
- **Cause:** Rare edge case in igraph calculation
- **Solution:** Re-run with different random seed or update igraph package
- **Alternative:** Use `betweenness(network, normalized = TRUE)`

**5. "Error in pheatmap()" - package not found**
- **Solution:** Install missing package: `install.packages("pheatmap")`

**6. File name mismatch errors**
- **Cause:** Script 04 has typo in filename (`transmissio` instead of `transmission`)
- **Solution:** Use exact filename `04_Carcass_transmissio_Risk.r` or rename file to fix typo

---

## Workflow Validation Checklist

Before running analyses, verify:
- [ ] Script 00 has cleaned species names
- [ ] Script 01 has created `all_networks.rds` (check file exists)
- [ ] Script 02 has produced `network_structural_metrics.csv`
- [ ] Script 03 validation completed (check `null_model_comparison.csv` exists)
- [ ] All required packages installed (see Software Requirements)

After Script 03 validation:
- [ ] Review `null_model_comparison.png` - are networks significantly different from random?
- [ ] Check console output - how many networks passed validation?
- [ ] Read interpretation notes - understand what centrality metrics mean given validation results

---

## Citation

If using these analyses, please cite:

> [Your Name]. (2025). Disease transmission risk in scavenger networks: A multi-scale network approach combining carcass-level risk assessment and species centrality analysis. *Journal Name*.

**Key Methods References:**
- **Null model approach:** Bascompte, J., et al. (2003). The nested assembly of plant-animal mutualistic networks. *Proceedings of the National Academy of Sciences*, 100(16), 9383-9387.
- **Configuration model:** Maslov, S., & Sneppen, K. (2002). Specificity and stability in topology of protein networks. *Science*, 296(5569), 910-913.
- **Centrality metrics:** Delmas, E., et al. (2019). Analyzing ecological networks of species interactions. *Biological Reviews*, 94(1), 16-36.
- **Temporal networks:** Blonder, B., et al. (2012). Temporal dynamics and network analysis. *Ecology Letters*, 15(8), 865-877.

---

## Contact

For questions about the workflow or methods:
- **Author:** [Your Name]
- **Email:** [your.email@sydney.edu.au]
- **Institution:** The University of Sydney
- **GitHub:** [your-repo-link] (if applicable)

**Report issues:**
- Data processing questions → [email]
- Statistical methods → [email]
- Network interpretation → [email]

---

## Version History

- **v1.0** (Dec 2025): Initial analysis pipeline
  - Scripts 00-02: Data preparation and network construction
  
- **v1.1** (Dec 2025): Added disease transmission analyses
  - Script 04: Carcass-level risk assessment
  - Script 05: Species-level centrality analysis with functional groups
  
- **v1.2** (Dec 2025): Added validation framework
  - Script 03: Null model validation (Erdős-Rényi + Configuration model)
  - Reordered scripts to validate before interpretation
  - Updated README to reflect stochastic assembly findings

- **v1.3** (Dec 2025): Documentation improvements
  - Comprehensive README with troubleshooting
  - Clarified script dependencies and execution order
  - Added validation checklist

---

## Acknowledgments

This work was supported by [Funding Source]. We thank [Collaborators] for data collection and [Advisors] for methodological guidance.

Camera trap data collection protocols followed [Ethics Approval Number].

---

*Last updated: December 9, 2025*