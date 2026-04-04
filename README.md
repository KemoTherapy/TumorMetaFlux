# Metabolic Flux Analysis of Tumors

## Project Overview

This repository contains a computational pipeline for analyzing metabolic flux and metabolite correlations across multiple cancer types using the Cancer Atlas of Metabolic Profiles (cAMP) dataset. The pipeline integrates transcriptomics and metabolomics data through constraint-based metabolic modeling to predict reaction-level fluxes and identify flux-metabolite relationships.

## Dataset

**Source**: Cancer Atlas of Metabolic Profiles (cAMP) dataset  
**Citation**: Benedetti E, et al. A Multimodal Atlas of Tumor Metabolism Reveals the Architecture of Gene-Metabolite Co-regulation. bioRxiv 2022. doi:10.1101/2022.11.23.517549

## Pipeline Structure

The analysis pipeline consists of 7 sequential R scripts:

```
1_InstallMETAFlux.R          → Install METAFlux and dependencies
2_RawData.R                  → Download and organize raw data files
3_DataSetup.R                → Load and structure data for analysis
4_PreprocessingData.R        → Quality control and data filtering
5_FluxCalculations.R         → Run METAFlux to predict metabolic fluxes
6_Correlations.R             → Calculate flux-metabolite correlations
7_Subsystem.R                → Pathway-level aggregation and analysis
```

### Execution Order

Scripts are intended to run sequentially as each depends on outputs from previous steps:

```bash
Rscript 1_InstallMETAFlux.R
Rscript 2_RawData.R
Rscript 3_DataSetup.R
Rscript 4_PreprocessingData.R
Rscript 5_FluxCalculations.R
Rscript 6_Correlations.R
Rscript 7_Subsystem.R
```

## Computational Requirements

- **R version**: ≥ 4.0
- **Memory**: ≥ 16 GB RAM (32 GB recommended for large datasets)
- **Storage**: ~50 GB for raw data and outputs
- **Runtime**: ~6-12 hours for complete pipeline (varies by system)

## Key Dependencies

```r
# Core packages
library(METAFlux)      # Flux balance analysis
library(tidyverse)     # Data manipulation
library(readxl)        # Excel file reading

# Bioconductor packages
BiocManager::install(c("sybil", "sybilSBML"))

# See 1_InstallMETAFlux.R for complete dependency list
```

## Methods Summary

### Flux Calculation (METAFlux)

Metabolic flux is predicted using constraint-based modeling with the Human-GEM genome-scale metabolic model:

**Mathematical formulation**:
- **Objective**: Minimize Σ(vᵢ - aᵢ)² where vᵢ is flux and aᵢ is gene expression-derived activity
- **Constraint**: S·v = 0 (steady-state mass balance)
- **GPR logic**: AND (min) for enzyme complexes, OR (max) for isozymes
- **Medium**: Human blood profile constraints (64 metabolites)

### Statistical Analysis

**Flux-metabolite correlations**:
- **Method**: Pearson correlation across matched samples
- **Testing**: t-test with H₀: ρ = 0
- **Multiple testing**: Benjamini-Hochberg FDR correction (α = 0.05)
- **Quality control**: Minimum 10 samples, remove zero-variance features

**Pathway analysis**:
- **Aggregation**: Mean flux across reactions in each subsystem
- **Tumor vs. normal**: Welch's t-test with log₂ fold change
- **Heterogeneity**: Coefficient of variation (CV) across samples

## Outputs

### Flux Results
- `flux_results_<dataset>.csv`: Predicted fluxes for all reactions × samples
- Dimensions: 13,082 reactions × N samples
- Units: Cubic-root normalized for variance stabilization

### Correlation Results
- `correlation_results_<cancer_type>.xlsx`: Flux-metabolite correlations
  - Sheet 1 (All): All computed correlations
  - Sheet 2 (Significant): FDR < 0.05 only
- Columns: Reaction, Metabolite, n, r, p, FDR

### Pathway Analysis
- `subsystem_activity_<cancer_type>.csv`: Pathway-level flux summaries
- `subsystem_tumor_vs_normal_<cancer_type>.csv`: Differential activity analysis
- `subsystem_heatmap_<cancer_type>.pdf`: Visualization of pathway activities

## Quick Start

### 1. Clone Repository
```bash
git clone <whatever this url is>
cd TumorMetaFlux
```

### 2. Install Dependencies
```r
source("1_InstallMETAFlux.R")
```

### 3. Download Data
```r
source("2_RawData.R")
# Follow prompts to download cAMP data from Zenodo
```

### 4. Run Pipeline
```r
source("3_DataSetup.R")
source("4_PreprocessingData.R")
source("5_FluxCalculations.R")
source("6_Correlations.R")
source("7_Subsystem.R")
```

## Key Results

### Flux Predictions
- Successfully predicted fluxes for 13,082 reactions across 988 tumor samples
- Flux values show expected patterns: high glycolysis, variable OXPHOS across cancer types
- Validates Warburg effect in most cancer types

### Flux-Metabolite Correlations
- Identified thousands of significant flux-metabolite relationships (FDR < 0.05)
- Strong correlations between glycolytic reactions and lactate/pyruvate
- Cancer-type-specific patterns in amino acid and lipid metabolism

### Pathway Analysis
- Core energy metabolism (glycolysis, TCA cycle) shows high activity across all cancers
- Specialized pathways (nucleotide synthesis, lipid metabolism) vary by cancer type
- High inter-patient heterogeneity in amino acid metabolism pathways

### Common Issues

**METAFlux installation fails**:
- Ensure R ≥ 4.0 and Bioconductor is installed
- Try installing sybil and sybilSBML manually first

**Out of memory errors**:
- Increase R memory limit: `memory.limit(size = 32000)` (Windows)
- Process datasets sequentially rather than in parallel

**Sample matching errors**:
- Verify master mapping file exists in `data/`
- Check that sample IDs match across transcriptomics and metabolomics files

**Zero correlations found**:
- Check that transcriptomic and metabolomic data are from same patients
- Verify sample ID mapping using CommonID, MetabID, RNAID columns

## Citation
If you use this pipeline, please cite:

**cAMP dataset**:  
Benedetti E, et al. (2022). A Multimodal Atlas of Tumor Metabolism Reveals the Architecture of Gene-Metabolite Co-regulation. bioRxiv. doi:10.1101/2022.11.23.517549

**Human-GEM**:  
Robinson JL, et al. (2020). An atlas of human metabolism. Science Signaling, 13(624), eaaz1482.

**METAFlux**:  
Huang Y, et al. (2023). Characterizing cancer metabolism from bulk and single-cell RNA-seq data using METAFlux. Nature Communications, 14, 4883.

## Contributing
This is currently a private repository. For questions or collaboration inquiries, please contact the repository owner.
## License
Nothing yet....

## Contact

- **Directory Author**: Kamaleldin Kamaleldin
- **Institution**: University of South Florida & Harvard Medical School / Massachusetts General Hospital
- **Email**: KKamaleldin@usf.edu / KKamaleldin@mgh.harvard.edu

## Acknowledgments
- cAMP consortium for data generation and sharing
- Human-GEM team for metabolic model curation
- METAFlux developers for computational tools
- HMS/MGH for computational resources
