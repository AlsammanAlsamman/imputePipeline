# Genomic Imputation Pipeline

A Snakemake-based pipeline for genomic imputation using Michigan Imputation Server and TOPMed databases. This pipeline handles pre-imputation quality control, harmonization, and post-imputation processing to produce clean PLINK format datasets.

## Features

- **Pre-imputation QC**: Missing data filtering, MAF filtering, HWE testing
- **Reference harmonization**: Strand checking and allele alignment
- **VCF preparation**: Conversion to imputation-ready format
- **Post-imputation processing**: Quality filtering and format conversion
- **Chromosome parallelization**: SLURM-optimized parallel processing
- **Independent rules**: Each step can be run standalone via sbatch

## Project Structure

```
configs/
├── software.yml     # Software modules and parameters
└── analysis.yml     # Analysis settings and filter thresholds

rules/               # Independent Snakemake rules (to be added)
scripts/             # Helper scripts (to be added)
utils/
└── bioconfigme.py   # Configuration file loader
inputs/              # Input data directory
```

## Quick Start

1. **Configure your analysis**: Edit `configs/analysis.yml` with your dataset paths and parameters
2. **Set software modules**: Update `configs/software.yml` with your cluster's module names
3. **Run individual rules**: Use sbatch to run specific pipeline steps independently

## Configuration

### Analysis Configuration (`configs/analysis.yml`)
- Dataset information and paths
- Quality control thresholds
- Imputation server settings
- Resource requirements

### Software Configuration (`configs/software.yml`)
- Module names for PLINK, bcftools, etc.
- Software parameters and versions
- SLURM resource defaults

## Pipeline Stages

1. **Pre-QC**: Basic quality control filtering
2. **Harmonize**: Reference alignment and strand checking  
3. **VCF Convert**: Preparation for imputation servers
4. **Post-Imputation**: Quality filtering of imputed data
5. **Merge Final**: Chromosome merging and final QC

## Requirements

- SLURM cluster environment
- PLINK (v1.9 or v2.0)
- bcftools, tabix, bgzip
- Python 3.9+ with PyYAML
- Reference genome files

## Output

Final output: Clean PLINK format files (.bed/.bim/.fam) ready for downstream analysis.

---

**DO NOT UPDATE THIS README UNLESS MENTIONED**