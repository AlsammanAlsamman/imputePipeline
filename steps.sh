#!/bin/bash

# Dataset name from config
dataset="MEX"
./submit.sh --snakefile rules/input_qc.smk results/01_pre_qc/${dataset}_input_qc.done

# Remove Ambiguous SNPs
./submit.sh --snakefile rules/remove_ambiguous.smk results/02_remove_ambiguous/${dataset}_remove_ambiguous.done

# VCF Convert (all chromosomes 1-22)
./submit.sh --snakefile rules/vcf_convert.smk results/03_vcf/${dataset}_convertvcf.done

# Harmonize VCF against reference (parallel by chromosome)
./submit.sh --snakefile rules/harmonize_vcf.smk --jobs 22 --cores 44 results/04_harmonize/${dataset}_harmonize.done

# Download imputed results from server
./submit.sh --snakefile rules/download_imputed.smk results/05_imputed_data/${dataset}_download_imputed.done

# Extract password-protected imputed files (parallel by chromosome)
./submit.sh --snakefile rules/extract_imputed.smk --jobs 22 --cores 22 results/06_extracted_data/${dataset}_extract_imputed.done

# Filter imputed data from Michigan server (parallel by chromosome)
./submit.sh --snakefile rules/filter_imputed_michigan.smk --jobs 22 results/07_filtered/michigan_filtering_complete.done

# Prepare VCF files for PLINK conversion (parallel by chromosome)
./submit.sh --snakefile rules/prepare_vcf_for_plink.smk --jobs 22 results/08_plink_ready/vcf_preparation_complete.done

# Convert VCF to PLINK1 binary format (parallel by chromosome)
./submit.sh --snakefile rules/convert_to_plink_all.smk --jobs 22 results/09_plink_format/plink_conversion_complete.done

# Merge chromosome PLINK binaries into a single dataset
./submit.sh --snakefile rules/merge_plink2.smk results/10_plink_merged/${dataset}_merged.done

# Convert merged PLINK2 dataset into PLINK1 binary format
./submit.sh --snakefile rules/convert_plink2_to_bed.smk results/10_plink_merged/${dataset}_bed.done
