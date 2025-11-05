#!/bin/bash

# Input QC
./submit.sh --snakefile rules/input_qc.smk results/01_pre_qc/MEX123_input_qc.done

# Remove Ambiguous SNPs
./submit.sh --snakefile rules/remove_ambiguous.smk results/02_remove_ambiguous/MEX123_remove_ambiguous.done

# VCF Convert (all chromosomes 1-22)
./submit.sh --snakefile rules/vcf_convert.smk results/03_vcf/MEX123_convertvcf.done

# Harmonize VCF against reference (parallel by chromosome)
./submit.sh --snakefile rules/harmonize_vcf.smk --jobs 22 --cores 44 results/04_harmonize/MEX123_harmonize.done

# Download imputed results from server
./submit.sh --snakefile rules/download_imputed.smk results/05_imputed_data/MEX123_download_imputed.done