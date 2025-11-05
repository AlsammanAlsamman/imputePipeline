#!/bin/bash
# Submit script for VCF preparation for PLINK conversion

# Run the VCF preparation pipeline
./submit.sh --snakefile rules/prepare_vcf_for_plink.smk results/08_plink_ready/vcf_preparation_complete.done