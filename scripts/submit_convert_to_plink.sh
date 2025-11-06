#!/bin/bash
# Submit script for VCF to PLINK1 conversion

# Run the VCF to PLINK conversion pipeline
./submit.sh --snakefile rules/convert_to_plink.smk results/09_plink_format/plink_conversion_complete.done