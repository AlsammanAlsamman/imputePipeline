""""""""""""

Snakemake rule for converting PLINK-ready VCF files to PLINK1 binary format.

"""Snakemake rule for converting PLINK-ready VCF files to PLINK1 binary format.



import osSnakemake rule for converting PLINK-ready VCF files to PLINK1 binary format.Snakemake rule for converting filtered Michigan VCF files to PLINK1 format.

import sys

This rule converts filtered and PLINK-ready VCF files to PLINK1 binary format

# Add utils directory to Python path

utils_path = "utils"(.bed/.bim/.fam) with post-imputation quality control filters applied during conversion.

if utils_path not in sys.path:

    sys.path.insert(0, utils_path)



from bioconfigme import get_config_value, get_software_moduleAuthor: Genomics PipelineThis rule converts filtered and PLINK-ready VCF files to PLINK1 binary formatThis rule handles the conversion of filtered imputed VCF files to PLINK1 format



# Get configuration valuesDate: November 2025

DATASET_NAME = get_config_value("configs/analysis.yml", "dataset.name")

N_CHROMOSOMES = get_config_value("configs/analysis.yml", "dataset.chromosomes")"""(.bed/.bim/.fam) with post-imputation quality control filters applied during conversion.with robust handling of multi-allelic variants and duplicate IDs.

INFO_SCORE = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.info_score")

MAF = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.maf")

MISSING_RATE = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.missing_rate")

PLINK1_MODULE = get_software_module("plink1")import os



# Define directoriesimport sys

PLINK_READY_DIR = "results/08_plink_ready"

PLINK_FORMAT_DIR = "results/09_plink_format"Author: Genomics PipelineAuthor: Genomics Pipeline



rule convert_to_plink_chr:# Add utils directory to Python path

    """Convert individual chromosome VCF to PLINK1 binary format."""

    input:utils_path = "utils"Date: November 2025Date: November 2025

        plink_ready_vcf = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz",

        plink_ready_index = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz.csi"if utils_path not in sys.path:

    output:

        plink_bed = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed",    sys.path.insert(0, utils_path)""""""

        plink_bim = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim", 

        plink_fam = f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam",

        conversion_log = f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log"

    params:from bioconfigme import get_config_value, get_software_module

        chromosome = "{chr}",

        dataset = DATASET_NAME,

        info_score = INFO_SCORE,

        maf = MAF,# Load configuration valuesimport osimport os

        missing_rate = MISSING_RATE,

        plink1_module = PLINK1_MODULE,def get_dataset_info():

        plink_prefix = f"{PLINK_FORMAT_DIR}/chr{{chr}}"

    log:    """Extract dataset information from configuration"""import sysimport sys

        f"results/log/convert_to_plink_chr{{chr}}_{DATASET_NAME}.log"

    threads: 2    dataset_name = get_config_value("configs/analysis.yml", "dataset.name")

    resources:

        mem_mb = 32000,    chromosomes = get_config_value("configs/analysis.yml", "dataset.chromosomes")

        time = "00:30:00"

    shell:    return dataset_name, chromosomes

        """

        module load {params.plink1_module}# Add utils directory to Python path# Add utils directory to Python path

        mkdir -p {PLINK_FORMAT_DIR}

        mkdir -p results/logdef get_post_imputation_filters():

        

        echo "Starting VCF to PLINK1 conversion for chromosome {params.chromosome}" > {log}    """Extract post-imputation QC filter values from configuration"""utils_path = "utils"utils_path = "utils"

        echo "Dataset: {params.dataset}" >> {log}

        echo "Loaded module: {params.plink1_module}" >> {log}    info_score = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.info_score")

        echo "QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}

        echo "Started at: $(date)" >> {log}    maf = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.maf")if utils_path not in sys.path:if utils_path not in sys.path:

        echo "" >> {log}

            missing_rate = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.missing_rate")

        bash scripts/convert_to_plink.sh \\

            {input.plink_ready_vcf} \\    return info_score, maf, missing_rate    sys.path.insert(0, utils_path)    sys.path.insert(0, utils_path)

            {params.plink_prefix} \\

            {output.conversion_log} \\

            {params.chromosome} \\

            {params.info_score} \\# Get configuration values

            {params.maf} \\

            {params.missing_rate} \\DATASET_NAME, N_CHROMOSOMES = get_dataset_info()

            {threads} 2>&1 | tee -a {log}

        INFO_SCORE, MAF, MISSING_RATE = get_post_imputation_filters()from bioconfigme import get_config_value, get_software_modulefrom bioconfigme import get_config_value, get_software_module

        echo "" >> {log}

        echo "PLINK conversion completed at: $(date)" >> {log}

        """

# Get software modules

rule convert_to_plink_all:

    """Convert all chromosomes and create completion marker."""PLINK1_MODULE = get_software_module("plink1")

    input:

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed", chr=range(1, N_CHROMOSOMES + 1)),# Load configuration values# Load configuration values

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim", chr=range(1, N_CHROMOSOMES + 1)),

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam", chr=range(1, N_CHROMOSOMES + 1)),# Define input and output directories

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log", chr=range(1, N_CHROMOSOMES + 1))

    output:PLINK_READY_DIR = "results/08_plink_ready"def get_dataset_info():def get_dataset_info():

        done_marker = f"{PLINK_FORMAT_DIR}/plink_conversion_complete.done"

    params:PLINK_FORMAT_DIR = "results/09_plink_format"

        dataset = DATASET_NAME,

        n_chr = N_CHROMOSOMES,    """Extract dataset information from configuration"""    """Extract dataset information from configuration"""

        info_score = INFO_SCORE,

        maf = MAF,rule convert_to_plink_chr:

        missing_rate = MISSING_RATE

    log:    """    dataset_name = get_config_value("configs/analysis.yml", "dataset.name")    dataset_name = get_config_value("configs/analysis.yml", "dataset.name")

        f"results/log/convert_to_plink_all_{DATASET_NAME}.log"

    resources:    Convert individual chromosome PLINK-ready VCF files to PLINK1 binary format.

        mem_mb = 32000,

        time = "00:30:00"        chromosomes = get_config_value("configs/analysis.yml", "dataset.chromosomes")    chromosomes = get_config_value("configs/analysis.yml", "dataset.chromosomes")

    shell:

        """    This rule:

        echo "Finalizing PLINK conversion for dataset {params.dataset}" > {log}

        echo "Processed {params.n_chr} chromosomes" >> {log}    1. Converts VCF to PLINK1 binary format (.bed/.bim/.fam)    return dataset_name, chromosomes    return dataset_name, chromosomes

        echo "Applied QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}

        echo "Started at: $(date)" >> {log}    2. Applies post-imputation QC filters during conversion

        

        echo "PLINK Conversion Complete" > {output.done_marker}    3. Validates the output PLINK files

        echo "========================" >> {output.done_marker}

        echo "Dataset: {params.dataset}" >> {output.done_marker}    """

        echo "Date: $(date)" >> {output.done_marker}

        echo "Chromosomes processed: {params.n_chr}" >> {output.done_marker}    input:def get_post_imputation_filters():# Get configuration values

        echo "" >> {output.done_marker}

        echo "QC Filters Applied:" >> {output.done_marker}        # PLINK-ready VCF file from previous preparation step

        echo "- INFO score threshold: {params.info_score}" >> {output.done_marker}

        echo "- MAF threshold: {params.maf}" >> {output.done_marker}        plink_ready_vcf = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz",    """Extract post-imputation QC filter values from configuration"""DATASET_NAME, N_CHROMOSOMES = get_dataset_info()

        echo "- Missing rate threshold: {params.missing_rate}" >> {output.done_marker}

        echo "" >> {output.done_marker}        plink_ready_index = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz.csi"

        

        echo "File Verification:" >> {output.done_marker}    output:    info_score = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.info_score")

        COMPLETE_COUNT=0

                # PLINK1 binary format files

        for chr in {{1..{params.n_chr}}}; do

            if [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bed ]] && \\        plink_bed = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed",    maf = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.maf")# Get software modules

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bim ]] && \\

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.fam ]]; then        plink_bim = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim", 

                echo "Chromosome ${{chr}}: Complete (.bed/.bim/.fam)" >> {output.done_marker}

                COMPLETE_COUNT=$((COMPLETE_COUNT + 1))        plink_fam = f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam",    missing_rate = get_config_value("configs/analysis.yml", "qc_filters.post_imputation.missing_rate")PLINK_MODULE = get_software_module("plink")

            else

                echo "Chromosome ${{chr}}: ERROR - Missing PLINK files" >> {output.done_marker}        # Conversion log

            fi

        done        conversion_log = f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log"    return info_score, maf, missing_rateBCFTOOLS_MODULE = get_software_module("bcftools")

        

        echo "" >> {output.done_marker}    params:

        echo "Summary: $COMPLETE_COUNT/{params.n_chr} chromosomes successfully converted" >> {output.done_marker}

        echo "Output directory: {PLINK_FORMAT_DIR}" >> {output.done_marker}        # Chromosome number

        echo "Ready for downstream analysis or merging." >> {output.done_marker}

                chromosome = "{chr}",

        echo "PLINK conversion finalization completed at: $(date)" >> {log}

        echo "Completion marker: {output.done_marker}" >> {log}        # Dataset name for logging# Get configuration values# Define input and output directories

        echo "Successfully converted: $COMPLETE_COUNT/{params.n_chr} chromosomes" >> {log}

        """        dataset = DATASET_NAME,

        # QC filter valuesDATASET_NAME, N_CHROMOSOMES = get_dataset_info()FILTERED_DIR = "results/07_filtered"

        info_score = INFO_SCORE,

        maf = MAF,INFO_SCORE, MAF, MISSING_RATE = get_post_imputation_filters()PLINK_DIR = "results/08_plink_conversion"

        missing_rate = MISSING_RATE,

        # Software module to load

        plink1_module = PLINK1_MODULE,

        # PLINK output prefix# Get software modulesrule convert_vcf_to_plink_chr:

        plink_prefix = f"{PLINK_FORMAT_DIR}/chr{{chr}}"

    log:PLINK1_MODULE = get_software_module("plink1")    """

        f"results/log/convert_to_plink_chr{{chr}}_{DATASET_NAME}.log"

    threads: 2    Convert individual chromosome VCF files to PLINK1 format.

    resources:

        mem_mb = 32000,# Define input and output directories    

        time = "00:30:00"

    shell:PLINK_READY_DIR = "results/08_plink_ready"    This rule:

        """

        # Load required software modulePLINK_FORMAT_DIR = "results/09_plink_format"    1. Removes multi-allelic variants from filtered VCF files

        module load {params.plink1_module}

            2. Converts clean VCF to PLINK1 (.bed/.bim/.fam) format

        # Create output directory

        mkdir -p {PLINK_FORMAT_DIR}rule convert_to_plink_chr:    3. Handles duplicate IDs and other PLINK conversion issues

        mkdir -p results/log

            """    """

        # Start logging

        echo "Starting VCF to PLINK1 conversion for chromosome {params.chromosome}" > {log}    Convert individual chromosome PLINK-ready VCF files to PLINK1 binary format.    input:

        echo "Dataset: {params.dataset}" >> {log}

        echo "Loaded module: {params.plink1_module}" >> {log}            # Filtered VCF file from previous filtering step

        echo "QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}

        echo "Started at: $(date)" >> {log}    This rule:        filtered_vcf = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz",

        echo "" >> {log}

            1. Converts VCF to PLINK1 binary format (.bed/.bim/.fam)        filtered_index = f"{FILTERED_DIR}/chr{{chr}}.filtered.vcf.gz.csi"

        # Use helper script for conversion

        bash scripts/convert_to_plink.sh \\    2. Applies post-imputation QC filters during conversion    output:

            {input.plink_ready_vcf} \\

            {params.plink_prefix} \\    3. Validates the output PLINK files        # PLINK1 format files

            {output.conversion_log} \\

            {params.chromosome} \\    """        plink_bed = f"{PLINK_DIR}/chr{{chr}}.bed",

            {params.info_score} \\

            {params.maf} \\    input:        plink_bim = f"{PLINK_DIR}/chr{{chr}}.bim", 

            {params.missing_rate} \\

            {threads} 2>&1 | tee -a {log}        # PLINK-ready VCF file from previous preparation step        plink_fam = f"{PLINK_DIR}/chr{{chr}}.fam",

        

        echo "" >> {log}        plink_ready_vcf = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz",        # Clean VCF intermediate file

        echo "PLINK conversion completed at: $(date)" >> {log}

        """        plink_ready_index = f"{PLINK_READY_DIR}/chr{{chr}}.plink_ready.vcf.gz.csi"        clean_vcf = f"{PLINK_DIR}/chr{{chr}}.clean.vcf.gz",



rule convert_to_plink_all:    output:        # Conversion summary

    """

    Convert all chromosomes from VCF to PLINK1 binary format.        # PLINK1 binary format files        conversion_summary = f"{PLINK_DIR}/chr{{chr}}.conversion_summary.txt"

    

    This rule coordinates the conversion of all autosomes (chromosomes 1-22)        plink_bed = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed",    params:

    and creates a completion marker when all conversions are done.

    """        plink_bim = f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim",         # Chromosome number

    input:

        # Wait for all individual chromosome conversions to complete        plink_fam = f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam",        chromosome = "{chr}",

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed", chr=range(1, N_CHROMOSOMES + 1)),

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim", chr=range(1, N_CHROMOSOMES + 1)),        # Conversion log        # Dataset name for logging

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam", chr=range(1, N_CHROMOSOMES + 1)),

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log", chr=range(1, N_CHROMOSOMES + 1))        conversion_log = f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log"        dataset = DATASET_NAME,

    output:

        # Completion marker    params:        # Software modules to load

        done_marker = f"{PLINK_FORMAT_DIR}/plink_conversion_complete.done"

    params:        # Chromosome number        plink_module = PLINK_MODULE,

        dataset = DATASET_NAME,

        n_chr = N_CHROMOSOMES,        chromosome = "{chr}",        bcftools_module = BCFTOOLS_MODULE

        info_score = INFO_SCORE,

        maf = MAF,        # Dataset name for logging    log:

        missing_rate = MISSING_RATE

    log:        dataset = DATASET_NAME,        f"results/log/convert_plink_chr{{chr}}_{DATASET_NAME}.log"

        f"results/log/convert_to_plink_all_{DATASET_NAME}.log"

    resources:        # QC filter values    threads: 2

        mem_mb = 32000,

        time = "00:30:00"        info_score = INFO_SCORE,    resources:

    shell:

        """        maf = MAF,        mem_mb = 32000,

        echo "Finalizing PLINK conversion for dataset {params.dataset}" > {log}

        echo "Processed {params.n_chr} chromosomes" >> {log}        missing_rate = MISSING_RATE,        time = "00:30:00"

        echo "Applied QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}

        echo "Started at: $(date)" >> {log}        # Software module to load    shell:

        

        # Create completion marker with summary information        plink1_module = PLINK1_MODULE,        """

        echo "PLINK Conversion Complete" > {output.done_marker}

        echo "========================" >> {output.done_marker}        # PLINK output prefix        # Load required software modules

        echo "Dataset: {params.dataset}" >> {output.done_marker}

        echo "Date: $(date)" >> {output.done_marker}        plink_prefix = f"{PLINK_FORMAT_DIR}/chr{{chr}}"        module load {params.bcftools_module}

        echo "Chromosomes processed: {params.n_chr}" >> {output.done_marker}

        echo "" >> {output.done_marker}    log:        module load {params.plink_module}

        echo "QC Filters Applied:" >> {output.done_marker}

        echo "- INFO score threshold: {params.info_score}" >> {output.done_marker}        f"results/log/convert_to_plink_chr{{chr}}_{DATASET_NAME}.log"        

        echo "- MAF threshold: {params.maf}" >> {output.done_marker}

        echo "- Missing rate threshold: {params.missing_rate}" >> {output.done_marker}    threads: 2        # Create output directory

        echo "" >> {output.done_marker}

            resources:        mkdir -p {PLINK_DIR}

        # Verify all chromosomes have complete PLINK files

        echo "File Verification:" >> {output.done_marker}        mem_mb = 32000,        mkdir -p results/log

        COMPLETE_COUNT=0

                time = "00:30:00"        

        for chr in {{1..{params.n_chr}}}; do

            if [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bed ]] && \\    shell:        # Start logging

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bim ]] && \\

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.fam ]]; then        """        echo "Starting VCF to PLINK conversion for chromosome {params.chromosome}" > {log}

                echo "Chromosome ${{chr}}: Complete (.bed/.bim/.fam)" >> {output.done_marker}

                COMPLETE_COUNT=$((COMPLETE_COUNT + 1))        # Load required software module        echo "Dataset: {params.dataset}" >> {log}

            else

                echo "Chromosome ${{chr}}: ERROR - Missing PLINK files" >> {output.done_marker}        module load {params.plink1_module}        echo "Loaded modules: {params.bcftools_module}, {params.plink_module}" >> {log}

            fi

        done                echo "Started at: $(date)" >> {log}

        

        echo "" >> {output.done_marker}        # Create output directory        echo "" >> {log}

        echo "Summary: $COMPLETE_COUNT/{params.n_chr} chromosomes successfully converted" >> {output.done_marker}

        echo "Output directory: {PLINK_FORMAT_DIR}" >> {output.done_marker}        mkdir -p {PLINK_FORMAT_DIR}        

        echo "Ready for downstream analysis or merging." >> {output.done_marker}

                mkdir -p results/log        # Use helper script for conversion

        echo "PLINK conversion finalization completed at: $(date)" >> {log}

        echo "Completion marker: {output.done_marker}" >> {log}                bash scripts/convert_chromosome_to_plink.sh \\

        echo "Successfully converted: $COMPLETE_COUNT/{params.n_chr} chromosomes" >> {log}

        """        # Start logging            {input.filtered_vcf} \\

        echo "Starting VCF to PLINK1 conversion for chromosome {params.chromosome}" > {log}            {output.clean_vcf} \\

        echo "Dataset: {params.dataset}" >> {log}            {PLINK_DIR}/chr{params.chromosome} \\

        echo "Loaded module: {params.plink1_module}" >> {log}            {params.chromosome} \\

        echo "QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}            {threads} 2>&1 | tee -a {log}

        echo "Started at: $(date)" >> {log}        

        echo "" >> {log}        echo "" >> {log}

                echo "Conversion completed at: $(date)" >> {log}

        # Use helper script for conversion        """

        bash scripts/convert_to_plink.sh \\

            {input.plink_ready_vcf} \\rule convert_vcf_to_plink_all:

            {params.plink_prefix} \\    """

            {output.conversion_log} \\    Convert all chromosomes from VCF to PLINK1 format in parallel.

            {params.chromosome} \\    

            {params.info_score} \\    This rule coordinates the conversion of all autosomes (chromosomes 1-22)

            {params.maf} \\    and creates a completion marker when all conversions are done.

            {params.missing_rate} \\    """

            {threads} 2>&1 | tee -a {log}    input:

                # Wait for all individual chromosome conversions to complete

        echo "" >> {log}        expand(f"{PLINK_DIR}/chr{{chr}}.bed", chr=range(1, N_CHROMOSOMES + 1)),

        echo "PLINK conversion completed at: $(date)" >> {log}        expand(f"{PLINK_DIR}/chr{{chr}}.bim", chr=range(1, N_CHROMOSOMES + 1)),

        """        expand(f"{PLINK_DIR}/chr{{chr}}.fam", chr=range(1, N_CHROMOSOMES + 1)),

        expand(f"{PLINK_DIR}/chr{{chr}}.conversion_summary.txt", chr=range(1, N_CHROMOSOMES + 1))

rule convert_to_plink_all:    output:

    """        # Completion marker and summary report

    Convert all chromosomes from VCF to PLINK1 binary format.        done_marker = f"{PLINK_DIR}/plink_conversion_complete.done",

            summary_report = f"{PLINK_DIR}/conversion_summary_report.txt"

    This rule coordinates the conversion of all autosomes (chromosomes 1-22)    params:

    and creates a completion marker when all conversions are done.        dataset = DATASET_NAME,

    """        n_chr = N_CHROMOSOMES

    input:    log:

        # Wait for all individual chromosome conversions to complete        f"results/log/convert_plink_all_{DATASET_NAME}.log"

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bed", chr=range(1, N_CHROMOSOMES + 1)),    resources:

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.bim", chr=range(1, N_CHROMOSOMES + 1)),        mem_mb = 32000,

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.fam", chr=range(1, N_CHROMOSOMES + 1)),        time = "00:30:00"

        expand(f"{PLINK_FORMAT_DIR}/chr{{chr}}.conversion.log", chr=range(1, N_CHROMOSOMES + 1))    shell:

    output:        """

        # Completion marker        echo "Generating PLINK conversion summary report for dataset {params.dataset}" > {log}

        done_marker = f"{PLINK_FORMAT_DIR}/plink_conversion_complete.done"        echo "Processed {params.n_chr} chromosomes" >> {log}

    params:        echo "Started at: $(date)" >> {log}

        dataset = DATASET_NAME,        

        n_chr = N_CHROMOSOMES,        # Create comprehensive summary report

        info_score = INFO_SCORE,        echo "PLINK Conversion Summary Report" > {output.summary_report}

        maf = MAF,        echo "==============================" >> {output.summary_report}

        missing_rate = MISSING_RATE        echo "Dataset: {params.dataset}" >> {output.summary_report}

    log:        echo "Date: $(date)" >> {output.summary_report}

        f"results/log/convert_to_plink_all_{DATASET_NAME}.log"        echo "Chromosomes processed: {params.n_chr}" >> {output.summary_report}

    resources:        echo "" >> {output.summary_report}

        mem_mb = 32000,        

        time = "00:30:00"        echo "Per-Chromosome Results:" >> {output.summary_report}

    shell:        echo "======================" >> {output.summary_report}

        """        

        echo "Finalizing PLINK conversion for dataset {params.dataset}" > {log}        # Aggregate statistics from all chromosomes

        echo "Processed {params.n_chr} chromosomes" >> {log}        TOTAL_VARIANTS=0

        echo "Applied QC filters - INFO: {params.info_score}, MAF: {params.maf}, Missing: {params.missing_rate}" >> {log}        TOTAL_SAMPLES=0

        echo "Started at: $(date)" >> {log}        

                for chr in {{1..{params.n_chr}}}; do

        # Create completion marker with summary information            if [[ -f {PLINK_DIR}/chr${{chr}}.conversion_summary.txt ]]; then

        echo "PLINK Conversion Complete" > {output.done_marker}                echo "Chromosome ${{chr}}:" >> {output.summary_report}

        echo "========================" >> {output.done_marker}                cat {PLINK_DIR}/chr${{chr}}.conversion_summary.txt | sed 's/^/  /' >> {output.summary_report}

        echo "Dataset: {params.dataset}" >> {output.done_marker}                echo "" >> {output.summary_report}

        echo "Date: $(date)" >> {output.done_marker}                

        echo "Chromosomes processed: {params.n_chr}" >> {output.done_marker}                # Extract variant count if available

        echo "" >> {output.done_marker}                if grep -q "Final variants in PLINK:" {PLINK_DIR}/chr${{chr}}.conversion_summary.txt; then

        echo "QC Filters Applied:" >> {output.done_marker}                    CHR_VARIANTS=$(grep "Final variants in PLINK:" {PLINK_DIR}/chr${{chr}}.conversion_summary.txt | awk '{{print $NF}}')

        echo "- INFO score threshold: {params.info_score}" >> {output.done_marker}                    if [[ "$CHR_VARIANTS" =~ ^[0-9]+$ ]]; then

        echo "- MAF threshold: {params.maf}" >> {output.done_marker}                        TOTAL_VARIANTS=$((TOTAL_VARIANTS + CHR_VARIANTS))

        echo "- Missing rate threshold: {params.missing_rate}" >> {output.done_marker}                    fi

        echo "" >> {output.done_marker}                fi

                        

        # Verify all chromosomes have complete PLINK files                # Extract sample count from first chromosome

        echo "File Verification:" >> {output.done_marker}                if [[ $chr -eq 1 ]] && grep -q "Samples:" {PLINK_DIR}/chr${{chr}}.conversion_summary.txt; then

        COMPLETE_COUNT=0                    TOTAL_SAMPLES=$(grep "Samples:" {PLINK_DIR}/chr${{chr}}.conversion_summary.txt | awk '{{print $NF}}')

                        fi

        for chr in {{1..{params.n_chr}}}; do            else

            if [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bed ]] && \\                echo "WARNING: Missing conversion summary for chromosome ${{chr}}" >> {log}

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.bim ]] && \\                echo "Chromosome ${{chr}}: ERROR - Missing conversion summary" >> {output.summary_report}

               [[ -f {PLINK_FORMAT_DIR}/chr${{chr}}.fam ]]; then            fi

                echo "Chromosome ${{chr}}: Complete (.bed/.bim/.fam)" >> {output.done_marker}        done

                COMPLETE_COUNT=$((COMPLETE_COUNT + 1))        

            else        # Overall summary statistics

                echo "Chromosome ${{chr}}: ERROR - Missing PLINK files" >> {output.done_marker}        echo "Overall Summary:" >> {output.summary_report}

            fi        echo "===============" >> {output.summary_report}

        done        echo "Total variants across all chromosomes: $TOTAL_VARIANTS" >> {output.summary_report}

                echo "Total samples: $TOTAL_SAMPLES" >> {output.summary_report}

        echo "" >> {output.done_marker}        echo "" >> {output.summary_report}

        echo "Summary: $COMPLETE_COUNT/{params.n_chr} chromosomes successfully converted" >> {output.done_marker}        

        echo "Output directory: {PLINK_FORMAT_DIR}" >> {output.done_marker}        echo "Output Files:" >> {output.summary_report}

        echo "Ready for downstream analysis or merging." >> {output.done_marker}        echo "============" >> {output.summary_report}

                for chr in {{1..{params.n_chr}}}; do

        echo "PLINK conversion finalization completed at: $(date)" >> {log}            if [[ -f {PLINK_DIR}/chr${{chr}}.bed ]]; then

        echo "Completion marker: {output.done_marker}" >> {log}                echo "Chromosome ${{chr}}:" >> {output.summary_report}

        echo "Successfully converted: $COMPLETE_COUNT/{params.n_chr} chromosomes" >> {log}                echo "  PLINK files: {PLINK_DIR}/chr${{chr}}.bed/.bim/.fam" >> {output.summary_report}

        """                echo "  Clean VCF: {PLINK_DIR}/chr${{chr}}.clean.vcf.gz" >> {output.summary_report}
                echo "  Summary: {PLINK_DIR}/chr${{chr}}.conversion_summary.txt" >> {output.summary_report}
            else
                echo "Chromosome ${{chr}}: ERROR - Missing PLINK files" >> {output.summary_report}
            fi
        done
        
        # Create completion marker
        echo "PLINK conversion completed successfully" > {output.done_marker}
        echo "Dataset: {params.dataset}" >> {output.done_marker}
        echo "Chromosomes: {params.n_chr}" >> {output.done_marker}
        echo "Completion time: $(date)" >> {output.done_marker}
        echo "Total variants: $TOTAL_VARIANTS" >> {output.done_marker}
        echo "Total samples: $TOTAL_SAMPLES" >> {output.done_marker}
        
        echo "Conversion summary completed at: $(date)" >> {log}
        echo "Summary report: {output.summary_report}" >> {log}
        echo "Completion marker: {output.done_marker}" >> {log}
        """