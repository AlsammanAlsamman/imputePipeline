import os

from snakemake.io import expand

from utils.bioconfigme import get_plink_conversion_params

include: "convert_to_plink_chr.smk"

cfg = get_plink_conversion_params()

DATASET = cfg["dataset"]
CHROMOSOMES = cfg["chromosomes"]
OUTPUT_DIR = cfg["output_dir"]
LOG_DIR = cfg["log_dir"]

os.makedirs(LOG_DIR, exist_ok=True)

rule convert_to_plink_all:
    """Aggregate per-chromosome PLINK conversions into a single completion marker."""
    input:
        expand(f"{OUTPUT_DIR}/chr{{chr}}.done", chr=CHROMOSOMES)
    output:
        done=f"{OUTPUT_DIR}/plink_conversion_complete.done"
    params:
        dataset=DATASET,
        chr_list=",".join(CHROMOSOMES)
    log:
        os.path.join(LOG_DIR, f"convert_to_plink_all_{DATASET}.log")
    threads: 2
    resources:
        mem_mb=32000,
        time="00:30:00"
    shell:
        """
        mkdir -p "$(dirname {output.done})"
        printf "Dataset: {params.dataset}\n" > {log}
        printf "Chromosomes: {params.chr_list}\n" >> {log}
        printf "Completed files:\n" >> {log}
        for path in {input}; do
            printf "%s\n" "$path" >> {log}
        done
        touch {output.done}
        """
