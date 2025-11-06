import os

from snakemake.io import expand

from utils.bioconfigme import get_plink_merge_params
from utils.bioconfigme import get_software_module

cfg = get_plink_merge_params()
plinkmodule = get_software_module("plink2")

DATASET = cfg["dataset"]
CHROMOSOMES = cfg["chromosomes"]
INPUT_DIR = cfg["input_dir"]
OUTPUT_DIR = cfg["output_dir"]
LOG_DIR = cfg["log_dir"]
PLINK2_MODULE = plinkmodule

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

rule merge_plink2:
    """Convert per-chromosome PLINK1 binaries to PLINK2 format and merge them."""
    input:
        bed=expand(f"{INPUT_DIR}/chr{{chr}}.bed", chr=CHROMOSOMES),
        bim=expand(f"{INPUT_DIR}/chr{{chr}}.bim", chr=CHROMOSOMES),
        fam=expand(f"{INPUT_DIR}/chr{{chr}}.fam", chr=CHROMOSOMES)
    output:
        pgen=f"{OUTPUT_DIR}/{DATASET}_merged.pgen",
        pvar=f"{OUTPUT_DIR}/{DATASET}_merged.pvar",
        psam=f"{OUTPUT_DIR}/{DATASET}_merged.psam",
        done=f"{OUTPUT_DIR}/{DATASET}_merged.done"
    params:
        dataset=DATASET,
        chromosome_list=" ".join(CHROMOSOMES),
        input_dir=INPUT_DIR,
        output_prefix=f"{OUTPUT_DIR}/{DATASET}_merged",
        module=PLINK2_MODULE
    log:
        os.path.join(LOG_DIR, f"merge_plink2_{DATASET}.log")
    threads: 2
    resources:
        mem_mb=32000,
        time="00:30:00"
    shell:
        """
        bash scripts/merge_plink2.sh \
            {params.dataset} \
            "{params.chromosome_list}" \
            {params.input_dir} \
            {params.output_prefix} \
            {log} \
            {params.module} \
            {output.pgen} \
            {output.pvar} \
            {output.psam} \
            {output.done}
        """