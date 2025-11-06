import os
sys.path.append("utils")
from  bioconfigme import get_plink_merge_params, get_software_module

cfg = get_plink_merge_params()

DATASET = cfg["dataset"]
OUTPUT_DIR = cfg["output_dir"]
LOG_DIR = cfg["log_dir"]
PLINK2_MODULE = get_software_module("plink2")

MERGED_PREFIX = f"{OUTPUT_DIR}/{DATASET}_merged"
BED_PREFIX = f"{OUTPUT_DIR}/{DATASET}_bed"

os.makedirs(LOG_DIR, exist_ok=True)

rule convert_plink2_to_bed:
    """Convert merged PLINK2 dataset to PLINK1 binary format."""
    input:
        pgen=f"{MERGED_PREFIX}.pgen",
        pvar=f"{MERGED_PREFIX}.pvar",
        psam=f"{MERGED_PREFIX}.psam",
        done=f"{MERGED_PREFIX}.done"
    output:
        bed=f"{BED_PREFIX}.bed",
        bim=f"{BED_PREFIX}.bim",
        fam=f"{BED_PREFIX}.fam",
        done=f"{BED_PREFIX}.done"
    params:
        dataset=DATASET,
        module=PLINK2_MODULE,
        pfile_prefix=MERGED_PREFIX,
        out_prefix=BED_PREFIX
    log:
        os.path.join(LOG_DIR, f"convert_plink2_to_bed_{DATASET}.log")
    threads: 2
    resources:
        mem_mb=32000,
        time="00:30:00"
    shell:
        """
        bash scripts/convert_plink2_to_bed.sh \
            {params.dataset} \
            {params.module} \
            {params.pfile_prefix} \
            {params.out_prefix} \
            {log} \
            {output.bed} \
            {output.bim} \
            {output.fam} \
            {output.done}
        """
