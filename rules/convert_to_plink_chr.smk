import os
sys.path.append("utils")
from bioconfigme import get_plink_conversion_params

cfg = get_plink_conversion_params()

DATASET = cfg["dataset"]
CHROMOSOMES = cfg["chromosomes"]
INPUT_DIR = cfg["input_dir"]
OUTPUT_DIR = cfg["output_dir"]
LOG_DIR = cfg["log_dir"]
INFO_SCORE = cfg["info_score"]
MAF_THRESHOLD = cfg["maf"]
MISSING_RATE = cfg["missing_rate"]
PLINK1_MODULE = cfg["plink1_module"]

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

wildcard_constraints:
    chr="|".join(CHROMOSOMES)

rule convert_to_plink_chr:
    """Convert chromosome-level VCF to PLINK1 format and record completion."""
    input:
        vcf=f"{INPUT_DIR}/chr{{chr}}.plink_ready.vcf.gz",
        index=f"{INPUT_DIR}/chr{{chr}}.plink_ready.vcf.gz.csi"
    output:
        done=f"{OUTPUT_DIR}/chr{{chr}}.done"
    params:
        dataset=DATASET,
        plink_prefix=f"{OUTPUT_DIR}/chr{{chr}}",
        conversion_log=f"{OUTPUT_DIR}/chr{{chr}}.conversion.log",
        info_score=INFO_SCORE,
        maf=MAF_THRESHOLD,
        missing_rate=MISSING_RATE,
        module=PLINK1_MODULE
    log:
        os.path.join(LOG_DIR, f"convert_to_plink_chr{{chr}}_{DATASET}.log")
    threads: 2
    resources:
        mem_mb=32000,
        time="00:30:00"
    shell:
        """
        bash scripts/convert_to_plink_chr.sh \
            {params.dataset} \
            {wildcards.chr} \
            {input.vcf} \
            {params.plink_prefix} \
            {params.conversion_log} \
            {params.info_score} \
            {params.maf} \
            {params.missing_rate} \
            {threads} \
            {params.module} \
            {output.done}
        """
