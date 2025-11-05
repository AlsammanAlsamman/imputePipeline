import sys
sys.path.append("utils")
from bioconfigme import get_qc_params

# Get QC parameters from config
qc_params = get_qc_params()

rule input_qc:
    """
    Input quality control filtering of PLINK files.
    Applies basic QC filters: missing data, MAF, and HWE.
    """
    output:
        done="results/01_pre_qc/{dataset}_input_qc.done",
        bed="results/01_pre_qc/{dataset}_qc.bed",
        bim="results/01_pre_qc/{dataset}_qc.bim", 
        fam="results/01_pre_qc/{dataset}_qc.fam",
        stats="results/01_pre_qc/{dataset}_qc_stats.txt"
    log:
        "results/logs/{dataset}_input_qc.log"
    resources:
        mem_mb=32000,
        cores=2,
        time="00:30:00"
    shell:
        """
        bash scripts/input_qc_filter.sh {qc_params[geno]} {qc_params[maf]} {qc_params[hwe]} {qc_params[input_plink]} {wildcards.dataset} {output.bed} {output.bim} {output.fam} {output.stats} {output.done} > {log} 2>&1
        """