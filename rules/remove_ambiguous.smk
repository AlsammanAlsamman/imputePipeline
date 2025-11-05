import sys
sys.path.append("utils")
from bioconfigme import get_harmonize_params

# Get harmonization parameters from config
harmonize_params = get_harmonize_params()

rule remove_ambiguous:
    """
    Remove ambiguous A/T and G/C SNPs that cause strand issues.
    First step in harmonization process.
    """
    input:
        done="results/01_pre_qc/{dataset}_input_qc.done",
        bed="results/01_pre_qc/{dataset}_qc.bed",
        bim="results/01_pre_qc/{dataset}_qc.bim",
        fam="results/01_pre_qc/{dataset}_qc.fam"
    output:
        done="results/02_remove_ambiguous/{dataset}_remove_ambiguous.done",
        bed="results/02_remove_ambiguous/{dataset}_no_ambiguous.bed",
        bim="results/02_remove_ambiguous/{dataset}_no_ambiguous.bim", 
        fam="results/02_remove_ambiguous/{dataset}_no_ambiguous.fam",
        stats="results/02_remove_ambiguous/{dataset}_remove_ambiguous_stats.txt"
    log:
        "results/logs/{dataset}_remove_ambiguous.log"
    resources:
        mem_mb=24000,
        cores=2,
        time="02:00:00"
    shell:
        """
        bash scripts/remove_ambiguous_filter.sh {harmonize_params[reference_fasta]} {input.bed} {input.bim} {input.fam} {wildcards.dataset} {output.bed} {output.bim} {output.fam} {output.stats} {output.done} > {log} 2>&1
        """