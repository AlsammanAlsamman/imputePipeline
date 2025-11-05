import sys
sys.path.append("utils")
from bioconfigme import get_vcf_convert_params

# Get VCF conversion parameters from config
vcf_params = get_vcf_convert_params()

rule vcf_convert:
    """
    Convert PLINK to VCF for all chromosomes and perform harmonization.
    Processes all chromosomes 1-22 in parallel and creates single done file.
    """
    input:
        done="results/02_remove_ambiguous/{dataset}_remove_ambiguous.done",
        bed="results/02_remove_ambiguous/{dataset}_no_ambiguous.bed",
        bim="results/02_remove_ambiguous/{dataset}_no_ambiguous.bim",
        fam="results/02_remove_ambiguous/{dataset}_no_ambiguous.fam"
    output:
        done="results/03_vcf/{dataset}_convertvcf.done",
        vcf=expand("results/03_vcf/{{dataset}}_chr{chr}.vcf.gz", chr=range(1, 23)),
        tbi=expand("results/03_vcf/{{dataset}}_chr{chr}.vcf.gz.tbi", chr=range(1, 23)),
        stats="results/03_vcf/{dataset}_convertvcf_stats.txt"
    log:
        "results/logs/{dataset}_convertvcf.log"
    resources:
        mem_mb=32000,
        cores=4,
        time="02:00:00"
    shell:
        """
        bash scripts/vcf_convert_harmonize.sh {vcf_params[reference_fasta]} {input.bed} {input.bim} {input.fam} {wildcards.dataset} {output.done} {output.stats} > {log} 2>&1
        """