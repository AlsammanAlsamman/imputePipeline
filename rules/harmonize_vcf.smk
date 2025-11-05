import sys
sys.path.append("utils")
from bioconfigme import get_harmonize_params

# Get harmonization parameters from config
harmonize_params = get_harmonize_params()

rule harmonize_vcf_chr:
    """
    Harmonize VCF files against reference genome for individual chromosomes.
    Handles chromosome naming based on genome build (hg19 vs hg38).
    """
    input:
        vcf=f"results/03_vcf/{harmonize_params['dataset']}_chr{{chr}}.vcf.gz",
        tbi=f"results/03_vcf/{harmonize_params['dataset']}_chr{{chr}}.vcf.gz.tbi"
    output:
        vcf="results/04_harmonize/chr{chr}.vcf.gz",
        tbi="results/04_harmonize/chr{chr}.vcf.gz.tbi"
    params:
        reference_fasta=harmonize_params['reference_fasta'],
        build=harmonize_params['build'],
        dataset=harmonize_params['dataset']
    resources:
        mem_mb=16000,
        cores=2,
        time="01:00:00"
    log:
        "results/log/chr{chr}_harmonize.log"
    shell:
        """
        bash scripts/harmonize_vcf_chr.sh \
            {params.dataset} \
            {wildcards.chr} \
            {params.reference_fasta} \
            {params.build} \
            {input.vcf} \
            {output.vcf} \
            2>&1 | tee {log}
        """

rule harmonize_vcf_all:
    """
    Aggregate rule to harmonize all chromosomes and create summary.
    """
    input:
        vcf_files=expand("results/04_harmonize/chr{chr}.vcf.gz", chr=range(1, 23)),
        tbi_files=expand("results/04_harmonize/chr{chr}.vcf.gz.tbi", chr=range(1, 23))
    output:
        done_file="results/04_harmonize/{dataset}_harmonize.done",
        stats="results/04_harmonize/{dataset}_harmonization_stats.txt"
    params:
        build=harmonize_params['build']
    log:
        "results/log/{dataset}_harmonize_summary.log"
    shell:
        """
        bash scripts/harmonize_summary.sh \
            {wildcards.dataset} \
            {params.build} \
            {output.stats} \
            {output.done_file} \
            2>&1 | tee {log}
        """