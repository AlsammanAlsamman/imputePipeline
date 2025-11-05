import sys
sys.path.append("utils")
from bioconfigme import get_extract_params

# Get extraction parameters from config
extract_params = get_extract_params()

rule extract_imputed_chr:
    """
    Extract password-protected imputed files for individual chromosomes.
    Processes each chromosome independently for parallel execution.
    """
    output:
        done_file="results/06_extracted_data/chr{chr}_extract_done.txt",
        extracted_files=directory("results/06_extracted_data/chr{chr}")
    params:
        dataset=extract_params['dataset'],
        input_dir=f"{extract_params['results_dir']}/05_imputed_data",
        pass_file=extract_params['pass_file']
    resources:
        mem_mb=32000,
        cores=2,
        time="00:30:00"
    log:
        "results/log/chr{chr}_extract_imputed.log"
    shell:
        """
        bash scripts/extract_imputed_chr.sh \
            {params.dataset} \
            {wildcards.chr} \
            {params.input_dir} \
            {output.extracted_files} \
            {params.pass_file} \
            {output.done_file} \
            2>&1 | tee {log}
        """

rule extract_imputed_all:
    """
    Aggregate rule to extract all chromosomes and create final done file.
    """
    input:
        done_files=expand("results/06_extracted_data/chr{chr}_extract_done.txt", chr=range(1, 23)),
        extracted_dirs=expand("results/06_extracted_data/chr{chr}", chr=range(1, 23))
    output:
        done_file="results/06_extracted_data/{dataset}_extract_imputed.done"
    params:
        dataset=extract_params['dataset'],
        chromosomes=extract_params['chromosomes']
    resources:
        mem_mb=32000,
        cores=2,
        time="00:30:00"
    log:
        "results/log/{dataset}_extract_aggregate.log"
    shell:
        """
        bash scripts/extract_aggregate.sh \
            {params.dataset} \
            {params.chromosomes} \
            {output.done_file} \
            2>&1 | tee {log}
        """