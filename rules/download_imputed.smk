import sys
sys.path.append("utils")
from bioconfigme import get_download_params

# Get download parameters from config
download_params = get_download_params()

rule download_imputed:
    """
    Download imputed chromosome files from Michigan/TOPMed server.
    Downloads files in batches of 5 with validation.
    """
    output:
        done_file="results/05_imputed_data/{dataset}_download_imputed.done",
        stats="results/05_imputed_data/{dataset}_download_stats.txt"
    params:
        dataset=download_params['dataset'],
        results_dir=download_params['results_dir'],
        links_file=download_params['links_file']
    resources:
        mem_mb=8000,
        cores=5,
        time="02:00:00"
    log:
        "results/log/{dataset}_download_imputed.log"
    shell:
        """
        bash scripts/download_imputed_batch.sh \
            {params.dataset} \
            {params.links_file} \
            {params.results_dir} \
            {output.stats} \
            {output.done_file} \
            2>&1 | tee {log}
        """