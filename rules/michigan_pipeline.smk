"""
Master Snakemake workflow for Michigan Imputation Server processing.

This workflow combines:
1. Extraction of password-protected imputation results
2. Quality filtering based on Michigan .info.gz files

Author: Genomics Pipeline  
Date: November 2025
"""

# Include individual rule files
include: "extract_imputed.smk"
include: "filter_imputed_michigan.smk"

# Master rule that runs the complete Michigan pipeline
rule michigan_pipeline_complete:
    """
    Complete Michigan Imputation Server processing pipeline.
    
    This rule coordinates:
    1. Extraction of all chromosome files
    2. Quality filtering of extracted files
    3. Generation of final summary reports
    """
    input:
        # Ensure extraction completes first
        extraction_done = "results/05_imputed_data/imputed_extraction_complete.done",
        # Then ensure filtering completes
        filtering_done = "results/07_filtered/michigan_filtering_complete.done"
    output:
        # Final completion marker
        pipeline_complete = "results/michigan_pipeline_complete.done"
    log:
        "logs/michigan_pipeline_complete.log"
    shell:
        """
        echo "Michigan Imputation Server pipeline completed successfully!" > {log}
        echo "Completion time: $(date)" >> {log}
        echo "" >> {log}
        
        echo "Pipeline Summary:" >> {log}
        echo "================" >> {log}
        
        # Extraction summary
        if [[ -f {input.extraction_done} ]]; then
            echo "✓ Extraction completed successfully" >> {log}
            echo "  Completion marker: {input.extraction_done}" >> {log}
        else
            echo "✗ Extraction failed" >> {log}
        fi
        
        # Filtering summary  
        if [[ -f {input.filtering_done} ]]; then
            echo "✓ Filtering completed successfully" >> {log}
            echo "  Completion marker: {input.filtering_done}" >> {log}
        else
            echo "✗ Filtering failed" >> {log}
        fi
        
        echo "" >> {log}
        echo "Output directories:" >> {log}
        echo "  Extracted files: results/05_imputed_data/" >> {log}
        echo "  Filtered files: results/07_filtered/" >> {log}
        
        # Create final completion marker
        echo "Michigan Imputation Server Pipeline - COMPLETE" > {output.pipeline_complete}
        echo "Extraction: DONE" >> {output.pipeline_complete}
        echo "Filtering: DONE" >> {output.pipeline_complete}
        echo "Completion time: $(date)" >> {output.pipeline_complete}
        
        echo "Pipeline completed successfully at: $(date)" >> {log}
        """