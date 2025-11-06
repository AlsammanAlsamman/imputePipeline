#!/usr/bin/env python3
"""
bioconfigme.py - Simple YAML configuration loader for imputation pipeline

This module provides basic functionality to load and access YAML configuration
files for the genomic imputation pipeline. It's designed to be lightweight and
easily extensible.

Functions:
    load_config(file_path): Load a single YAML config file
    get(config, key, default=None): Get a value from config dict
    
Usage:
    from utils.bioconfigme import load_config, get
    
    # Load configuration
    software_config = load_config("configs/software.yml")
    analysis_config = load_config("configs/analysis.yml")
    
    # Access values
    plink_module = get(software_config, "plink.module")
    maf_threshold = get(analysis_config, "qc_filters.post_imputation.maf")
"""

import yaml
import os
from typing import Any, Dict, Optional


def load_config(file_path: str) -> Dict[str, Any]:
    """
    Load a YAML configuration file.
    
    Args:
        file_path (str): Path to the YAML configuration file
        
    Returns:
        Dict[str, Any]: Parsed YAML configuration as dictionary
        
    Raises:
        FileNotFoundError: If the config file doesn't exist
        yaml.YAMLError: If the YAML file is malformed
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Configuration file not found: {file_path}")
    
    try:
        with open(file_path, 'r') as file:
            config = yaml.safe_load(file)
            return config if config is not None else {}
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Error parsing YAML file {file_path}: {e}")


def get(config: Dict[str, Any], key: str, default: Any = None) -> Any:
    """
    Get a value from configuration dictionary using dot notation.
    
    Args:
        config (Dict[str, Any]): Configuration dictionary
        key (str): Key path using dot notation (e.g., "plink.module")
        default (Any): Default value to return if key not found
        
    Returns:
        Any: Value at the specified key path, or default if not found
        
    Examples:
        >>> config = {"plink": {"module": "plink2/1.90b3w"}}
        >>> get(config, "plink.module")
        'plink2/1.90b3w'
        >>> get(config, "nonexistent.key", "default_value")
        'default_value'
    """
    keys = key.split('.')
    current = config
    
    for k in keys:
        if isinstance(current, dict) and k in current:
            current = current[k]
        else:
            return default
    
    return current


def load_all_configs(config_dir: str = "configs") -> Dict[str, Dict[str, Any]]:
    """
    Load all configuration files from the configs directory.
    
    Args:
        config_dir (str): Directory containing configuration files
        
    Returns:
        Dict[str, Dict[str, Any]]: Dictionary with config names as keys
        
    Example:
        >>> configs = load_all_configs()
        >>> plink_module = get(configs["software"], "plink.module")
    """
    configs = {}
    
    # Standard config files
    config_files = {
        "software": os.path.join(config_dir, "software.yml"),
        "analysis": os.path.join(config_dir, "analysis.yml")
    }
    
    for name, path in config_files.items():
        if os.path.exists(path):
            configs[name] = load_config(path)
        else:
            print(f"Warning: Configuration file not found: {path}")
            configs[name] = {}
    
    return configs


def get_qc_params() -> Dict[str, Any]:
    """
    Get QC parameters for input filtering.
    
    Returns:
        Dict[str, Any]: Dictionary with QC parameters and paths
    """
    analysis_config = load_config("configs/analysis.yml")
    
    return {
        "geno": get(analysis_config, "qc_filters.pre_imputation.geno", 0.05),
        "maf": get(analysis_config, "qc_filters.pre_imputation.maf", 0.01),
        "hwe": get(analysis_config, "qc_filters.pre_imputation.hwe", 1e-6),
        "input_plink": get(analysis_config, "dataset.input_plink", "")
    }


def get_harmonize_params() -> Dict[str, Any]:
    """
    Get harmonization parameters for ambiguous SNP removal.
    
    Returns:
        Dict[str, Any]: Dictionary with harmonization parameters and paths
    """
    analysis_config = load_config("configs/analysis.yml")
    
    return {
        "reference_fasta": get(analysis_config, "harmonization.reference_fasta", ""),
        "results_dir": get(analysis_config, "results_dir", "results")
    }


def get_vcf_convert_params() -> Dict[str, Any]:
    """
    Get VCF conversion and harmonization parameters.
    
    Returns:
        Dict[str, Any]: Dictionary with VCF conversion parameters and paths
    """
    analysis_config = load_config("configs/analysis.yml")
    
    return {
        "reference_fasta": get(analysis_config, "harmonization.reference_fasta", ""),
        "chromosomes": get(analysis_config, "dataset.chromosomes", 22),
        "results_dir": get(analysis_config, "results_dir", "results")
    }


def get_harmonize_params() -> Dict[str, Any]:
    """
    Get harmonization parameters from analysis config.
    
    Returns:
        Dict[str, Any]: Harmonization parameters including reference paths
    """
    analysis_config = load_config("configs/analysis.yml")
    
    params = {
        'dataset': get(analysis_config, 'dataset.name'),
        'results_dir': get(analysis_config, 'results_dir'),
        'reference_fasta': get(analysis_config, 'harmonization.reference_fasta'),
        'build': get(analysis_config, 'dataset.build'),
    }
    
    return params


def get_download_params() -> Dict[str, Any]:
    """
    Get download parameters from analysis config.
    
    Returns:
        Dict[str, Any]: Download parameters including links file path
    """
    analysis_config = load_config("configs/analysis.yml")
    
    params = {
        'dataset': get(analysis_config, 'dataset.name'),
        'results_dir': get(analysis_config, 'results_dir'),
        'links_file': get(analysis_config, 'dataset.imputation_results_links'),
    }
    
    return params


def get_extract_params() -> Dict[str, Any]:
    """
    Get extraction parameters from analysis config.
    
    Returns:
        Dict[str, Any]: Extraction parameters including password file path
    """
    analysis_config = load_config("configs/analysis.yml")
    
    results_dir = get(analysis_config, 'results_dir', 'results')
    
    params = {
        'dataset': get(analysis_config, 'dataset.name'),
        'results_dir': results_dir,
        'pass_file': get(analysis_config, 'dataset.pass_file'),
        'extracted_dir': f"{results_dir}/06_extracted_data",  # Constructed from results_dir
        'chromosomes': get(analysis_config, 'dataset.chromosomes', 22),
    }
    
    return params


def _require_config_value(name: str, value: Any) -> Any:
    """Ensure a config value exists, raising a clear error otherwise."""
    if value in (None, ""):
        raise KeyError(f"Required configuration value '{name}' is missing")
    return value


def get_plink_conversion_params() -> Dict[str, Any]:
    """Collect PLINK conversion settings with validation."""
    analysis_config = load_config("configs/analysis.yml")
    software_config = load_config("configs/software.yml")

    dataset = _require_config_value(
        "dataset.name",
        get(analysis_config, "dataset.name"),
    )

    results_dir = _require_config_value(
        "results_dir",
        get(analysis_config, "results_dir"),
    )

    plink_ready_dir = f"{results_dir}/08_plink_ready"
    plink_format_dir = f"{results_dir}/09_plink_format"
    log_dir = f"{results_dir}/log"

    raw_chromosomes = _require_config_value(
        "dataset.chromosomes",
        get(analysis_config, "dataset.chromosomes"),
    )

    if isinstance(raw_chromosomes, int):
        chromosomes = [str(c) for c in range(1, raw_chromosomes + 1)]
    elif isinstance(raw_chromosomes, (list, tuple)):
        chromosomes = [str(c) for c in raw_chromosomes]
    else:
        raise ValueError("dataset.chromosomes must be integer or list")

    info_score = _require_config_value(
        "qc_filters.post_imputation.info_score",
        get(analysis_config, "qc_filters.post_imputation.info_score"),
    )

    maf = _require_config_value(
        "qc_filters.post_imputation.maf",
        get(analysis_config, "qc_filters.post_imputation.maf"),
    )

    missing_rate = _require_config_value(
        "qc_filters.post_imputation.missing_rate",
        get(analysis_config, "qc_filters.post_imputation.missing_rate"),
    )

    plink1_module = _require_config_value(
        "software.plink1.module",
        get(software_config, "plink1.module"),
    )

    return {
        "dataset": dataset,
        "chromosomes": chromosomes,
        "results_dir": results_dir,
        "input_dir": plink_ready_dir,
        "output_dir": plink_format_dir,
        "log_dir": log_dir,
        "info_score": info_score,
        "maf": maf,
        "missing_rate": missing_rate,
        "plink1_module": plink1_module,
    }


def get_plink_merge_params() -> Dict[str, Any]:
    """Get parameters required for merging per-chromosome PLINK datasets."""
    analysis_config = load_config("configs/analysis.yml")
    software_config = load_config("configs/software.yml")

    dataset = _require_config_value(
        "dataset.name",
        get(analysis_config, "dataset.name"),
    )

    results_dir = _require_config_value(
        "results_dir",
        get(analysis_config, "results_dir"),
    )

    raw_chromosomes = _require_config_value(
        "dataset.chromosomes",
        get(analysis_config, "dataset.chromosomes"),
    )

    if isinstance(raw_chromosomes, int):
        chromosomes = [str(c) for c in range(1, raw_chromosomes + 1)]
    elif isinstance(raw_chromosomes, (list, tuple)):
        chromosomes = [str(c) for c in raw_chromosomes]
    else:
        raise ValueError("dataset.chromosomes must be integer or list")

    plink_module = _require_config_value(
        "software.plink2.module",
        get(software_config, "plink2.module"),
    )

    input_dir = f"{results_dir}/09_plink_format"
    output_dir = f"{results_dir}/10_plink_merged"
    log_dir = f"{results_dir}/log"

    return {
        "dataset": dataset,
        "chromosomes": chromosomes,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "log_dir": log_dir,
        "plink_module": plink_module,
    }


def get_software_module(software_name: str) -> str:
    """
    Get the module name for a specific software tool.
    
    Args:
        software_name (str): Name of the software (e.g., "bcftools", "plink", "plink2")
        
    Returns:
        str: Module name to load (e.g., "bcftools/1.15")
        
    Raises:
        KeyError: If software is not found in configuration
        
    Examples:
        >>> get_software_module("bcftools")
        'bcftools/1.15'
        >>> get_software_module("plink")
        'plink2/1.90b3w'
    """
    software_config = load_config("configs/software.yml")
    module_name = get(software_config, f"{software_name}.module")
    
    if module_name is None:
        raise KeyError(f"Software '{software_name}' not found in software configuration")
    
    return module_name


def get_software_params(software_name: str) -> Dict[str, Any]:
    """
    Get the parameters for a specific software tool.
    
    Args:
        software_name (str): Name of the software (e.g., "bcftools", "plink")
        
    Returns:
        Dict[str, Any]: Dictionary with software parameters
        
    Examples:
        >>> get_software_params("bcftools")
        {'output_type': '-O z', 'threads': '--threads 2'}
    """
    software_config = load_config("configs/software.yml")
    params = get(software_config, f"{software_name}.params", {})
    return params


def get_config_value(config_file: str, key_path: str, default: Any = None) -> Any:
    """
    Get a configuration value from a specific config file using dot notation.
    
    This is a convenience function that combines load_config and get.
    
    Args:
        config_file (str): Path to the configuration file
        key_path (str): Key path using dot notation (e.g., "dataset.name")
        default (Any): Default value to return if key not found
        
    Returns:
        Any: Value at the specified key path, or default if not found
        
    Examples:
        >>> get_config_value("configs/analysis.yml", "dataset.name")
        'MEX123'
        >>> get_config_value("configs/software.yml", "bcftools.module")
        'bcftools/1.15'
    """
    config = load_config(config_file)
    return get(config, key_path, default)


# TODO: Add these functions as the pipeline develops
def merge_configs(*configs: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge multiple configuration dictionaries.
    Later configurations override earlier ones.
    
    This function will be implemented when needed for complex config merging.
    """
    pass


def validate_config(config: Dict[str, Any], schema: Dict[str, Any]) -> bool:
    """
    Validate configuration against a schema.
    
    This function will be implemented when config validation is needed.
    """
    pass


def substitute_variables(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Substitute variables in configuration values (e.g., {base_dir}/data).
    
    This function will be implemented when variable substitution is needed.
    """
    pass


if __name__ == "__main__":
    # Simple test of the module
    print("Testing bioconfigme module...")
    
    # Test loading configs if they exist
    try:
        software_config = load_config("configs/software.yml")
        print("✓ Successfully loaded software.yml")
        
        analysis_config = load_config("configs/analysis.yml") 
        print("✓ Successfully loaded analysis.yml")
        
        # Test getting values
        plink_module = get(software_config, "plink.module", "not found")
        print(f"✓ PLINK module: {plink_module}")
        
        maf_threshold = get(analysis_config, "qc_filters.post_imputation.maf", 0.01)
        print(f"✓ MAF threshold: {maf_threshold}")
        
    except Exception as e:
        print(f"✗ Error testing module: {e}")