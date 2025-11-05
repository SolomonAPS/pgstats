"""
Data validation utilities for genetic datasets.
"""

import xarray as xr
import numpy as np
from typing import List, Dict, Any


def validate_dataset(ds: xr.Dataset, required_vars: List[str] = None) -> Dict[str, Any]:
    """
    Validate that a dataset contains required variables and has proper structure.
    
    Args:
        ds: xarray Dataset to validate
        required_vars: List of required variable names
        
    Returns:
        Dictionary with validation results
    """
    if required_vars is None:
        required_vars = ["call_genotype", "variant_allele"]
    
    results = {
        "valid": True,
        "errors": [],
        "warnings": [],
        "info": {}
    }
    
    # Check required variables
    for var in required_vars:
        if var not in ds.data_vars:
            results["valid"] = False
            results["errors"].append(f"Required variable '{var}' not found")
    
    # Check dimensions
    if "variants" not in ds.dims:
        results["valid"] = False
        results["errors"].append("Required dimension 'variants' not found")
    
    if "samples" not in ds.dims:
        results["valid"] = False
        results["errors"].append("Required dimension 'samples' not found")
    
    # Check genotype data
    if "call_genotype" in ds.data_vars:
        gt = ds["call_genotype"]
        if "ploidy" not in gt.dims:
            results["warnings"].append("Genotype data missing 'ploidy' dimension")
        
        # Check for missing data
        missing_count = np.sum(gt.values < 0)
        if missing_count > 0:
            results["info"]["missing_genotypes"] = int(missing_count)
            results["warnings"].append(f"Found {missing_count} missing genotypes")
    
    return results


def check_missing_data(ds: xr.Dataset, call_genotype: str = "call_genotype") -> Dict[str, Any]:
    """
    Check for missing data patterns in the dataset.
    
    Args:
        ds: xarray Dataset
        call_genotype: Name of the genotype variable
        
    Returns:
        Dictionary with missing data statistics
    """
    if call_genotype not in ds.data_vars:
        return {"error": f"Variable '{call_genotype}' not found"}
    
    gt = ds[call_genotype].values
    n_variants, n_samples, ploidy = gt.shape
    
    # Count missing data
    missing_per_variant = np.sum(gt < 0, axis=(1, 2))
    missing_per_sample = np.sum(gt < 0, axis=(0, 2))
    
    results = {
        "total_missing": int(np.sum(gt < 0)),
        "missing_per_variant": missing_per_variant.tolist(),
        "missing_per_sample": missing_per_sample.tolist(),
        "variants_with_missing": int(np.sum(missing_per_variant > 0)),
        "samples_with_missing": int(np.sum(missing_per_sample > 0)),
        "missing_rate": float(np.sum(gt < 0) / (n_variants * n_samples * ploidy))
    }
    
    return results
