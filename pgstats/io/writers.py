"""
Data writing utilities for various output formats.
"""

import sgkit as sg
import xarray as xr
import numpy as np
from typing import Optional, Dict, Any


def save_results(ds: xr.Dataset, 
                 output_path: str,
                 format: str = "zarr",
                 **kwargs) -> None:
    """
    Save analysis results to various formats.
    
    Args:
        ds: Dataset to save
        output_path: Output file path
        format: Output format ('zarr', 'csv', 'tsv')
        **kwargs: Additional arguments passed to writer
    """
    if format == "zarr":
        save_zarr(ds, output_path, **kwargs)
    elif format == "csv":
        save_csv(ds, output_path, **kwargs)
    elif format == "tsv":
        save_tsv(ds, output_path, **kwargs)
    else:
        raise ValueError(f"Unsupported output format: {format}")


def save_zarr(ds: xr.Dataset, output_path: str, **kwargs) -> None:
    """
    Save dataset in Zarr format using sgkit.
    
    Args:
        ds: Dataset to save
        output_path: Output Zarr path
        **kwargs: Additional arguments passed to sgkit.save_dataset
    """
    sg.save_dataset(ds, output_path, **kwargs)


def save_csv(ds: xr.Dataset, output_path: str, **kwargs) -> None:
    """
    Save dataset statistics to CSV format.
    
    Args:
        ds: Dataset to save
        output_path: Output CSV path
        **kwargs: Additional arguments
    """
    # Extract statistical variables
    stat_vars = []
    for var_name in ds.data_vars:
        if var_name in ['tajima_d', 'fu_li_d', 'fu_li_f', 'nucleotide_diversity', 
                       'watterson_theta', 'fay_wu_theta', 'hka_statistic', 'hka_p_value']:
            stat_vars.append(var_name)
    
    if not stat_vars:
        raise ValueError("No statistical variables found in dataset")
    
    # Create DataFrame
    import pandas as pd
    
    data = {}
    for var_name in stat_vars:
        if var_name in ds.data_vars:
            data[var_name] = ds[var_name].values
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False, **kwargs)


def save_tsv(ds: xr.Dataset, output_path: str, **kwargs) -> None:
    """
    Save dataset statistics to TSV format.
    
    Args:
        ds: Dataset to save
        output_path: Output TSV path
        **kwargs: Additional arguments
    """
    # Extract statistical variables
    stat_vars = []
    for var_name in ds.data_vars:
        if var_name in ['tajima_d', 'fu_li_d', 'fu_li_f', 'nucleotide_diversity', 
                       'watterson_theta', 'fay_wu_theta', 'hka_statistic', 'hka_p_value']:
            stat_vars.append(var_name)
    
    if not stat_vars:
        raise ValueError("No statistical variables found in dataset")
    
    # Create DataFrame
    import pandas as pd
    
    data = {}
    for var_name in stat_vars:
        if var_name in ds.data_vars:
            data[var_name] = ds[var_name].values
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False, **kwargs)


def save_plink(ds: xr.Dataset, output_prefix: str, **kwargs) -> None:
    """
    Save dataset in PLINK format using sgkit.
    
    Args:
        ds: Dataset to save
        output_prefix: Output file prefix (without extension)
        **kwargs: Additional arguments passed to sgkit.write_plink
    """
    sg.write_plink(ds, output_prefix, **kwargs)
