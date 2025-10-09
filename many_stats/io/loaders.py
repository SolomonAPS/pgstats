"""
Data loading utilities for various genetic data formats.
"""

import sgkit as sg
import xarray as xr
from typing import Optional, Union


def load_vcf(vcf_path: str, **kwargs) -> xr.Dataset:
    """
    Load VCF data using sgkit (via bio2zarr).
    
    Args:
        vcf_path: Path to VCF file
        **kwargs: Additional arguments passed to sgkit
        
    Returns:
        sgkit Dataset
    """
    # Note: This would typically use bio2zarr for VCF files
    # For now, we'll provide a placeholder
    raise NotImplementedError(
        "VCF loading requires bio2zarr. "
        "Please convert VCF to Zarr format first using bio2zarr."
    )


def load_plink(bed_path: str, 
               bim_path: Optional[str] = None,
               fam_path: Optional[str] = None,
               **kwargs) -> xr.Dataset:
    """
    Load PLINK data using sgkit.
    
    Args:
        bed_path: Path to .bed file
        bim_path: Path to .bim file (optional, inferred from bed_path)
        fam_path: Path to .fam file (optional, inferred from bed_path)
        **kwargs: Additional arguments passed to sgkit
        
    Returns:
        sgkit Dataset
    """
    return sg.read_plink(bed_path, bim_path, fam_path, **kwargs)


def load_zarr(zarr_path: str, **kwargs) -> xr.Dataset:
    """
    Load Zarr dataset using sgkit.
    
    Args:
        zarr_path: Path to Zarr dataset
        **kwargs: Additional arguments passed to sgkit
        
    Returns:
        sgkit Dataset
    """
    return sg.load_dataset(zarr_path, **kwargs)


def load_dataset(data_path: str, 
                 format: Optional[str] = None,
                 **kwargs) -> xr.Dataset:
    """
    Load genetic dataset from various formats.
    
    Args:
        data_path: Path to data file
        format: Data format ('plink', 'zarr', 'vcf')
        **kwargs: Additional arguments passed to loader
        
    Returns:
        sgkit Dataset
    """
    if format is None:
        # Infer format from file extension
        if data_path.endswith('.zarr'):
            format = 'zarr'
        elif data_path.endswith('.bed'):
            format = 'plink'
        elif data_path.endswith('.vcf') or data_path.endswith('.vcf.gz'):
            format = 'vcf'
        else:
            raise ValueError(f"Could not infer format from path: {data_path}")
    
    if format == 'zarr':
        return load_zarr(data_path, **kwargs)
    elif format == 'plink':
        return load_plink(data_path, **kwargs)
    elif format == 'vcf':
        return load_vcf(data_path, **kwargs)
    else:
        raise ValueError(f"Unsupported format: {format}")
