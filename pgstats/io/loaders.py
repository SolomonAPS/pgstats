"""
Data loading utilities for various genetic data formats.
"""

import os
import tempfile
import numpy as np
import sgkit as sg
import xarray as xr
import bio2zarr.vcf as v2z
from typing import Optional, Union, List
from pathlib import Path
import hashlib
import time


def _get_vcf_hash(vcf_paths: List[str]) -> str:
    """Generate a hash of VCF file(s) for cache validation."""
    hasher = hashlib.md5()
    for vcf_path in sorted(vcf_paths):
        path = Path(vcf_path)
        if path.exists():
            # Include file size and modification time
            stat = path.stat()
            hasher.update(f"{path.name}:{stat.st_size}:{stat.st_mtime}".encode())
    return hasher.hexdigest()


def _is_zarr_up_to_date(zarr_path: Path, vcf_paths: List[str]) -> bool:
    """Check if Zarr file is up-to-date with VCF file(s)."""
    if not zarr_path.exists():
        return False
    
    # Check if metadata file exists and contains VCF hash
    metadata_file = zarr_path / ".vcf_metadata"
    if not metadata_file.exists():
        return False
    
    try:
        with open(metadata_file, 'r') as f:
            stored_hash = f.read().strip()
        current_hash = _get_vcf_hash(vcf_paths)
        return stored_hash == current_hash
    except:
        return False


def _save_zarr_metadata(zarr_path: Path, vcf_paths: List[str]):
    """Save VCF metadata to Zarr directory for cache validation."""
    metadata_file = zarr_path / ".vcf_metadata"
    vcf_hash = _get_vcf_hash(vcf_paths)
    with open(metadata_file, 'w') as f:
        f.write(vcf_hash)


def load_vcf(vcf_path: Union[str, List[str]], 
             temp_dir: Optional[str] = None,
             keep_zarr: bool = False,
             output_dir: Optional[str] = None,
             variants_chunk_size: Optional[int] = None,
             samples_chunk_size: Optional[int] = None,
             worker_processes: int = 0,
             show_progress: bool = True,
             **kwargs) -> xr.Dataset:
    """
    Load VCF data using bio2zarr Python API.
    
    This function converts VCF to Zarr format using bio2zarr's Python API
    and then loads the resulting Zarr dataset with sgkit.
    
    Args:
        vcf_path: Path to VCF file (.vcf or .vcf.gz) or list of VCF files
        temp_dir: Directory for temporary files (default: system temp)
        keep_zarr: Whether to keep the intermediate Zarr files
        variants_chunk_size: Chunk size for variants dimension
        samples_chunk_size: Chunk size for samples dimension
        worker_processes: Number of worker processes for parallel conversion
        show_progress: Whether to show progress during conversion
        **kwargs: Additional arguments passed to sgkit.load_dataset
        
    Returns:
        sgkit Dataset
        
    Raises:
        FileNotFoundError: If VCF file doesn't exist
        ImportError: If bio2zarr is not installed
        RuntimeError: If conversion fails
    """
    
    # Handle single VCF path or list of paths
    if isinstance(vcf_path, str):
        vcf_paths = [vcf_path]
    else:
        vcf_paths = vcf_path
    
    # Validate VCF files exist
    for path in vcf_paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"VCF file not found: {path}")
    
    # Set up temporary directory
    if temp_dir is None:
        if keep_zarr:
            if output_dir:
                # If keeping Zarr files and output directory specified, use output directory
                temp_dir = output_dir
            else:
                # If keeping Zarr files but no output directory, use directory next to VCF
                vcf_dir = Path(vcf_paths[0]).parent
                temp_dir = vcf_dir
        else:
            # If not keeping Zarr files, use system temp
            temp_dir = tempfile.mkdtemp(prefix="vcf2zarr_")
    else:
        os.makedirs(temp_dir, exist_ok=True)
    
    temp_dir = Path(temp_dir)
    
    # Generate Zarr output path with proper extension handling
    if len(vcf_paths) == 1:
        # Handle .vcf.gz, .vcf, .bcf.gz, .bcf properly
        vcf_name = Path(vcf_paths[0]).name
        # Remove .gz/.bgz first if present
        if vcf_name.endswith('.gz') or vcf_name.endswith('.bgz'):
            vcf_name = vcf_name.rsplit('.', 1)[0]
        # Remove .vcf/.bcf extension
        if vcf_name.endswith('.vcf') or vcf_name.endswith('.bcf'):
            vcf_name = vcf_name.rsplit('.', 1)[0]
        zarr_path = temp_dir / f"{vcf_name}.vcz"
    else:
        zarr_path = temp_dir / "combined.vcz"
    
    # Check if we can reuse existing Zarr file
    if _is_zarr_up_to_date(zarr_path, vcf_paths):
        print(f"Reusing existing Zarr file: {zarr_path}", flush=True)
        print(f"Loading Zarr dataset...", flush=True)
        # Use auto-chunking for memory efficiency
        dataset = sg.load_dataset(str(zarr_path), chunks='auto', **kwargs)
        
        # BUGFIX: Cast variant_contig to int16 to prevent int8 overflow
        if 'variant_contig' in dataset.data_vars:
            if dataset.variant_contig.dtype == np.int8:
                print(f"Converting variant_contig from int8 to int16 to prevent overflow...", flush=True)
                dataset = dataset.assign(
                    variant_contig=dataset.variant_contig.astype(np.int16)
                )
        
        print(f"Dataset loaded: {len(dataset.variants)} variants, {len(dataset.samples)} samples", flush=True)
        return dataset
    
    try:
        # Convert VCF to Zarr using bio2zarr Python API
        print(f"Converting VCF to Zarr: {zarr_path}", flush=True)
        v2z.convert(
            vcf_paths,
            str(zarr_path),
            variants_chunk_size=variants_chunk_size,
            samples_chunk_size=samples_chunk_size,
            worker_processes=worker_processes,
            show_progress=show_progress
        )
        print(f"VCF conversion complete", flush=True)
        
        # Save metadata for future cache validation
        _save_zarr_metadata(zarr_path, vcf_paths)
        
        # Load Zarr dataset with sgkit (try without dask chunking first for reliability)
        print(f"Loading Zarr dataset...", flush=True)
        # Use auto-chunking for memory efficiency
        dataset = sg.load_dataset(str(zarr_path), chunks='auto', **kwargs)
        
        # BUGFIX: Cast variant_contig to int16 to prevent int8 overflow
        # bio2zarr encodes variant_contig as int8 by default, which overflows at 128 contigs
        # or can get corrupted during dask operations with many BED regions
        if 'variant_contig' in dataset.data_vars:
            if dataset.variant_contig.dtype == np.int8:
                print(f"Converting variant_contig from int8 to int16 to prevent overflow...", flush=True)
                dataset = dataset.assign(
                    variant_contig=dataset.variant_contig.astype(np.int16)
                )
        
        # Validate dataset loaded correctly (without triggering computation)
        print(f"Dataset loaded: {len(dataset.variants)} variants, {len(dataset.samples)} samples", flush=True)
        
        return dataset
        
    except Exception as e:
        raise RuntimeError(f"VCF conversion failed: {e}")
    
    finally:
        # Clean up Zarr files unless requested to keep them
        if not keep_zarr and zarr_path.exists():
            import shutil
            shutil.rmtree(zarr_path)


def load_vcf_simple(vcf_path: Union[str, List[str]], keep_zarr: bool = False, temp_dir: Optional[str] = None, output_dir: Optional[str] = None, **kwargs) -> xr.Dataset:
    """
    Simple VCF loading function with automatic cleanup.
    
    This is a convenience wrapper around load_vcf that automatically
    handles temporary files and cleanup.
    
    Args:
        vcf_path: Path to VCF file (.vcf or .vcf.gz) or list of VCF files
        keep_zarr: Whether to keep the intermediate Zarr files
        temp_dir: Custom directory for Zarr files (overrides default logic)
        output_dir: Output directory for analysis (used for Zarr location with --keep-zarr)
        **kwargs: Additional arguments passed to load_vcf
        
    Returns:
        sgkit Dataset
    """
    return load_vcf(vcf_path, keep_zarr=keep_zarr, temp_dir=temp_dir, output_dir=output_dir, **kwargs)


def load_zarr(zarr_path: str, **kwargs) -> xr.Dataset:
    """
    Load Zarr dataset using sgkit.
    
    Args:
        zarr_path: Path to Zarr dataset
        **kwargs: Additional arguments passed to sgkit
        
    Returns:
        sgkit Dataset
    """
    dataset = sg.load_dataset(zarr_path, **kwargs)
    
    # BUGFIX: Cast variant_contig to int16 to prevent int8 overflow
    if 'variant_contig' in dataset.data_vars:
        if dataset.variant_contig.dtype == np.int8:
            print(f"Converting variant_contig from int8 to int16 to prevent overflow...", flush=True)
            dataset = dataset.assign(
                variant_contig=dataset.variant_contig.astype(np.int16)
            )
    
    return dataset


def check_bio2zarr_available() -> bool:
    """
    Check if bio2zarr package is available.
    
    Returns:
        True if bio2zarr is available, False otherwise
    """
    try:
        import bio2zarr.vcf
        return True
    except ImportError:
        return False


def load_dataset(data_path: str, 
                 format: Optional[str] = None,
                 **kwargs) -> xr.Dataset:
    """
    Load genetic dataset from various formats.
    
    Args:
        data_path: Path to data file
        format: Data format ('zarr', 'vcf')
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
    elif format == 'vcf':
        return load_vcf_simple(data_path, **kwargs)
    else:
        raise ValueError(f"Unsupported format: {format}")
