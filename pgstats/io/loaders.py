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

# Try to import scikit-allel for lightweight VCF loading
try:
    import allel
    ALLEL_AVAILABLE = True
except ImportError:
    ALLEL_AVAILABLE = False


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
             max_variants_in_memory: int = 10000,
             **kwargs) -> xr.Dataset:
    """
    Load VCF data, using in-memory loading for small files or bio2zarr for large files.
    
    For small VCFs (<= max_variants_in_memory), uses sgkit.vcf_to_zarr with in-memory
    conversion which is much faster and uses minimal memory. For larger VCFs, falls back
    to bio2zarr's two-step workflow.
    
    Args:
        vcf_path: Path to VCF file (.vcf or .vcf.gz) or list of VCF files
        temp_dir: Directory for temporary files (default: system temp)
        keep_zarr: Whether to keep the intermediate Zarr files
        variants_chunk_size: Chunk size for variants dimension
        samples_chunk_size: Chunk size for samples dimension
        worker_processes: Number of worker processes for parallel conversion
        show_progress: Whether to show progress during conversion
        max_variants_in_memory: Max variants to load in-memory (default: 10000)
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
    
    # DEBUG: Print what we're actually getting
    print(f"DEBUG load_vcf: vcf_path type={type(vcf_path)}, value={vcf_path}", flush=True)
    print(f"DEBUG load_vcf: vcf_paths={vcf_paths}", flush=True)
    
    # Validate VCF files exist
    for path in vcf_paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"VCF file not found: {path}")
    
    # NOTE: sgkit.vcf_to_zarr was deprecated - bio2zarr is now the only option
    # We use bio2zarr's two-step workflow which is more memory-efficient than the old single-pass
    
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
        print(f"DEBUG: Initial vcf_name from path: {vcf_name}", flush=True)
        # Remove .gz/.bgz first if present
        if vcf_name.endswith('.gz') or vcf_name.endswith('.bgz'):
            vcf_name = vcf_name.rsplit('.', 1)[0]
            print(f"DEBUG: After removing .gz: {vcf_name}", flush=True)
        # Remove .vcf/.bcf extension
        if vcf_name.endswith('.vcf') or vcf_name.endswith('.bcf'):
            vcf_name = vcf_name.rsplit('.', 1)[0]
            print(f"DEBUG: After removing .vcf: {vcf_name}", flush=True)
        zarr_path = temp_dir / f"{vcf_name}.vcz"
        print(f"DEBUG: Final zarr_path: {zarr_path}", flush=True)
    else:
        zarr_path = temp_dir / "combined.vcz"
        print(f"DEBUG: Multiple VCFs, zarr_path: {zarr_path}", flush=True)
    
    # Check if we can reuse existing Zarr file
    if _is_zarr_up_to_date(zarr_path, vcf_paths):
        print(f"Reusing existing Zarr file: {zarr_path}", flush=True)
        print(f"Loading Zarr dataset...", flush=True)
        # Load without chunking to avoid dask corruption
        print(f"DEBUG: Loading without dask chunking to prevent corruption...", flush=True)
        dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
        
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
        # Use two-step workflow (explode → encode) for better memory efficiency
        # This is the "medium dataset" workflow recommended by bio2zarr docs
        print(f"Converting VCF to Zarr: {zarr_path}", flush=True)
        
        # Step 1: Explode VCF to intermediate columnar format (ICF)
        icf_path = zarr_path.parent / f"{zarr_path.stem}.icf"
        print(f"Step 1/2: Exploding VCF to ICF format...", flush=True)
        
        # bio2zarr.vcf.explode signature: explode(icf_path, vcfs, **kwargs)
        # Note: icf_path is FIRST, vcfs is SECOND
        v2z.explode(
            str(icf_path),
            vcf_paths,
            worker_processes=worker_processes,
            show_progress=show_progress
        )
        
        # Step 2: Encode ICF to final Zarr format
        print(f"Step 2/2: Encoding ICF to Zarr...", flush=True)
        
        # bio2zarr.vcf.encode expects: encode(icf_path, zarr_path, **kwargs)
        v2z.encode(
            str(icf_path),
            str(zarr_path),
            variants_chunk_size=variants_chunk_size,
            samples_chunk_size=samples_chunk_size,
            worker_processes=worker_processes,
            show_progress=show_progress
        )
        
        # Clean up intermediate ICF directory
        import shutil
        if icf_path.exists():
            shutil.rmtree(icf_path)
        
        print(f"VCF conversion complete", flush=True)
        
        # Save metadata for future cache validation
        _save_zarr_metadata(zarr_path, vcf_paths)
        
        # DEBUG: Check Zarr file directly before sgkit loads it
        print(f"DEBUG: Checking raw Zarr file before sgkit loading...", flush=True)
        import zarr
        try:
            z = zarr.open(str(zarr_path), mode='r')
            if 'variant_position' in z:
                raw_positions = z['variant_position'][:]
                print(f"DEBUG: Raw Zarr positions - min: {raw_positions.min()}, max: {raw_positions.max()}, dtype: {raw_positions.dtype}", flush=True)
                print(f"DEBUG: Raw Zarr first 10: {raw_positions[:10]}", flush=True)
                print(f"DEBUG: Raw Zarr last 10: {raw_positions[-10:]}", flush=True)
            else:
                print(f"DEBUG: variant_position not found in Zarr root, checking variant_POS...", flush=True)
                if 'variant_POS' in z:
                    raw_positions = z['variant_POS'][:]
                    print(f"DEBUG: Raw Zarr variant_POS - min: {raw_positions.min()}, max: {raw_positions.max()}", flush=True)
        except Exception as e:
            print(f"DEBUG: Error reading raw Zarr: {e}", flush=True)
        
        # Load Zarr dataset with sgkit
        # CRITICAL: Load without chunking to avoid dask corruption issues
        # For small-medium datasets, loading into memory is more reliable than lazy dask evaluation
        print(f"Loading Zarr dataset...", flush=True)
        print(f"DEBUG: Loading without dask chunking to prevent corruption...", flush=True)
        dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
        
        # DEBUG: Check positions immediately after loading
        positions = dataset.variant_position.values
        print(f"DEBUG: Loaded positions - min: {positions.min()}, max: {positions.max()}, dtype: {positions.dtype}", flush=True)
        print(f"DEBUG: First 10 positions: {positions[:10]}", flush=True)
        print(f"DEBUG: Last 10 positions: {positions[-10:]}", flush=True)
        
        # CRITICAL FIX: Force dask arrays to be computed and cached immediately
        # This prevents corruption that occurs during lazy evaluation
        print(f"DEBUG: Computing and caching all arrays to prevent dask corruption...", flush=True)
        if hasattr(dataset.variant_position, 'compute'):
            print(f"DEBUG: variant_position is a dask array, computing...", flush=True)
            dataset = dataset.assign(variant_position=dataset.variant_position.compute())
        if hasattr(dataset.variant_contig, 'compute'):
            print(f"DEBUG: variant_contig is a dask array, computing...", flush=True)
            dataset = dataset.assign(variant_contig=dataset.variant_contig.compute())
        if hasattr(dataset.call_genotype, 'compute'):
            print(f"DEBUG: call_genotype is a dask array, computing...", flush=True)
            dataset = dataset.assign(call_genotype=dataset.call_genotype.compute())
        print(f"DEBUG: All critical arrays computed and cached", flush=True)
        
        # BUGFIX: Check for bio2zarr memory-related corruption
        # If bio2zarr runs out of memory during encode phase, variant_contig can contain invalid values
        if 'variant_contig' in dataset.data_vars:
            variant_contig_values = dataset.variant_contig.values
            if hasattr(variant_contig_values, 'compute'):
                variant_contig_values = variant_contig_values.compute()
            
            num_contigs = len(dataset.contig_id.values)
            if np.any(variant_contig_values < 0) or np.any(variant_contig_values >= num_contigs):
                print(f"ERROR: Detected corrupted variant_contig values!", flush=True)
                print(f"  This indicates bio2zarr ran out of memory during the encode phase.", flush=True)
                print(f"  Found {np.sum(variant_contig_values < 0)} negative values", flush=True)
                print(f"  Found {np.sum(variant_contig_values >= num_contigs)} values >= {num_contigs} (max should be {num_contigs-1})", flush=True)
                print(f"", flush=True)
                print(f"  SOLUTION: Increase memory allocation (currently insufficient)", flush=True)
                print(f"  For VCFs with applied missingness patterns, try 256GB or higher.", flush=True)
                raise ValueError("Corrupted variant_contig detected - insufficient memory during VCF->Zarr conversion")
        
        # BUGFIX: Check for corrupted variant_position values
        # Positions should be positive and in ascending order (or at least reasonable)
        if 'variant_position' in dataset.data_vars:
            variant_positions = dataset.variant_position.values
            if hasattr(variant_positions, 'compute'):
                # Only compute first and last few positions to check
                first_positions = variant_positions[:min(10, len(variant_positions))].compute()
                last_positions = variant_positions[-min(10, len(variant_positions)):].compute()
                
                # Check for obviously wrong values (negative, or > 1 billion which is larger than any chromosome)
                if np.any(first_positions < 0) or np.any(last_positions < 0):
                    print(f"ERROR: Detected negative variant positions!", flush=True)
                    print(f"  First positions: {first_positions}", flush=True)
                    print(f"  Last positions: {last_positions}", flush=True)
                    raise ValueError("Corrupted variant_position detected - negative values found")
                
                if np.any(first_positions > 1_000_000_000) or np.any(last_positions > 1_000_000_000):
                    print(f"ERROR: Detected impossibly large variant positions (>1 billion)!", flush=True)
                    print(f"  First positions: {first_positions}", flush=True)
                    print(f"  Last positions: {last_positions}", flush=True)
                    print(f"  This usually indicates bio2zarr ran out of memory.", flush=True)
                    print(f"  Please increase memory allocation and re-run.", flush=True)
                    raise ValueError("Corrupted variant_position detected - likely due to insufficient memory during VCF->Zarr conversion")
                
                # Check if positions are in descending order (first > last), which is wrong
                if first_positions[0] > last_positions[-1]:
                    print(f"ERROR: Variant positions are in wrong order!", flush=True)
                    print(f"  First position: {first_positions[0]}", flush=True)
                    print(f"  Last position: {last_positions[-1]}", flush=True)
                    print(f"  This usually indicates bio2zarr ran out of memory.", flush=True)
                    print(f"  Please increase memory allocation and re-run.", flush=True)
                    raise ValueError("Corrupted variant_position detected - positions in wrong order, likely due to insufficient memory during VCF->Zarr conversion")
        
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
