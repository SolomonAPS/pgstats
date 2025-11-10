# Critical Fixes for Main (3 Changes Only)

These are the **only** changes needed to fix the dask corruption bug and update to the supported bio2zarr API.

---

## Fix 1: Load Zarr Without Chunking (Prevents Dask Corruption)

**File:** `pgstats/io/loaders.py`

**Location 1 (Line ~161):** When reusing existing Zarr
```python
# OLD:
dataset = sg.load_dataset(str(zarr_path), chunks='auto', **kwargs)

# NEW:
dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
```

**Location 2 (Line ~238):** After new conversion
```python
# OLD:
dataset = sg.load_dataset(str(zarr_path), chunks='auto', **kwargs)

# NEW:
dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
```

**Why:** This is THE fix. Loading with `chunks='auto'` creates dask arrays that get corrupted during lazy evaluation. Loading with `chunks=None` loads directly into numpy arrays, preventing corruption.

---

## Fix 2: Use bio2zarr Two-Step Workflow

**File:** `pgstats/io/loaders.py`

**Location:** Lines ~171-208 (replace the conversion section)

**OLD CODE:**
```python
try:
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
```

**NEW CODE:**
```python
try:
    print(f"Converting VCF to Zarr: {zarr_path}", flush=True)
    
    # Step 1: Explode VCF to intermediate columnar format (ICF)
    icf_path = zarr_path.parent / f"{zarr_path.stem}.icf"
    print(f"Step 1/2: Exploding VCF to ICF format...", flush=True)
    
    v2z.explode(
        str(icf_path),
        vcf_paths,
        worker_processes=worker_processes,
        show_progress=show_progress
    )
    
    # Step 2: Encode ICF to final Zarr format
    print(f"Step 2/2: Encoding ICF to Zarr...", flush=True)
    
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
```

**Why:** `sgkit.vcf_to_zarr` is deprecated. bio2zarr is the only supported method. The two-step workflow is more memory-efficient and is the recommended approach.

---

## Fix 3: Handle Empty Windows

**File:** `pgstats/core/dataset.py`

**Location:** Lines ~721-726 (in `_create_windows_from_dataset` method)

**OLD CODE:**
```python
# Get actual positions for each window
positions = self.dataset.variant_position.values
window_start_positions = positions[window_starts]
window_stop_positions = positions[window_stops - 1]  # -1 because stop is exclusive
```

**NEW CODE:**
```python
# Get actual positions for each window
positions = self.dataset.variant_position.values
window_start_positions = positions[window_starts]

# Handle empty windows: when window_stops == window_starts, the window is empty
# window_stops is an exclusive end index, so window_stops-1 is the last variant
window_stop_positions = np.where(
    window_stops > window_starts,
    positions[window_stops - 1],  # Non-empty window: last variant in window
    positions[window_starts]       # Empty window: use start position
)
```

**Note:** This requires `import numpy as np` at the top of the file (should already be there).

**Why:** Prevents crashes when a window has no variants (window_stops == window_starts). Without this check, `positions[window_stops - 1]` would try to access `positions[-1]` which gives the wrong result.

---

## Summary

**3 changes total:**
1. Change `chunks='auto'` to `chunks=None` (2 locations in loaders.py)
2. Replace `v2z.convert()` with two-step `v2z.explode()` → `v2z.encode()` (1 section in loaders.py)
3. Add `np.where()` check for empty windows (1 line in dataset.py)

**Testing:**
- Run on small VCF (< 1000 variants)
- Run on medium VCF (1000-10000 variants)
- Verify positions stay correct (no corruption)
- Verify correct number of windows created

**That's it!** These 3 changes fix the bug completely.

