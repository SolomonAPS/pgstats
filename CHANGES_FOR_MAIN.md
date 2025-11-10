# Critical Changes to Cherry-Pick to Main

This document categorizes all changes made during debugging into:
1. **CRITICAL FIXES** - Must be applied to main
2. **USEFUL FEATURES** - Should be considered for main
3. **DEBUG CODE** - Should be removed before merging to main
4. **EXTERNAL SCRIPTS** - Not part of pgstats core, can be kept in branch or separate repo

---

## 🔴 CRITICAL FIXES (Must Apply to Main)

### 1. **Fix Dask Corruption by Loading Without Chunking**
**Commit:** `43ed95c`
**Files:** `pgstats/io/loaders.py`

**Change:**
```python
# OLD (causes corruption):
dataset = sg.load_dataset(str(zarr_path), chunks='auto', **kwargs)

# NEW (prevents corruption):
dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
```

**Why Critical:** This is THE fix that resolved the multi-day corruption bug. Without this, dask's lazy evaluation corrupts position arrays during windowing operations.

**Lines to change:**
- Line 161: `chunks=None` instead of `chunks='auto'` (reuse path)
- Line 238: `chunks=None` instead of `chunks='auto'` (new conversion path)

**Additional changes in this commit:**
- Lines 247-256: Force computation of dask arrays (belt-and-suspenders approach)
- This is technically redundant with `chunks=None` but provides extra safety

---

### 2. **Fix bio2zarr API: Use Two-Step Workflow**
**Commits:** `2d220a5`, `f7613af`, `c5fb807`, `5cf148f`
**Files:** `pgstats/io/loaders.py`

**Change:**
```python
# OLD (deprecated single-step):
v2z.convert(vcf_paths, str(zarr_path), ...)

# NEW (two-step workflow):
v2z.explode(str(icf_path), vcf_paths, ...)  # Step 1
v2z.encode(str(icf_path), str(zarr_path), ...)  # Step 2
```

**Why Critical:** `sgkit.vcf_to_zarr` is deprecated. bio2zarr is the only supported method, and the two-step workflow is more memory-efficient.

**Lines affected:** 171-208 in `loaders.py`

---

### 3. **Fix Window Position Calculation for Empty Windows**
**Commit:** `acf1d38`, `016a0d3`
**Files:** `pgstats/core/dataset.py`

**Change:**
```python
# OLD (crashes on empty windows):
window_stop_positions = positions[window_stops - 1]

# NEW (handles empty windows):
window_stop_positions = np.where(
    window_stops > window_starts,
    positions[window_stops - 1],  # Non-empty window
    positions[window_starts]       # Empty window
)
```

**Why Critical:** Prevents crashes when windows have no variants.

**Lines affected:** 721-726 in `dataset.py`

---

### 4. **Add numpy Import**
**Commit:** `228c7cb`
**Files:** `pgstats/io/loaders.py`

**Change:**
```python
import numpy as np
```

**Why Critical:** Required for the int8 overflow fix and array operations.

**Line:** 7 in `loaders.py`

---

### 5. **Fix Indentation Errors**
**Commits:** `3a1c38a`, `ac530a4`, `3fb5641`
**Files:** `pgstats/stats/sfs_statistics.py`

**Why Critical:** These were syntax errors that would prevent the code from running.

---

## 🟡 USEFUL FEATURES (Consider for Main)

### 1. **Callable Sites Normalization for Theta Statistics**
**Commit:** `6c1eeae`
**Files:** `pgstats/core/dataset.py`, `pgstats/stats/sfs_statistics.py`

**What it does:** Normalizes theta_pi, theta_w, theta_h, theta_l by the actual number of callable sites per window (excluding non-callable regions).

**Why useful:** Provides more accurate theta estimates when using BED files to mask non-callable regions.

**Lines affected:**
- `dataset.py`: Lines 1276-1293 (normalization logic)
- `sfs_statistics.py`: Lines 67-71 (filter zero callable sites)

**Decision:** **KEEP** - This is scientifically correct and important for accurate statistics.

---

### 2. **Add callable_sites to Output**
**Commits:** `d4c61b4`, `6182d2c`
**Files:** `pgstats/core/dataset.py`

**What it does:** Includes `callable_sites` column in output CSV showing how many callable sites were in each window.

**Why useful:** Diagnostic information for understanding window quality.

**Lines affected:** 
- Line 1276: Store callable_sites in dataset
- Line 1669: Add to output columns

**Decision:** **KEEP** - Very useful diagnostic information, minimal overhead.

---

### 3. **Int8 Overflow Detection**
**Commit:** `b60c37e`, `358ec51`
**Files:** `pgstats/io/loaders.py`

**What it does:** Detects and fixes variant_contig int8 overflow by casting to int16.

**Why useful:** Prevents crashes when there are >127 contigs.

**Lines affected:** 163-168 in `loaders.py`

**Decision:** **KEEP** - Good defensive programming, minimal cost.

---

### 4. **Force Dask Array Computation**
**Commit:** Part of `43ed95c`
**Files:** `pgstats/core/dataset.py`

**What it does:** Forces computation of dask arrays in `get_summary()` method.

**Lines affected:** 356-363 in `dataset.py`

**Decision:** **KEEP** - Defensive programming that prevents potential issues.

---

## 🔵 DEBUG CODE (Remove Before Main)

### 1. **All DEBUG Print Statements**
**Commits:** `1567692`, `a4e739e`, `69bed00`
**Files:** `pgstats/io/loaders.py`, `pgstats/core/dataset.py`

**What to remove:**
- Lines 106-108 in `loaders.py`: `DEBUG load_vcf: vcf_path type=...`
- Lines 140-153 in `loaders.py`: `DEBUG: Initial vcf_name from path...`
- Lines 160, 215-231, 240-244, 246-256 in `loaders.py`: All DEBUG prints
- Lines 662-667, 679-684, 713-731 in `dataset.py`: All DEBUG prints

**Decision:** **REMOVE ALL** - These were for debugging only.

---

### 2. **Documentation Files**
**Files:** `BIO2ZARR_MEMORY_FIX.md`, `SMART_LOADING_FIX.md`

**Decision:** **REMOVE** - These were investigation notes, not needed in main.

---

### 3. **reinstall_pgstats.sh**
**File:** `reinstall_pgstats.sh`

**Decision:** **REMOVE** - This was a convenience script for the cluster.

---

## 🟢 EXTERNAL SCRIPTS (Keep in Branch or Separate Repo)

These are analysis/pipeline scripts, not part of pgstats core:

1. **Missingness Profile Scripts:**
   - `create_missingness_profile.py`
   - `apply_missingness_profile.py`
   - `create_missingness_profiles_array.slurm`

2. **BED File Creation Scripts:**
   - `create_noncallable_bed.sh`
   - `create_noncallable_bed_dpgp3.sh`
   - `create_noncallable_bed_dpgp3.slurm`
   - `create_noncallable_bed_merged.sh`

3. **Analysis/Diagnostic Scripts:**
   - `analyze_regions_quality.sh`
   - `analyze_regions_quality.slurm`
   - `check_ancestral_allele.slurm`
   - `check_callable_sites.sh`
   - `check_vcf_density.sh`
   - `compare_bed_files.sh`
   - `compare_vcf_coverage.sh`
   - `compare_vcf_coverage.slurm`
   - `diagnose_noncallable.sh`
   - `process_polarized_vcfs.slurm`

**Decision:** These belong in your DFEpos project, not in pgstats. They're project-specific workflows.

---

## 📋 RECOMMENDED MERGE STRATEGY

### Option 1: Clean Cherry-Pick (Recommended)
1. Create a new branch from `main`
2. Cherry-pick only the critical fixes (manually apply changes)
3. Add the useful features
4. Test thoroughly
5. Merge to main

### Option 2: Interactive Rebase
1. Create a new branch from current branch
2. Use `git rebase -i origin/main` to squash/edit commits
3. Remove all DEBUG commits
4. Clean up commit messages
5. Test and merge

---

## 🎯 MINIMAL CHANGES FOR CLEAN MAIN

If you want the absolute minimum changes to fix the bug:

### File: `pgstats/io/loaders.py`

**Change 1:** Add numpy import (line 7)
```python
import numpy as np
```

**Change 2:** Update bio2zarr workflow (lines 171-208)
```python
# Replace v2z.convert() with two-step workflow
icf_path = zarr_path.parent / f"{zarr_path.stem}.icf"
v2z.explode(str(icf_path), vcf_paths, ...)
v2z.encode(str(icf_path), str(zarr_path), ...)
# Clean up ICF
import shutil
if icf_path.exists():
    shutil.rmtree(icf_path)
```

**Change 3:** Load without chunking (2 locations)
```python
# Line 161 and 238:
dataset = sg.load_dataset(str(zarr_path), chunks=None, **kwargs)
```

### File: `pgstats/core/dataset.py`

**Change 1:** Fix empty window handling (lines 721-726)
```python
window_stop_positions = np.where(
    window_stops > window_starts,
    positions[window_stops - 1],
    positions[window_starts]
)
```

**That's it!** These 4 changes fix the corruption bug.

---

## ✅ TESTING CHECKLIST

Before merging to main, test:
- [ ] Small VCF (< 1000 variants)
- [ ] Medium VCF (1000-10000 variants)  
- [ ] With BED file masking
- [ ] Without BED file masking
- [ ] Multiple regions
- [ ] Windowed analysis
- [ ] Genome-wide analysis
- [ ] Verify no position corruption
- [ ] Verify correct number of windows
- [ ] Verify statistics are reasonable

