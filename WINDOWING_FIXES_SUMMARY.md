# Windowing Fixes Summary

## Key Changes Made

### 1. **dataset.py** - 3 Bug Fixes

#### Line 933-936: Fixed variant counting
**Before:**
```python
window_starts = self.windowed_dataset.window_start.values  # Wrong: positions
window_stops = self.windowed_dataset.window_stop.values    # Wrong: positions
```

**After:**
```python
window_starts = self.windowed_dataset.window_start_idx.values  # Correct: indices
window_stops = self.windowed_dataset.window_stop_idx.values    # Correct: indices
```

#### Line 866-868: Fixed window filtering
Same fix as above - use `_idx` variables to get variant indices.

#### Line 786-787: Fixed boundary checking
Same fix as above - use `_idx` variables to get variant indices.

### 2. **sfs_statistics.py** - Fixed 5 Statistics

Fixed: `theta_pi`, `theta_w`, `theta_h`, `theta_l`, `tajima_d`

**Pattern Applied:**
```python
# OLD: Loop through variants
for i in range(n_variants):
    pi_values[i] = calculate_pi(variant_matrix[i:i+1, :])
result["theta_pi"] = (["variants"], pi_values)  # Wrong dimension

# NEW: Loop through windows
window_starts = ds.window_start_idx.values
window_stops = ds.window_stop_idx.values
pi_values = np.zeros(n_windows)

for w_idx in range(n_windows):
    window_variant_matrix = variant_matrix[window_starts[w_idx]:window_stops[w_idx], :]
    pi_values[w_idx] = calculate_pi(window_variant_matrix)
result["theta_pi"] = (["windows"], pi_values)  # Correct dimension
```

### 3. **ld_statistics.py** - Fixed window dimension

**Line 580-588:** Changed from `(['variants'], ...)` to `(['windows'], ...)` for all LD statistics.

## Expected Input/Output

### Input to Statistic Functions

```python
# Dataset passed to stat functions has:
ds.windows                              # Shape: (n_windows,)
ds.window_start_idx.values             # Variant indices: [0, 5, 10, ...]
ds.window_stop_idx.values              # Variant indices: [5, 10, 15, ...]
ds[call_genotype]                      # Shape: (n_variants, n_samples, ploidy)
```

### Output from Statistic Functions

```python
# Returns dataset with:
result[stat_name]                      # Dimension: ['windows'], not ['variants']
result[stat_name].shape                # (n_windows,)
```

### Flow in calculate_windowed_stats()

```python
# 1. Call statistic function
stat_ds = tajima_d(self.windowed_dataset)

# 2. Merge results
result_dataset = result_dataset.merge(stat_ds)

# 3. Verify dimensions
assert 'windows' in result_dataset.tajima_d.dims
assert result_dataset.tajima_d.shape == (n_windows,)
```

## Testing

Run `test_tajima_d_windowing.py` to verify:

1. ✅ Results have `['windows']` dimension (not `['variants']`)
2. ✅ Shape matches number of windows
3. ✅ Can convert to DataFrame without errors
4. ✅ Works with single genome-wide window
5. ✅ Works with multiple windows

## Still To Fix

- [ ] fu_li_d, fu_li_f, fu_li_d_unfolded, fu_li_f_unfolded
- [ ] zeng_e, fay_wu_h
- [ ] haplotype_diversity, garud_h_statistics
- [ ] Clean up debug statements

## Next Steps

1. Run test to verify tajima_d works
\|2. Apply same pattern to remaining statistics
3. Remove debug print statements
4. Final testing with all statistics


