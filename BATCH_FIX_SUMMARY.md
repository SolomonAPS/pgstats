# Batch Window Fix Summary

## Functions Still Needing Fix (in sfs_statistics.py)

1. **fu_li_d** (line 1141) - returns `(["variants"], fu_li_d_values)`
2. **fu_li_f** (line 1271) - returns `(["variants"], fu_li_f_values)`  
3. **zeng_e** (line 1420) - returns `(["variants"], zeng_e_values)`
4. **singletons** (line 1492/1494) - returns `(["variants"], singleton_values)`
5. **fay_wu_h** (line 1596) - returns `(["variants"], fay_wu_h_values)`

## Pattern to Apply

All functions need the same fix as tajima_d:

**Replace:**
```python
# Calculate for each variant
for i in range(n_variants):
    stat_values[i] = calculate_stat(variant_matrix[i:i+1, :])
result["stat_name"] = (["variants"], stat_values)
```

**With:**
```python
# Get window information
n_windows = len(ds.windows)
window_starts = ds.window_start_idx.values
window_stops = ds.window_stop_idx.values

# Calculate for each window
stat_values = np.zeros(n_windows)
for w_idx in range(n_windows):
    window_start = window_starts[w_idx]
    window_stop = window_stops[w_idx]
    window_variant_matrix = variant_matrix[window_start:window_stop, :]
    stat_values[w_idx] = calculate_stat(window_variant_matrix)
    
result["stat_name"] = (["windows"], stat_values)
```

## Status
Ready to apply - all functions follow similar patterns to tajima_d which is already working.

