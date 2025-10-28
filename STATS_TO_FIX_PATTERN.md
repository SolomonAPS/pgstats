# Remaining Stats to Fix - Exact Pattern

## Pattern (from tajima_d/fu_li_d)

Replace the loop section in each function:

**OLD:**
```python
# Calculate max sample size across all variants
max_n = 0
for i in range(n_variants):
    site_n = ...
    if site_n > max_n: max_n = site_n

# Calculate harmonic numbers once using max_n
a1 = calculate_a1(max_n)
...

# Calculate statistic for each variant
stat_values = np.zeros(n_vari relie ants)

for i in range(n_variants):
    # calculate on variant_matrix[i:i+1, :]
    stat_values[i] = ...

result["stat_name"] = (["variants"], stat_values)
```

**NEW:**
```python
# Get window information (windows are always defined)
n_windows = len(ds.windows)
window_starts = ds.window_start_idx.values
window_stops = ds.window_stop_idx.values

# Calculate statistic for each window
stat_values = np.zeros(n_windows)

for w_idx in range(n_windows):
    # Extract variants in this window
    window_start = window_starts[w_idx]
    window_stop = window_stops[w_idx]
    window_variant_matrix = variant_matrix[window_start:window_stop, :]
    
    # Calculate max_n for THIS window (may differ per window)
    window_max_n = 0
    for i in range(len(window_variant_matrix)):
        site_n = 0
        for j in range(n_samples):
            if window_variant_matrix[i, j] != -1:
                site_n += 1
        if site_n > window_max_n:
            window_max_n = site_n
    
    if window_max_n <= 1:
        stat_values[w_idx] = 0.0
        continue
    
    # Calculate harmonic numbers for THIS window
    a1 = calculate_a1(window_max_n)
    ...
    
    # Calculate statistic on window_variant_matrix (entire window)
    stat_values[w_idx] = ...

result["stat_name"] = (["windows"], stat_values)
```

## Functions Still Needing Fix

1. **fu_li_f** (line ~1256) - Similar to fu_li_d (DONE IN CURRENT BRANCH)
2. **zeng_e** (line ~1417) - Similar to tajima_d  
3. **fay_wu_h** (line ~1597) - Similar to tajima_d
4. **singletons** (line ~1492/1494) - Simple theta style

## Status

- ✅ tajima_d - Working
- ✅ fu_li_d - Working (just pushed)
- 🔄 fu_li_f - Needs fix (same as fu_li_d pattern)
- 🔄 zeng_e - Needs fix
- 🔄 fay_wu_h - Needs fix  
- 🔄 singletons - Needs fix

## Next Steps

Apply the pattern to each function. The key differences:
- Loop through `w_idx` instead of `i` 
- Use `window_variant_matrix` instead of `variant_matrix[i:i+1, :]`
- Calculate `window_max_n` per window instead of global `max_n`
- Use `window_max_n` instead of `max_n` in calculations
- Return `(["windows"], ...)` instead of `(["variants"], ...)`

