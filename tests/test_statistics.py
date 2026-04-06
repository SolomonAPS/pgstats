"""
Unit tests for population genetics statistics.
"""

import pytest
import numpy as np
import sgkit as sg
import xarray as xr

from pgstats.stats.sfs_statistics import (
    # High-level statistics
    tajima_d,
    fu_li_d,
    fu_li_f,
    theta_pi,
    theta_w,
    theta_h,
    fay_wu_h,
    # Low-level functions for testing
    calculate_a1,
    calculate_a2,
    calculate_pi,
    calculate_S,
    calculate_theta_w,
    calculate_theta_h,
    get_unfolded_sfs,
    get_folded_sfs,
    windowed_sfs,
)


class TestHelperFunctions:
    """Tests for helper functions used in statistics calculations."""
    
    def test_calculate_a1(self):
        """Test calculation of a1 (harmonic number)."""
        # For n=5, a1 = 1 + 1/2 + 1/3 + 1/4 = 2.083333...
        a1 = calculate_a1(5)
        assert np.isclose(a1, 2.083333, atol=0.0001)
        
        # For n=2, a1 = 1
        a1 = calculate_a1(2)
        assert np.isclose(a1, 1.0)
    
    def test_calculate_a2(self):
        """Test calculation of a2 (sum of 1/i^2)."""
        # For n=5, a2 = 1 + 1/4 + 1/9 + 1/16 = 1.4236...
        a2 = calculate_a2(5)
        assert np.isclose(a2, 1.4236, atol=0.001)
        
        # For n=2, a2 = 1
        a2 = calculate_a2(2)
        assert np.isclose(a2, 1.0)
    
    def test_get_unfolded_sfs(self):
        """Test calculation of unfolded site frequency spectrum."""
        # Create simple variant matrix
        # 3 samples, 4 variants
        variant_matrix = np.array([
            [0, 0, 1],  # 1 derived allele
            [0, 1, 1],  # 2 derived alleles
            [1, 1, 1],  # 3 derived alleles (monomorphic)
            [0, 0, 0],  # 0 derived alleles (monomorphic)
        ], dtype=np.int8)
        
        sfs, n = get_unfolded_sfs(variant_matrix)
        
        assert n == 3
        # SFS should be [1, 1] (one variant with 1 derived, one with 2 derived)
        # Variants with 0 or n derived alleles are excluded
        assert sfs[0] == 1  # 1 variant with 1 derived allele
        assert sfs[1] == 1  # 1 variant with 2 derived alleles
    
    def test_get_folded_sfs(self):
        """Test calculation of folded site frequency spectrum."""
        # Create variant matrix
        variant_matrix = np.array([
            [0, 0, 1, 1],  # 2 derived (minor allele count = 2)
            [1, 1, 1, 0],  # 3 derived (minor allele count = 1)
            [1, 0, 0, 0],  # 1 derived (minor allele count = 1)
        ], dtype=np.int8)
        
        folded_sfs, n = get_folded_sfs(variant_matrix)
        
        assert n == 4
        # Should have 2 variants with minor allele count 1
        assert folded_sfs[0] == 2
        # Should have 1 variant with minor allele count 2
        assert folded_sfs[1] == 1
    
    def test_calculate_S(self):
        """Test calculation of number of segregating sites."""
        sfs = np.array([5, 3, 2, 1])
        S = calculate_S(sfs)
        assert S == 11  # Sum of all counts
    
    def test_calculate_pi(self):
        """Test calculation of nucleotide diversity (pi)."""
        # Simple test case
        variant_matrix = np.array([
            [0, 0, 1],  # 1 derived allele
            [0, 1, 1],  # 2 derived alleles
        ], dtype=np.int8)
        
        pi = calculate_pi(variant_matrix)
        
        # Pi should be > 0 for polymorphic data
        assert pi > 0


class TestWindowedSfs:
    """Per-window site frequency spectrum (full bin vector)."""

    def test_windowed_sfs_unfolded_one_window(self):
        import xarray as xr

        call_genotype = np.array(
            [
                [[0, 0], [0, 0], [1, 1]],
                [[0, 0], [1, 1], [1, 1]],
            ],
            dtype=np.int8,
        )
        ds = xr.Dataset(
            {
                "call_genotype": (["variants", "samples", "ploidy"], call_genotype),
                "window_start_idx": ("windows", np.array([0], dtype=np.int64)),
                "window_stop_idx": ("windows", np.array([2], dtype=np.int64)),
            },
            coords={
                "variants": np.arange(2),
                "samples": np.arange(3),
                "ploidy": np.arange(2),
                "windows": np.array([0]),
            },
        )
        out = windowed_sfs(ds, folded=False)
        assert "sfs_unfolded" in out.data_vars
        assert out.sfs_unfolded.shape == (1, 2)
        np.testing.assert_array_equal(out.sfs_unfolded.values[0], [1, 1])
        np.testing.assert_array_equal(out.sfs_bin.values, [1, 2])

    def test_windowed_sfs_folded_one_window(self):
        import xarray as xr

        call_genotype = np.array(
            [
                [[0, 0], [0, 0], [1, 1], [1, 1]],
                [[1, 1], [1, 1], [1, 1], [0, 0]],
                [[1, 1], [0, 0], [0, 0], [0, 0]],
            ],
            dtype=np.int8,
        )
        ds = xr.Dataset(
            {
                "call_genotype": (["variants", "samples", "ploidy"], call_genotype),
                "window_start_idx": ("windows", np.array([0], dtype=np.int64)),
                "window_stop_idx": ("windows", np.array([3], dtype=np.int64)),
            },
            coords={
                "variants": np.arange(3),
                "samples": np.arange(4),
                "ploidy": np.arange(2),
                "windows": np.array([0]),
            },
        )
        out = windowed_sfs(ds, folded=True)
        assert out.sfs_folded.shape == (1, 2)
        np.testing.assert_array_equal(out.sfs_folded.values[0], [2, 1])


class TestTajimaDStatistic:
    """Tests for Tajima's D calculation."""
    
    def test_tajima_d_basic(self):
        """Test basic Tajima's D calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.0)
        
        result = tajima_d(ds)
        
        assert 'tajima_d' in result.data_vars
        assert len(result.tajima_d) == 50
        # Values should be finite
        assert np.all(np.isfinite(result.tajima_d.values))
    
    def test_tajima_d_with_missing(self):
        """Test Tajima's D with missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.2)
        
        result = tajima_d(ds)
        
        assert 'tajima_d' in result.data_vars
        # Should handle missing data gracefully
        assert np.all(np.isfinite(result.tajima_d.values) | np.isnan(result.tajima_d.values))
    
    def test_tajima_d_monomorphic(self):
        """Test Tajima's D with monomorphic sites."""
        # Create dataset with all same genotype
        genotypes = np.zeros((10, 20, 2), dtype=np.int8)
        
        import xarray as xr
        ds = xr.Dataset({
            'call_genotype': (['variants', 'samples', 'ploidy'], genotypes),
            'variant_position': (['variants'], np.arange(10)),
            'variant_contig': (['variants'], np.zeros(10, dtype=int)),
            'variant_allele': (['variants', 'alleles'], np.array([['A', 'T']] * 10))
        })
        ds = ds.assign_coords({
            'variants': np.arange(10),
            'samples': [f"sample_{i}" for i in range(20)],
            'ploidy': [0, 1],
            'alleles': [0, 1],
            'contigs': [0],
            'contig_id': ('contigs', ['chr1'])
        })
        
        result = tajima_d(ds)
        
        # Tajima's D should be 0 or NaN for monomorphic sites
        assert np.all((result.tajima_d.values == 0) | np.isnan(result.tajima_d.values))


class TestFuLiStatistics:
    """Tests for Fu and Li's D* and F* statistics."""
    
    def test_fu_li_d_basic(self):
        """Test basic Fu and Li's D* calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = fu_li_d(ds)
        
        assert 'fu_li_d' in result.data_vars
        assert len(result.fu_li_d) == 50
    
    def test_fu_li_f_basic(self):
        """Test basic Fu and Li's F* calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = fu_li_f(ds)
        
        assert 'fu_li_f' in result.data_vars
        assert len(result.fu_li_f) == 50
    
    def test_fu_li_with_missing(self):
        """Test Fu and Li statistics with missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.2)
        
        result_d = fu_li_d(ds)
        result_f = fu_li_f(ds)
        
        # Should handle missing data
        assert 'fu_li_d' in result_d.data_vars
        assert 'fu_li_f' in result_f.data_vars


class TestDiversityStatistics:
    """Tests for diversity measures."""
    
    def test_theta_pi(self):
        """Test theta_pi (nucleotide diversity) calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = theta_pi(ds)
        
        assert 'theta_pi' in result.data_vars
        assert len(result.theta_pi) == 50
        # Pi should be non-negative
        assert np.all(result.theta_pi.values >= 0)
    
    def test_theta_w(self):
        """Test theta_w (Watterson's estimator) calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = theta_w(ds)
        
        assert 'theta_w' in result.data_vars
        assert len(result.theta_w) == 50
        # Theta_w should be non-negative
        assert np.all(result.theta_w.values >= 0)
    
    def test_theta_h(self):
        """Test theta_h (Fay and Wu's estimator) calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = theta_h(ds)
        
        assert 'theta_h' in result.data_vars
        assert len(result.theta_h) == 50
        # Theta_h should be non-negative
        assert np.all(result.theta_h.values >= 0)
    
    def test_backward_compatibility_aliases(self):
        """Test that old names still work for backward compatibility."""
        ds = sg.simulate_genotype_call_dataset(n_variant=20, n_sample=10)
        
        # Test aliases
        result_pi_old = nucleotide_diversity(ds)
        result_w_old = watterson_theta(ds)
        result_h_old = fay_wu_theta(ds)
        
        result_pi_new = theta_pi(ds)
        result_w_new = theta_w(ds)
        result_h_new = theta_h(ds)
        
        # Output variable names should be new names
        assert 'theta_pi' in result_pi_old.data_vars
        assert 'theta_w' in result_w_old.data_vars
        assert 'theta_h' in result_h_old.data_vars
    
    def test_diversity_with_missing(self):
        """Test diversity statistics with missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.2)
        
        result_pi = theta_pi(ds)
        result_theta_w = theta_w(ds)
        result_theta_h = theta_h(ds)
        
        assert 'theta_pi' in result_pi.data_vars
        assert 'theta_w' in result_theta_w.data_vars
        assert 'theta_h' in result_theta_h.data_vars


class TestStatisticsConsistency:
    """Tests for consistency between related statistics."""
    
    def test_theta_calculations_consistent(self):
        """Test that different theta calculations are in reasonable ranges."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=30)
        
        result_pi = theta_pi(ds)
        result_theta_w = theta_w(ds)
        result_theta_h = theta_h(ds)
        
        pi_vals = result_pi.theta_pi.values
        theta_w_vals = result_theta_w.theta_w.values
        theta_h_vals = result_theta_h.theta_h.values
        
        # All should be non-negative
        assert np.all(pi_vals >= 0)
        assert np.all(theta_w_vals >= 0)
        assert np.all(theta_h_vals >= 0)
        
        # For neutral evolution, these should be in similar ranges
        # (though not necessarily equal)
        mean_pi = np.mean(pi_vals)
        mean_theta_w = np.mean(theta_w_vals)
        
        # They should be within an order of magnitude for simulated neutral data
        if mean_pi > 0 and mean_theta_w > 0:
            ratio = mean_pi / mean_theta_w
            assert 0.1 < ratio < 10.0


class TestStatisticsEdgeCases:
    """Tests for edge cases in statistics calculations."""
    
    def test_empty_dataset(self):
        """Test statistics with very small dataset."""
        # Minimum viable dataset
        ds = sg.simulate_genotype_call_dataset(n_variant=5, n_sample=5)
        
        # Should not raise errors
        result = tajima_d(ds)
        assert 'tajima_d' in result.data_vars
    
    def test_high_missing_data(self):
        """Test statistics with high proportion of missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.8)
        
        # Should handle gracefully
        result = tajima_d(ds)
        result_pi = theta_pi(ds)
        
        assert 'tajima_d' in result.data_vars
        assert 'theta_pi' in result_pi.data_vars
    
    def test_large_sample_size(self):
        """Test statistics with large sample size."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=100)
        
        result = tajima_d(ds)
        result_pi = theta_pi(ds)
        
        assert 'tajima_d' in result.data_vars
        assert 'theta_pi' in result_pi.data_vars


class TestNumbaFunctions:
    """Tests for numba-compiled helper functions."""
    
    def test_calculate_theta_w_direct(self):
        """Test direct calculation of Watterson's theta."""
        sfs = np.array([10, 5, 3, 2], dtype=np.int64)
        n = 5
        
        theta_w = calculate_theta_w(sfs, n)
        
        # Should be positive
        assert theta_w > 0
        # Should equal S/a1
        S = calculate_S(sfs)
        a1 = calculate_a1(n)
        expected = S / a1
        assert np.isclose(theta_w, expected)
    
    def test_calculate_theta_h_direct(self):
        """Test direct calculation of Fay and Wu's theta."""
        sfs = np.array([10, 5, 3, 2], dtype=np.int64)
        n = 5
        
        theta_h = calculate_theta_h(sfs, n)
        
        # Should be non-negative
        assert theta_h >= 0


class TestFayWuH:
    """Tests for Fay and Wu's H statistic."""
    
    def test_fay_wu_h_basic(self):
        """Test basic Fay and Wu's H calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.0)
        
        result = fay_wu_h(ds)
        
        assert 'fay_wu_h' in result.data_vars
        assert len(result.fay_wu_h) == 50
        # Values should be finite
        assert np.all(np.isfinite(result.fay_wu_h.values))
    
    def test_fay_wu_h_with_missing(self):
        """Test Fay and Wu's H with missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.2)
        
        result = fay_wu_h(ds)
        
        assert 'fay_wu_h' in result.data_vars
        # Should handle missing data gracefully
        assert np.all(np.isfinite(result.fay_wu_h.values) | np.isnan(result.fay_wu_h.values))
    
    def test_fay_wu_h_monomorphic(self):
        """Test Fay and Wu's H with monomorphic sites."""
        # Create dataset with all same genotype
        genotypes = np.zeros((10, 20, 2), dtype=np.int8)
        
        ds = xr.Dataset({
            'call_genotype': (['variants', 'samples', 'ploidy'], genotypes),
            'variant_position': (['variants'], np.arange(10)),
            'variant_contig': (['variants'], np.zeros(10, dtype=int)),
            'variant_allele': (['variants', 'alleles'], np.array([['A', 'T']] * 10))
        })
        
        result = fay_wu_h(ds)
        
        assert 'fay_wu_h' in result.data_vars
        # All values should be 0 for monomorphic sites
        assert np.all(result.fay_wu_h.values == 0)


class TestLdDecay:
    """Tests for LD decay (mean r-squared by distance bin)."""

    def _make_ld_dataset(self):
        """Two SNPs in perfect LD at distance 100 bp, one independent SNP at 600 bp."""
        genotypes = np.array([
            # SNP at pos 100: alleles perfectly correlated with SNP at pos 200
            [[0, 0], [0, 0], [1, 1], [1, 1]],
            # SNP at pos 200
            [[0, 0], [0, 0], [1, 1], [1, 1]],
            # SNP at pos 700: independent pattern
            [[1, 1], [0, 0], [1, 1], [0, 0]],
        ], dtype=np.int8)
        positions = np.array([100, 200, 700])
        ds = xr.Dataset({
            'call_genotype': (['variants', 'samples', 'ploidy'], genotypes),
            'variant_position': (['variants'], positions),
            'variant_contig': (['variants'], np.zeros(3, dtype=int)),
        })
        return ds

    def test_basic_output_shape(self):
        from pgstats.stats.ld_statistics import ld_decay
        ds = self._make_ld_dataset()
        result = ld_decay(ds, max_distance=1000, bin_size=500)
        assert 'mean_r_squared' in result.data_vars
        assert 'n_pairs' in result.data_vars
        assert 'bin_start' in result.data_vars
        assert 'bin_end' in result.data_vars
        assert 'bin_midpoint' in result.data_vars
        assert result.sizes['distance_bins'] == 2  # 1000/500 = 2 bins

    def test_bin_boundaries(self):
        from pgstats.stats.ld_statistics import ld_decay
        ds = self._make_ld_dataset()
        result = ld_decay(ds, max_distance=1000, bin_size=500)
        np.testing.assert_array_equal(result['bin_start'].values, [1, 501])
        np.testing.assert_array_equal(result['bin_end'].values, [500, 1000])
        np.testing.assert_array_equal(result['bin_midpoint'].values, [250.5, 750.5])

    def test_perfect_ld_pair(self):
        """SNPs at pos 100 and 200 (dist=100, bin 0) should have r2=1."""
        from pgstats.stats.ld_statistics import ld_decay
        ds = self._make_ld_dataset()
        result = ld_decay(ds, max_distance=1000, bin_size=500)
        # Bin 0 [1-500]: pairs (100,200)=dist 100, (200,700)=dist 500
        # Bin 1 [501-1000]: pair (100,700)=dist 600
        assert result['n_pairs'].values[0] == 2  # two pairs in bin 0
        assert result['n_pairs'].values[1] == 1  # one pair in bin 1

    def test_max_distance_filter(self):
        """Pairs beyond max_distance should be excluded."""
        from pgstats.stats.ld_statistics import ld_decay
        ds = self._make_ld_dataset()
        result = ld_decay(ds, max_distance=500, bin_size=500)
        assert result.sizes['distance_bins'] == 1
        # Only pair (100,200) dist=100 and (200,700) dist=500 are within range
        assert result['n_pairs'].values[0] == 2

    def test_maf_filter(self):
        """With high min_maf, monomorphic-ish variants get dropped."""
        from pgstats.stats.ld_statistics import ld_decay
        # All 3 SNPs have MAF=0.5 in our fixture, so min_maf=0.4 keeps all
        ds = self._make_ld_dataset()
        result_no_filter = ld_decay(ds, max_distance=1000, bin_size=500)
        result_with_filter = ld_decay(ds, max_distance=1000, bin_size=500, min_maf=0.4)
        np.testing.assert_array_equal(
            result_no_filter['n_pairs'].values,
            result_with_filter['n_pairs'].values,
        )

    def test_empty_bins_are_nan(self):
        """Bins with no pairs should have NaN mean_r_squared."""
        from pgstats.stats.ld_statistics import ld_decay
        ds = self._make_ld_dataset()
        # bin_size=50, max_distance=1000 → many empty bins
        result = ld_decay(ds, max_distance=1000, bin_size=50)
        empty_mask = result['n_pairs'].values == 0
        assert np.any(empty_mask)
        assert np.all(np.isnan(result['mean_r_squared'].values[empty_mask]))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

