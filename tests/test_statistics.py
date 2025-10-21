"""
Unit tests for population genetics statistics.
"""

import pytest
import numpy as np
import sgkit as sg
import xarray as xr

from many_stats.stats.sfs_statistics import (
    # High-level statistics
    tajima_d,
    fu_li_d,
    fu_li_f,
    theta_pi,
    theta_w,
    theta_h,
    fay_wu_h,
    dh_joint_test,
    # Low-level functions for testing
    calculate_a1,
    calculate_a2,
    calculate_pi,
    calculate_S,
    calculate_theta_w,
    calculate_theta_h,
    get_unfolded_sfs,
    get_folded_sfs
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


class TestDHJointTest:
    """Tests for DH joint test combining Tajima's D and Fay and Wu's H."""
    
    def test_dh_joint_test_basic(self):
        """Test basic DH joint test calculation."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.0)
        
        result = dh_joint_test(ds)
        
        # Should contain both statistics
        assert 'tajima_d' in result.data_vars
        assert 'fay_wu_h' in result.data_vars
        assert len(result.tajima_d) == 50
        assert len(result.fay_wu_h) == 50
        
        # Values should be finite
        assert np.all(np.isfinite(result.tajima_d.values))
        assert np.all(np.isfinite(result.fay_wu_h.values))
    
    def test_dh_joint_test_consistency(self):
        """Test that DH joint test gives same results as individual tests."""
        ds = sg.simulate_genotype_call_dataset(n_variant=30, n_sample=15, missing_pct=0.1)
        
        # Calculate individually
        tajima_result = tajima_d(ds)
        fay_wu_result = fay_wu_h(ds)
        
        # Calculate jointly
        joint_result = dh_joint_test(ds)
        
        # Results should be identical
        np.testing.assert_array_almost_equal(
            joint_result.tajima_d.values, 
            tajima_result.tajima_d.values
        )
        np.testing.assert_array_almost_equal(
            joint_result.fay_wu_h.values, 
            fay_wu_result.fay_wu_h.values
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

