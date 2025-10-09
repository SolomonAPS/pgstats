"""
Basic tests for many-stats functionality.
"""

import pytest
import numpy as np
import sgkit as sg
import many_stats as ms


def test_basic_functionality():
    """Test basic functionality of many-stats."""
    
    # Create a simple test dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20, missing_pct=0.1)
    
    # Test validation
    validation_results = ms.utils.validate_dataset(ds)
    assert validation_results["valid"] == True
    
    # Test missing data check
    missing_stats = ms.utils.check_missing_data(ds)
    assert "total_missing" in missing_stats
    assert "missing_rate" in missing_stats
    
    # Test data conversion
    binary_matrix = ms.utils.convert_to_variant_matrix(ds)
    assert binary_matrix.shape == (100, 20)
    
    # Test genotype index conversion
    ds_with_index = ms.utils.convert_call_to_index(ds)
    assert "call_genotype_index" in ds_with_index.data_vars


def test_population_statistics():
    """Test population statistics calculations."""
    
    # Create test dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.1)
    
    # Test Tajima's D
    ds_with_tajima = ms.stats.tajima_d(ds)
    assert "tajima_d" in ds_with_tajima.data_vars
    assert len(ds_with_tajima.tajima_d.values) == 50
    
    # Test Fu and Li's statistics
    ds_with_fu_li_d = ms.stats.fu_li_d(ds)
    assert "fu_li_d" in ds_with_fu_li_d.data_vars
    
    ds_with_fu_li_f = ms.stats.fu_li_f(ds)
    assert "fu_li_f" in ds_with_fu_li_f.data_vars


def test_diversity_statistics():
    """Test diversity statistics calculations."""
    
    # Create test dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.1)
    
    # Test nucleotide diversity
    ds_with_pi = ms.stats.nucleotide_diversity(ds)
    assert "nucleotide_diversity" in ds_with_pi.data_vars
    
    # Test Watterson's theta
    ds_with_theta_w = ms.stats.watterson_theta(ds)
    assert "watterson_theta" in ds_with_theta_w.data_vars
    
    # Test Fay and Wu's theta
    ds_with_theta_h = ms.stats.fay_wu_theta(ds)
    assert "fay_wu_theta" in ds_with_theta_h.data_vars


def test_selection_statistics():
    """Test selection statistics calculations."""
    
    # Create test dataset
    ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.1)
    
    # Test HKA analysis
    ds_with_hka = ms.stats.hudson_kreitman_aguade_analysis(ds)
    assert "hka_statistic" in ds_with_hka.data_vars
    assert "hka_p_value" in ds_with_hka.data_vars


if __name__ == "__main__":
    # Run basic tests
    test_basic_functionality()
    test_population_statistics()
    test_diversity_statistics()
    test_selection_statistics()
    print("All tests passed!")
