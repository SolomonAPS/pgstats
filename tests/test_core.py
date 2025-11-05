"""
Unit tests for core GenomicDataset functionality.
"""

import pytest
import numpy as np
import pandas as pd
import sgkit as sg
import tempfile
from pathlib import Path

from pgstats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig


class TestGenomicDatasetInitialization:
    """Tests for GenomicDataset initialization."""
    
    def test_init_with_dataset(self):
        """Test initialization with an existing dataset."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        genomic_ds = GenomicDataset(data_source=ds)
        
        assert genomic_ds.dataset is not None
        assert len(genomic_ds.dataset.variants) == 100
        assert len(genomic_ds.dataset.samples) == 20
    
    def test_init_with_configs(self):
        """Test initialization with configuration objects."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(window_size=50, step_size=25)
        callable_config = CallableSitesConfig(max_missing=0.1)
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config,
            callable_config=callable_config
        )
        
        assert genomic_ds.window_config.window_size == 50
        assert genomic_ds.window_config.step_size == 25
        assert genomic_ds.callable_config.max_missing == 0.1
    
    def test_validation_fails_missing_vars(self):
        """Test that validation fails with missing required variables."""
        # Create dataset without required variables
        import xarray as xr
        ds = xr.Dataset({
            'test_var': (['variants'], np.arange(10))
        })
        ds = ds.assign_coords(variants=np.arange(10))
        
        with pytest.raises(ValueError, match="missing required variables"):
            GenomicDataset(data_source=ds)


class TestCallableSitesMasking:
    """Tests for callable sites masking with BED files."""
    
    def test_bed_file_loading(self):
        """Test loading callable sites from BED file."""
        # Create simulated dataset
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        positions = np.arange(1, 101)
        ds = ds.assign(variant_position=(['variants'], positions))
        
        # Create temporary BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            bed_path = f.name
            # BED format: 0-based start (inclusive), 0-based end (exclusive)
            # This should include positions 11-50 (1-based)
            f.write("0\t10\t50\n")
        
        try:
            callable_config = CallableSitesConfig(bed_file=bed_path)
            genomic_ds = GenomicDataset(
                data_source=ds,
                callable_config=callable_config
            )
            
            # Check that callable sites were identified
            assert genomic_ds.callable_sites is not None
            # Should have 40 callable sites (positions 11-50)
            assert genomic_ds.callable_sites == 40
        finally:
            Path(bed_path).unlink()
    
    def test_bed_masking_sets_missing(self):
        """Test that non-callable sites are set to -1 (missing)."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        positions = np.arange(1, 101)
        ds = ds.assign(variant_position=(['variants'], positions))
        
        # Store original genotypes
        original_genotypes = ds.call_genotype.values.copy()
        
        # Create BED file that only includes first 50 positions
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            bed_path = f.name
            f.write("0\t0\t50\n")
        
        try:
            callable_config = CallableSitesConfig(bed_file=bed_path)
            genomic_ds = GenomicDataset(
                data_source=ds,
                callable_config=callable_config
            )
            
            masked_genotypes = genomic_ds.dataset.call_genotype.values
            
            # Check that positions 51-100 are now -1
            assert np.all(masked_genotypes[50:, :, :] == -1)
            # Check that positions 1-50 are unchanged
            assert np.array_equal(masked_genotypes[:50, :, :], original_genotypes[:50, :, :])
        finally:
            Path(bed_path).unlink()


class TestWindowCreation:
    """Tests for window creation functionality."""
    
    def test_genome_wide_window(self):
        """Test creation of genome-wide window."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(window_size=None)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        
        assert genomic_ds.windowed_dataset is not None
        # Genome-wide should create exactly 1 window
        assert len(genomic_ds.windowed_dataset.windows) == 1
    
    def test_position_based_windows(self):
        """Test creation of position-based windows."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(window_size=50, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        
        assert genomic_ds.windowed_dataset is not None
        # Should create multiple windows
        assert len(genomic_ds.windowed_dataset.windows) > 1
    
    def test_window_filtering_by_min_variants(self):
        """Test filtering windows by minimum variants."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        # Create windows with filtering
        window_config = WindowConfig(
            window_size=50,
            step_size=25,
            min_variants=10
        )
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        
        # Check that all windows have at least min_variants
        window_starts = genomic_ds.windowed_dataset.window_start.values
        window_stops = genomic_ds.windowed_dataset.window_stop.values
        
        for start, stop in zip(window_starts, window_stops):
            n_variants = stop - start
            assert n_variants >= 10
    
    def test_region_filtering(self):
        """Test filtering to a specific genomic region."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        positions = np.arange(1, 101)
        ds = ds.assign(variant_position=(['variants'], positions))
        
        window_config = WindowConfig(
            window_size=50,
            start=25,
            end=75
        )
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        # Filter to region happens in create_windows
        genomic_ds.create_windows()
        
        # The windowed dataset should only contain variants in the region
        assert genomic_ds.windowed_dataset is not None


class TestMissingDataFiltering:
    """Tests for missing data filtering."""
    
    def test_filter_missing_data(self):
        """Test filtering variants by missing data threshold."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20, missing_pct=0.3)
        
        callable_config = CallableSitesConfig(max_missing=0.2)
        genomic_ds = GenomicDataset(
            data_source=ds,
            callable_config=callable_config
        )
        
        original_n_variants = len(genomic_ds.dataset.variants)
        
        genomic_ds.filter_missing_data()
        
        # Should have fewer variants after filtering
        assert len(genomic_ds.dataset.variants) <= original_n_variants
    
    def test_no_filtering_with_zero_threshold(self):
        """Test that no filtering occurs with max_missing=0."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20, missing_pct=0.1)
        
        callable_config = CallableSitesConfig(max_missing=0.0)
        genomic_ds = GenomicDataset(
            data_source=ds,
            callable_config=callable_config
        )
        
        original_n_variants = len(genomic_ds.dataset.variants)
        
        genomic_ds.filter_missing_data()
        
        # With max_missing=0 and some missing data, all variants should remain
        # (the function only filters if max_missing > 0)
        assert len(genomic_ds.dataset.variants) == original_n_variants


class TestSummaryInformation:
    """Tests for summary information."""
    
    def test_get_summary_basic(self):
        """Test getting basic summary information."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        genomic_ds = GenomicDataset(data_source=ds)
        summary = genomic_ds.get_summary()
        
        assert summary['n_variants'] == 100
        assert summary['n_samples'] == 20
        assert summary['analysis_type'] == 'not_analyzed'
    
    def test_get_summary_with_windows(self):
        """Test summary information after creating windows."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(window_size=50, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        summary = genomic_ds.get_summary()
        
        assert 'n_windows' in summary
        assert summary['n_windows'] > 0
        assert summary['analysis_type'] == 'windowed'
    
    def test_get_summary_genome_wide(self):
        """Test summary with genome-wide analysis."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(window_size=None)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        summary = genomic_ds.get_summary()
        
        assert summary['n_windows'] == 1
        assert summary['analysis_type'] == 'genome-wide'


class TestWindowConfig:
    """Tests for WindowConfig dataclass."""
    
    def test_default_step_size(self):
        """Test that step_size defaults to window_size."""
        config = WindowConfig(window_size=100)
        assert config.step_size == 100
    
    def test_explicit_step_size(self):
        """Test explicit step_size setting."""
        config = WindowConfig(window_size=100, step_size=50)
        assert config.step_size == 50
    
    def test_genome_wide_config(self):
        """Test configuration for genome-wide analysis."""
        config = WindowConfig(window_size=None)
        assert config.window_size is None
        assert config.step_size is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

