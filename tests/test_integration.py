"""
Integration tests for end-to-end workflows.
"""

import pytest
import numpy as np
import pandas as pd
import sgkit as sg
import tempfile
from pathlib import Path

from pgstats.core.dataset import GenomicDataset, WindowConfig, CallableSitesConfig


class TestEndToEndWorkflows:
    """Tests for complete end-to-end workflows."""
    
    def test_simulated_data_windowed_analysis(self):
        """Test complete workflow with simulated data and windowed analysis."""
        # Create simulated dataset
        ds = sg.simulate_genotype_call_dataset(n_variant=200, n_sample=30, missing_pct=0.1)
        
        # Configure analysis
        window_config = WindowConfig(
            window_size=50,
            step_size=25,
            min_variants=5
        )
        
        callable_config = CallableSitesConfig(max_missing=0.2)
        
        # Initialize dataset
        genomic_ds = GenomicDataset(
            data_source=ds,
            callable_config=callable_config,
            window_config=window_config
        )
        
        # Filter missing data
        genomic_ds.filter_missing_data()
        
        # Create windows
        genomic_ds.create_windows()
        assert genomic_ds.windowed_dataset is not None
        
        # Calculate windowed statistics
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta']
        )
        
        # Verify results
        assert 'tajima_d' in window_stats.data_vars
        assert 'nucleotide_diversity' in window_stats.data_vars
        assert 'watterson_theta' in window_stats.data_vars
        assert len(window_stats.windows) > 0
    
    def test_genome_wide_analysis(self):
        """Test complete workflow with genome-wide analysis."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        # Genome-wide configuration
        window_config = WindowConfig(window_size=None)
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        # Create genome-wide window
        genomic_ds.create_windows()
        
        # Calculate genome-wide statistics
        genome_stats = genomic_ds.calculate_genome_wide_stats(
            stats=['tajima_d', 'nucleotide_diversity', 'watterson_theta']
        )
        
        # Verify results
        assert 'tajima_d' in genome_stats
        assert 'nucleotide_diversity' in genome_stats
        assert 'watterson_theta' in genome_stats
        
        # Values should be finite
        assert np.isfinite(genome_stats['tajima_d'])
        assert np.isfinite(genome_stats['nucleotide_diversity'])
        assert np.isfinite(genome_stats['watterson_theta'])
    
    def test_workflow_with_bed_masking(self):
        """Test complete workflow with BED file masking."""
        # Create simulated dataset
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        positions = np.arange(1, 101)
        ds = ds.assign(variant_position=(['variants'], positions))
        
        # Create temporary BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            bed_path = f.name
            # Include positions 1-50
            f.write("0\t0\t50\n")
        
        try:
            # Configure with BED masking
            callable_config = CallableSitesConfig(
                bed_file=bed_path,
                max_missing=0.1
            )
            
            window_config = WindowConfig(
                window_size=25,
                step_size=25,
                min_variants=3
            )
            
            # Run analysis
            genomic_ds = GenomicDataset(
                data_source=ds,
                callable_config=callable_config,
                window_config=window_config
            )
            
            # Verify callable sites were masked
            assert genomic_ds.callable_sites == 50
            
            # Continue with analysis
            genomic_ds.create_windows()
            window_stats = genomic_ds.calculate_windowed_stats(
                stats=['tajima_d', 'nucleotide_diversity']
            )
            
            assert 'tajima_d' in window_stats.data_vars
            assert 'nucleotide_diversity' in window_stats.data_vars
            
        finally:
            Path(bed_path).unlink()
    
    def test_region_specific_analysis(self):
        """Test workflow with region-specific analysis."""
        ds = sg.simulate_genotype_call_dataset(n_variant=200, n_sample=20)
        positions = np.arange(1, 201)
        ds = ds.assign(variant_position=(['variants'], positions))
        
        # Configure for specific region
        window_config = WindowConfig(
            window_size=25,
            start=50,
            end=150,
            min_variants=3
        )
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        # Create windows for region
        genomic_ds.create_windows()
        
        # Calculate statistics
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'nucleotide_diversity']
        )
        
        # Verify analysis was restricted to region
        assert window_stats is not None
        assert len(window_stats.windows) > 0
    
    def test_multiple_statistics(self):
        """Test calculation of multiple statistics at once."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=25)
        
        window_config = WindowConfig(window_size=50, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        
        # Calculate all available statistics
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'fu_li_d', 'fu_li_f', 
                   'nucleotide_diversity', 'watterson_theta', 'fay_wu_theta']
        )
        
        # Verify all statistics are present
        assert 'tajima_d' in window_stats.data_vars
        assert 'fu_li_d' in window_stats.data_vars
        assert 'fu_li_f' in window_stats.data_vars
        assert 'nucleotide_diversity' in window_stats.data_vars
        assert 'watterson_theta' in window_stats.data_vars
        assert 'fay_wu_theta' in window_stats.data_vars


class TestResultSaving:
    """Tests for saving analysis results."""
    
    def test_save_results_csv(self):
        """Test saving results to CSV format."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        window_config = WindowConfig(window_size=25, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        genomic_ds.calculate_windowed_stats(stats=['tajima_d', 'nucleotide_diversity'])
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            genomic_ds.save_results(csv_path, format='csv')
            
            # Verify file was created
            assert Path(csv_path).exists()
            
            # Verify contents can be read
            df = pd.read_csv(csv_path)
            assert len(df) > 0
            
        finally:
            if Path(csv_path).exists():
                Path(csv_path).unlink()
    
    def test_save_results_tsv(self):
        """Test saving results to TSV format."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        window_config = WindowConfig(window_size=25, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        genomic_ds.calculate_windowed_stats(stats=['tajima_d'])
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            tsv_path = f.name
        
        try:
            genomic_ds.save_results(tsv_path, format='tsv')
            
            assert Path(tsv_path).exists()
            
            # Verify TSV format
            df = pd.read_csv(tsv_path, sep='\t')
            assert len(df) > 0
            
        finally:
            if Path(tsv_path).exists():
                Path(tsv_path).unlink()


class TestDataQualityChecks:
    """Tests for data quality and validation."""
    
    def test_summary_information(self):
        """Test that summary information is accurate."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=25)
        
        window_config = WindowConfig(window_size=50, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        summary = genomic_ds.get_summary()
        
        # Verify basic counts
        assert summary['n_variants'] == 100
        assert summary['n_samples'] == 25
        
        # Create windows and check updated summary
        genomic_ds.create_windows()
        summary = genomic_ds.get_summary()
        
        assert 'n_windows' in summary
        assert summary['n_windows'] > 0
    
    def test_missing_data_handling(self):
        """Test that missing data is handled correctly throughout workflow."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20, missing_pct=0.3)
        
        callable_config = CallableSitesConfig(max_missing=0.2)
        window_config = WindowConfig(window_size=50, step_size=25)
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            callable_config=callable_config,
            window_config=window_config
        )
        
        original_n = len(genomic_ds.dataset.variants)
        
        # Filter by missing data
        genomic_ds.filter_missing_data()
        
        filtered_n = len(genomic_ds.dataset.variants)
        
        # Should have fewer variants after filtering
        assert filtered_n <= original_n
        
        # Statistics should still work
        genomic_ds.create_windows()
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'nucleotide_diversity']
        )
        
        assert window_stats is not None


class TestComplexScenarios:
    """Tests for complex realistic scenarios."""
    
    def test_multi_contig_dataset(self):
        """Test analysis with multiple contigs."""
        ds = sg.simulate_genotype_call_dataset(
            n_variant=150,
            n_sample=20,
            n_contig=3
        )
        
        window_config = WindowConfig(window_size=50, step_size=25)
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d', 'nucleotide_diversity']
        )
        
        # Should handle multiple contigs
        assert window_stats is not None
        assert len(window_stats.windows) > 0
    
    def test_large_window_small_variants(self):
        """Test with window size larger than variant count."""
        ds = sg.simulate_genotype_call_dataset(n_variant=20, n_sample=15)
        
        window_config = WindowConfig(
            window_size=1000,  # Large window
            min_variants=5
        )
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        
        # Should create at least one window
        assert genomic_ds.windowed_dataset is not None
    
    def test_overlapping_windows(self):
        """Test with overlapping windows (step < window_size)."""
        ds = sg.simulate_genotype_call_dataset(n_variant=100, n_sample=20)
        
        window_config = WindowConfig(
            window_size=50,
            step_size=10,  # Overlapping windows
            min_variants=5
        )
        
        genomic_ds = GenomicDataset(
            data_source=ds,
            window_config=window_config
        )
        
        genomic_ds.create_windows()
        window_stats = genomic_ds.calculate_windowed_stats(
            stats=['tajima_d']
        )
        
        # Should handle overlapping windows
        assert len(window_stats.windows) > 0
        # More windows due to smaller step size
        assert len(window_stats.windows) > 5


class TestErrorHandling:
    """Tests for error handling and edge cases."""
    
    def test_calculate_stats_before_windows(self):
        """Test that calculating stats before creating windows raises error."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        genomic_ds = GenomicDataset(data_source=ds)
        
        # Should raise ValueError
        with pytest.raises(ValueError, match="Windows not created"):
            genomic_ds.calculate_windowed_stats(stats=['tajima_d'])
    
    def test_save_before_calculating_stats(self):
        """Test that saving before calculating stats raises error."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        genomic_ds = GenomicDataset(data_source=ds)
        
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            # Should raise ValueError
            with pytest.raises(ValueError, match="No window statistics"):
                genomic_ds.save_results(csv_path)
        finally:
            if Path(csv_path).exists():
                Path(csv_path).unlink()
    
    def test_invalid_bed_file(self):
        """Test error handling with non-existent BED file."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        callable_config = CallableSitesConfig(
            bed_file="/nonexistent/path/to/file.bed"
        )
        
        # Should raise FileNotFoundError
        with pytest.raises(FileNotFoundError):
            GenomicDataset(
                data_source=ds,
                callable_config=callable_config
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])



class TestWindowsWithoutComputableLD:
    """A window where every SNP pair has NaN LD must stay in the output with NaN LD columns."""

    def _dataset(self):
        ds = sg.simulate_genotype_call_dataset(n_variant=20, n_sample=6, missing_pct=0.0, seed=0)
        ds["variant_position"] = ("variants", np.arange(1, 21, dtype=np.int64))
        gt = ds.call_genotype.values.copy()
        # Second window (variant indices ~10-19):
        #   variant 9/10 called only in samples 0-2, variant 10 only in samples 3-5
        #   -> no jointly called gamete; variants 12-19 monomorphic among the samples
        #   where they are called -> every pair has undefined r2 / D'.
        # Index 9 is included too so the construction holds whichever side of
        # the boundary the windowing code puts position 10 on.
        gt[9:20] = -1
        gt[9, 0:3] = [[0, 1], [1, 0], [0, 0]]
        gt[10, 3:6] = [[0, 1], [1, 1], [0, 0]]
        gt[11:20, 0:3] = 0
        ds["call_genotype"] = (("variants", "samples", "ploidy"), gt)
        ds["call_genotype_mask"] = (("variants", "samples", "ploidy"), gt < 0)
        return ds

    def test_rows_kept_with_nan_ld(self):
        genomic_ds = GenomicDataset(
            data_source=self._dataset(),
            callable_config=CallableSitesConfig(),
            window_config=WindowConfig(window_size=10, step_size=10, min_variants=5),
        )
        genomic_ds.create_windows()
        n_windows = len(genomic_ds.windowed_dataset.windows)
        result = genomic_ds.calculate_windowed_stats(stats=["ld_r2", "tajima_d"])
        assert len(result.windows) == n_windows == 2
        r2 = result["mean_r_squared"].values.astype(float)
        assert np.isnan(r2[1]) and not np.isnan(r2[0])
        assert np.isfinite(result["tajima_d"].values.astype(float)).all()
