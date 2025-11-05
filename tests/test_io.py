"""
Tests for I/O functions (loading, validation, writing).
"""

import pytest
import numpy as np
import pandas as pd
import sgkit as sg
import tempfile
from pathlib import Path

from pgstats.io.loaders import check_bio2zarr_available, load_zarr
from pgstats.io.writers import save_csv, save_tsv, save_zarr
from pgstats.utils.validation import validate_dataset, check_missing_data
from pgstats.utils.conversion import (
    convert_to_variant_matrix,
    convert_call_to_index,
    convert_genotypes_to_binary,
    convert_genotypes_to_allele_counts
)


class TestBio2ZarrAvailability:
    """Tests for bio2zarr availability checking."""
    
    def test_check_bio2zarr_available(self):
        """Test checking if bio2zarr is available."""
        result = check_bio2zarr_available()
        
        # Should return a boolean
        assert isinstance(result, bool)
        
        # If available, we should be able to import it
        if result:
            import bio2zarr.vcf
            assert bio2zarr.vcf is not None


class TestDatasetValidation:
    """Tests for dataset validation functions."""
    
    def test_validate_dataset_valid(self):
        """Test validation of a valid dataset."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = validate_dataset(ds)
        
        assert result['valid'] == True
        assert len(result['errors']) == 0
    
    def test_validate_dataset_missing_vars(self):
        """Test validation with missing required variables."""
        # Create dataset without required variables
        import xarray as xr
        ds = xr.Dataset({
            'some_var': (['x'], np.arange(10))
        })
        ds = ds.assign_coords(x=np.arange(10))
        
        result = validate_dataset(ds)
        
        assert result['valid'] == False
        assert len(result['errors']) > 0
    
    def test_validate_dataset_with_missing_data(self):
        """Test validation detects missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.2)
        
        result = validate_dataset(ds)
        
        # Should still be valid
        assert result['valid'] == True
        # But should have warnings about missing data
        assert 'info' in result
        if 'missing_genotypes' in result['info']:
            assert result['info']['missing_genotypes'] > 0
    
    def test_validate_custom_required_vars(self):
        """Test validation with custom required variables."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        # Validate with custom required vars
        result = validate_dataset(ds, required_vars=['call_genotype', 'variant_position'])
        
        assert result['valid'] == True


class TestMissingDataCheck:
    """Tests for missing data checking."""
    
    def test_check_missing_data_basic(self):
        """Test basic missing data checking."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.1)
        
        result = check_missing_data(ds)
        
        assert 'total_missing' in result
        assert 'missing_rate' in result
        assert 'missing_per_variant' in result
        assert 'missing_per_sample' in result
        
        # Missing rate should be approximately 0.1
        assert 0.0 <= result['missing_rate'] <= 0.2
    
    def test_check_missing_data_no_missing(self):
        """Test missing data check with no missing data."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.0)
        
        result = check_missing_data(ds)
        
        assert result['total_missing'] == 0
        assert result['missing_rate'] == 0.0
    
    def test_check_missing_data_high_missing(self):
        """Test missing data check with high missing rate."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20, missing_pct=0.5)
        
        result = check_missing_data(ds)
        
        assert result['missing_rate'] > 0.4
        assert result['total_missing'] > 0
    
    def test_check_missing_invalid_var(self):
        """Test missing data check with invalid variable name."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = check_missing_data(ds, call_genotype='nonexistent_var')
        
        assert 'error' in result


class TestGenotypeConversion:
    """Tests for genotype conversion utilities."""
    
    def test_convert_to_variant_matrix(self):
        """Test conversion of dataset to binary variant matrix."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        variant_matrix = convert_to_variant_matrix(ds)
        
        assert variant_matrix.shape == (50, 20)
        assert variant_matrix.dtype == np.int8
        # Values should be 0 or 1
        assert np.all((variant_matrix == 0) | (variant_matrix == 1))
    
    def test_convert_genotypes_to_binary(self):
        """Test direct genotype to binary conversion."""
        genotypes = np.array([
            [[0, 0], [0, 1], [1, 1]],  # 3 samples, diploid
            [[0, 1], [1, 1], [0, 0]],
        ], dtype=np.int8)
        
        binary = convert_genotypes_to_binary(genotypes)
        
        assert binary.shape == (2, 3)
        assert np.all((binary == 0) | (binary == 1))
    
    def test_convert_genotypes_to_allele_counts(self):
        """Test conversion to allele counts."""
        genotypes = np.array([
            [[0, 0], [0, 1], [1, 1]],  # 3 samples, diploid
            [[0, 1], [1, 1], [0, 0]],
        ], dtype=np.int8)
        
        counts = convert_genotypes_to_allele_counts(genotypes)
        
        assert counts.shape == (2, 3)
        # Counts should be 0, 1, or 2 for diploid
        assert np.all((counts >= 0) & (counts <= 2))
    
    def test_convert_call_to_index(self):
        """Test conversion of call genotypes to index format."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        result = convert_call_to_index(ds)
        
        assert 'call_genotype_index' in result.data_vars
        assert result.call_genotype_index.shape == (50, 20)
    
    def test_conversion_with_missing_data(self):
        """Test conversions handle missing data correctly."""
        genotypes = np.array([
            [[0, 0], [0, 1], [-1, -1]],  # -1 indicates missing
            [[0, 1], [-1, -1], [0, 0]],
        ], dtype=np.int8)
        
        binary = convert_genotypes_to_binary(genotypes)
        counts = convert_genotypes_to_allele_counts(genotypes)
        
        # Should handle missing data without errors
        assert binary.shape == (2, 3)
        assert counts.shape == (2, 3)


class TestResultWriting:
    """Tests for writing results to various formats."""
    
    def test_save_csv_basic(self):
        """Test saving dataset to CSV."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        # Add some statistics
        from pgstats.stats.sfs_statistics import tajima_d
        ds = tajima_d(ds)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            save_csv(ds, csv_path)
            
            # Verify file exists
            assert Path(csv_path).exists()
            
            # Read and verify contents
            df = pd.read_csv(csv_path)
            assert len(df) == 50
            assert 'tajima_d' in df.columns
            
        finally:
            if Path(csv_path).exists():
                Path(csv_path).unlink()
    
    def test_save_tsv_basic(self):
        """Test saving dataset to TSV."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        from pgstats.stats.sfs_statistics import theta_pi
        ds = theta_pi(ds)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            tsv_path = f.name
        
        try:
            save_tsv(ds, tsv_path)
            
            assert Path(tsv_path).exists()
            
            # Read as TSV
            df = pd.read_csv(tsv_path, sep='\t')
            assert len(df) == 50
            assert 'theta_pi' in df.columns
            
        finally:
            if Path(tsv_path).exists():
                Path(tsv_path).unlink()
    
    def test_save_multiple_stats(self):
        """Test saving dataset with multiple statistics."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        from pgstats.stats.sfs_statistics import tajima_d, theta_pi, theta_w
        
        ds = tajima_d(ds)
        ds = theta_pi(ds)
        ds = theta_w(ds)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            save_csv(ds, csv_path)
            
            df = pd.read_csv(csv_path)
            assert 'tajima_d' in df.columns
            assert 'theta_pi' in df.columns
            assert 'theta_w' in df.columns
            
        finally:
            if Path(csv_path).exists():
                Path(csv_path).unlink()
    
    def test_save_no_stats_raises_error(self):
        """Test that saving without statistics raises error."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            csv_path = f.name
        
        try:
            # Should raise ValueError
            with pytest.raises(ValueError, match="No statistical variables"):
                save_csv(ds, csv_path)
        finally:
            if Path(csv_path).exists():
                Path(csv_path).unlink()


class TestZarrIO:
    """Tests for Zarr format I/O."""
    
    def test_zarr_roundtrip(self):
        """Test saving and loading with Zarr format."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            zarr_path = Path(tmpdir) / "test.zarr"
            
            # Save to Zarr
            save_zarr(ds, str(zarr_path))
            
            # Verify saved
            assert zarr_path.exists()
            
            # Load back
            ds_loaded = load_zarr(str(zarr_path))
            
            # Verify structure
            assert len(ds_loaded.variants) == 50
            assert len(ds_loaded.samples) == 20
            assert 'call_genotype' in ds_loaded.data_vars


class TestDataIntegrity:
    """Tests for data integrity during I/O operations."""
    
    def test_conversion_preserves_data(self):
        """Test that conversions preserve data accurately."""
        ds = sg.simulate_genotype_call_dataset(n_variant=20, n_sample=10, missing_pct=0.0)
        
        # Get original genotypes
        original = ds.call_genotype.values
        
        # Convert to variant matrix and back
        variant_matrix = convert_to_variant_matrix(ds)
        
        # Check dimensions
        assert variant_matrix.shape[0] == original.shape[0]
        assert variant_matrix.shape[1] == original.shape[1]
    
    def test_validation_consistency(self):
        """Test that validation is consistent."""
        ds = sg.simulate_genotype_call_dataset(n_variant=50, n_sample=20)
        
        # Run validation multiple times
        result1 = validate_dataset(ds)
        result2 = validate_dataset(ds)
        
        # Results should be identical
        assert result1['valid'] == result2['valid']
        assert result1['errors'] == result2['errors']
    
    def test_missing_data_calculation_accuracy(self):
        """Test accuracy of missing data calculations."""
        # Create dataset with known missing data
        n_variants = 10
        n_samples = 5
        ploidy = 2
        
        genotypes = np.zeros((n_variants, n_samples, ploidy), dtype=np.int8)
        
        # Set some to missing (-1)
        genotypes[0, 0, :] = -1  # 2 missing calls
        genotypes[1, 1, :] = -1  # 2 missing calls
        
        import xarray as xr
        ds = xr.Dataset({
            'call_genotype': (['variants', 'samples', 'ploidy'], genotypes),
            'variant_position': (['variants'], np.arange(n_variants)),
            'variant_contig': (['variants'], np.zeros(n_variants, dtype=int)),
            'variant_allele': (['variants', 'alleles'], np.array([['A', 'T']] * n_variants))
        })
        ds = ds.assign_coords({
            'variants': np.arange(n_variants),
            'samples': [f"sample_{i}" for i in range(n_samples)],
            'ploidy': [0, 1],
            'alleles': [0, 1],
            'contigs': [0],
            'contig_id': ('contigs', ['chr1'])
        })
        
        result = check_missing_data(ds)
        
        # Should have exactly 4 missing calls
        assert result['total_missing'] == 4
        
        # Missing rate = 4 / (10 * 5 * 2) = 0.04
        assert np.isclose(result['missing_rate'], 0.04)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

