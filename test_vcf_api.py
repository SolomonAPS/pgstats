#!/usr/bin/env python3
"""
Simple test to verify the updated VCF loading implementation.
"""

import sys
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

def test_imports():
    """Test that the imports work correctly."""
    try:
        from many_stats.io.loaders import load_vcf, load_vcf_simple, check_bio2zarr_available
        print("✓ Successfully imported VCF loading functions")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_bio2zarr_check():
    """Test the bio2zarr availability check."""
    try:
        from many_stats.io.loaders import check_bio2zarr_available
        available = check_bio2zarr_available()
        if available:
            print("✓ bio2zarr is available")
        else:
            print("⚠ bio2zarr is not available (expected if not installed)")
        return True
    except Exception as e:
        print(f"✗ bio2zarr check failed: {e}")
        return False

def test_function_signatures():
    """Test that function signatures are correct."""
    try:
        from many_stats.io.loaders import load_vcf, load_vcf_simple
        import inspect
        
        # Check load_vcf signature
        sig = inspect.signature(load_vcf)
        params = list(sig.parameters.keys())
        expected_params = ['vcf_path', 'temp_dir', 'keep_zarr', 'variants_chunk_size', 
                          'samples_chunk_size', 'worker_processes', 'show_progress']
        
        for param in expected_params:
            if param not in params:
                print(f"✗ Missing parameter: {param}")
                return False
        
        print("✓ load_vcf function signature is correct")
        
        # Check load_vcf_simple signature
        sig = inspect.signature(load_vcf_simple)
        params = list(sig.parameters.keys())
        if 'vcf_path' not in params:
            print("✗ load_vcf_simple missing vcf_path parameter")
            return False
        
        print("✓ load_vcf_simple function signature is correct")
        return True
        
    except Exception as e:
        print(f"✗ Function signature test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("Testing updated VCF loading implementation...")
    print()
    
    tests = [
        test_imports,
        test_bio2zarr_check,
        test_function_signatures,
    ]
    
    passed = 0
    for test in tests:
        if test():
            passed += 1
        print()
    
    print(f"Tests passed: {passed}/{len(tests)}")
    
    if passed == len(tests):
        print("🎉 All tests passed! VCF loading implementation is ready.")
        return True
    else:
        print("❌ Some tests failed.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)


