from setuptools import setup, find_packages

setup(
    name="pgstats",
    version="0.1.0",
    description="Population Genetics Statistics Toolkit - A comprehensive statistical genomics toolkit built on sgkit",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    python_requires=">=3.8,<3.13",
    install_requires=[
        "sgkit>=0.10.0",
        "numba>=0.56.0",
        "numpy>=1.20.0,<3.0.0",
        "dask>=2021.0.0",
        "xarray>=2023.8.0",  # Required for sgkit 0.10.0 (needs xarray.namedarray module)
        "zarr>=2.10.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "bio2zarr>=0.1.6",
        "pyranges>=0.0.129",
    ],
    entry_points={
        'console_scripts': [
            'pgstats=pgstats.cli.main:main',
        ],
    },
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
            "nbsphinx",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
