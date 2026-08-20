"""Fallback setup.py for pip < 21 / Python 3.6."""
from setuptools import setup, find_packages

setup(
    name="coveragemaster",
    version="1.1.0",
    description="CNV caller from samtools depth data using Wavelet-HMM at nucleotide resolution",
    long_description="",
    packages=find_packages(),
    package_data={"coveragemaster": ["REF/*"]},
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.16",
        "scipy>=1.2",
        "matplotlib>=2.2",
        "sympy>=1.0",
        "PyWavelets>=1.0",
        "pandas>=1.0",
        "more-itertools>=8.0",
    ],
    entry_points={
        "console_scripts": [
            "coveragemaster=coveragemaster.cli:main",
            "coveragemaster-build-ref=coveragemaster.tools.create_ref:main",
            "coveragemaster-compute-stats=coveragemaster.tools.compute_stats:main",
        ]
    },
)
