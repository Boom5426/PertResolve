from setuptools import setup, find_packages

setup(
    name="alleleperturb",
    version="0.1.0",
    description="AllelePerturb: A Benchmark for Allele-Resolution Single-Cell Perturbation Prediction",
    author="Bob Zhang, Qianqian Song",
    author_email="boom5426@ufl.edu",
    url="https://github.com/Boom5426/AllelePerturb",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.23",
        "pandas>=1.5",
        "scipy>=1.10",
        "scikit-learn>=1.2",
    ],
    extras_require={
        "full": [
            "torch>=2.0",
            "scanpy>=1.9",
            "anndata>=0.9",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
