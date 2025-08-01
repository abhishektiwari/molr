"""Setup script for MolR - Molecular Representation package."""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read version from package
version_file = os.path.join(os.path.dirname(__file__), "molr", "__init__.py")
with open(version_file, "r") as f:
    for line in f:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip().strip('"').strip("'")
            break
    else:
        version = "0.1.0"

setup(
    name="molr",
    version=version,
    author="Abhishek Tiwari",
    author_email="",
    description="Molecular Realm for Spatial Indexed Structures - Fast spatial operations for molecular data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/abhishektiwari/molr",
    packages=find_packages(exclude=["tests", "tests.*", "experiments", "docs"]),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",  # For spatial indexing with KDTree
        "pyparsing>=3.0.0",  # For selection language parser
        "pdbreader>=0.1.0",  # For PDB file parsing
        "mmcif>=0.1.0",  # For mmCIF file parsing
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "mypy>=1.0.0",
            "flake8>=6.0.0",
            "isort>=5.0.0",
        ],
        "docs": [
            "sphinx>=5.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "sphinx-autodoc-typehints>=1.0.0",
        ],
    },
    package_data={
        "molr": [
            "constants/*.py",
            "py.typed",  # PEP 561 marker for type hints
        ],
    },
    include_package_data=True,
    zip_safe=False,
    entry_points={
        "console_scripts": [
            # Add any command-line scripts here if needed
        ],
    },
)
