![Molr](molr.png)

![GitHub Release](https://img.shields.io/github/v/release/abhishektiwari/molr)
![GitHub Actions Test Workflow Status](https://img.shields.io/github/actions/workflow/status/abhishektiwari/molr/test.yml?label=tests)
![PyPI - Version](https://img.shields.io/pypi/v/molr)
![Python Wheels](https://img.shields.io/pypi/wheel/molr)
![Python Versions](https://img.shields.io/pypi/pyversions/molr?logo=python&logoColor=white)
![GitHub last commit](https://img.shields.io/github/last-commit/abhishektiwari/molr)
![PyPI - Status](https://img.shields.io/pypi/status/molr)
![Conda Version](https://img.shields.io/conda/v/molr/molr)
![License](https://img.shields.io/github/license/abhishektiwari/molr)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/abhishektiwari/molr/total?label=GitHub%20Downloads)
![PyPI Downloads](https://img.shields.io/pepy/dt/molr?label=PyPI%20Downloads)
[![codecov](https://codecov.io/gh/abhishektiwari/molr/graph/badge.svg?token=QSKYLB3M1V)](https://codecov.io/gh/abhishektiwari/molr)
[![Socket](https://socket.dev/api/badge/pypi/package/molr/0.0.2?artifact_id=py3-none-any-whl)](https://socket.dev/pypi/package/molr/overview/0.0.2/py3-none-any-whl)
[![CodeFactor](https://www.codefactor.io/repository/github/abhishektiwari/molr/badge/main)](https://www.codefactor.io/repository/github/abhishektiwari/molr/overview/main)

# MolR - Molecular Realm for Spatial Indexed Structures

A high-performance Python package that creates a spatial realm for molecular structures, providing lightning-fast neighbor searches, geometric queries, and spatial operations through integrated KDTree indexing.


## Features

### High-Performance Structure Representation
- NumPy-based Structure class with Structure of Arrays (SoA)
- Efficient spatial indexing with scipy KDTree integration for O(log n) neighbor queries
- Memory-efficient trajectory handling with StructureEnsemble
- Lazy initialization of optional annotations to minimize memory usage

### Comprehensive Bond Detection

- Hierarchical bond detection with multiple providers:
  - File-based bonds from PDB CONECT records and mmCIF data
  - Template-based detection using standard residue topologies
  - Chemical Component Dictionary (CCD) lookup for ligands
  - Distance-based detection with Van der Waals radii
- Intelligent fallback system ensures complete bond coverage
- Partial processing support for incremental bond detection

### Powerful Selection Language
- MDAnalysis/VMD-inspired syntax for complex atom queries
- Spatial selections with `within`, `around`, and center-of-geometry queries
- Boolean operations (and, or, not) for combining selections
- Residue-based selections with `byres` modifier

### Multi-Format I/O Support
- PDB format with multi-model support and CONECT record parsing
- mmCIF format with chemical bond information extraction
- Auto-detection of single structures vs. trajectories
- String-based parsing for in-memory structure creation

## Installation

```bash
pip install molr
```

For development installation:
```bash
git clone https://github.com/abhishektiwari/molr.git
cd molr
pip install -e .[dev]
```

## Requirements

- Python ≥3.8
- NumPy ≥1.20.0
- SciPy ≥1.7.0 (for spatial indexing)
- pyparsing ≥3.0.0 (for selection language)

## Usage

Please review [Molr documentation](https://hbat.abhishek-tiwari.com/) for more details on how to use Molr for various use cases.

### Quick Example

```python
import molr

# Load a structure
structure = molr.Structure.from_pdb("protein.pdb")

# Detect bonds
bonds = structure.detect_bonds()

# Use selection language
active_site = structure.select("within 5.0 of (resname HIS)")

# Fast spatial queries
neighbors = structure.get_neighbors_within(atom_idx=100, radius=5.0)
```


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing 

See our [contributing guide](CONTRIBUTING.md) and [development guide](https://hbat.abhishek-tiwari.com/development). At a high-level,

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request