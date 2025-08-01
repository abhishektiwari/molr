Changelog
=========

All notable changes to MolR are documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Unreleased
----------

Added
~~~~~

* Initial release of MolR package
* NumPy-based Structure class with Structure of Arrays (SoA) design
* Efficient spatial indexing with scipy KDTree integration
* Comprehensive bond detection system with multiple providers:

  * PDB CONECT record parser
  * Residue template-based detection
  * Chemical Component Dictionary (CCD) support
  * Distance-based detection with Van der Waals radii

* Selection language with MDAnalysis/VMD-inspired syntax
* Spatial selection expressions (within, around, cog)
* PDB and mmCIF file format support
* Multi-model trajectory support with StructureEnsemble
* Incremental bond detection with partial processing
* Comprehensive test suite with performance benchmarks
* Full type hints and mypy support
* Complete documentation with Sphinx

Changed
~~~~~~~

* Migrated from molr.space module to standalone molr package

Fixed
~~~~~

* Import paths updated for standalone package structure

0.1.0 - TBD
-----------

Initial public release featuring:

Core Features
~~~~~~~~~~~~~

* **High-Performance Structure Representation**
  
  * NumPy-based Structure class with SoA design
  * Efficient memory layout for vectorized operations
  * Built-in scipy KDTree integration for O(log n) spatial queries
  * Lazy initialization of optional annotations

* **Comprehensive Bond Detection**
  
  * Hierarchical detection system with intelligent fallback
  * Multiple bond providers: file-based, template-based, CCD lookup, distance-based
  * Support for PDB CONECT records and mmCIF bond information
  * Partial processing for incremental bond detection
  * Configurable Van der Waals scaling factors

* **Powerful Selection Language**
  
  * MDAnalysis/VMD-inspired syntax using pyparsing
  * Boolean operations (and, or, not) for complex queries
  * Spatial selections with within, around, and center-of-geometry
  * Residue-based selections with byres modifier
  * Predefined atom groups (protein, backbone, sidechain)

* **Multi-Format I/O Support**
  
  * PDB format with multi-model support and CONECT parsing
  * mmCIF format with chemical bond extraction
  * Automatic detection of single structures vs. trajectories
  * String-based parsing for in-memory structures

* **Advanced Spatial Operations**
  
  * Fast neighbor searches using KDTree indexing
  * Sphere-based atom selection around points or selections
  * Inter-selection contact analysis
  * Center of geometry calculations

Technical Features
~~~~~~~~~~~~~~~~~~

* **Development Tools**
  
  * Modern Python packaging with pyproject.toml
  * Complete type hints with py.typed marker
  * Code quality tools: black, mypy, flake8, isort
  * Comprehensive test suite with pytest
  * Performance benchmarks and integration tests

* **Documentation**
  
  * Complete Sphinx documentation
  * API reference with autodoc
  * User guide with detailed examples
  * Installation and quickstart guides

* **Performance Optimizations**
  
  * Structure of Arrays design for memory efficiency
  * Vectorized NumPy operations throughout
  * Lazy evaluation for optional data
  * Spatial indexing for O(log n) queries vs O(n²) brute force

Dependencies
~~~~~~~~~~~~

* Python ≥ 3.8
* NumPy ≥ 1.20.0
* SciPy ≥ 1.7.0 (KDTree spatial indexing)
* pyparsing ≥ 3.0.0 (selection language)
* pdbreader ≥ 0.1.0 (PDB parsing)
* mmcif ≥ 0.1.0 (mmCIF parsing)

Test Data
~~~~~~~~~

Comprehensive test suite includes:

* Standard protein structures (1crn, 2ptc, 4hhb)
* Structures with hydrogen bonds (6rsa, 2izf)
* Halogen bond examples (4x21, 4laz, 4ub7)
* Multi-model trajectories (1bq0, multi_model.pdb)
* Large structures for performance testing (3j3q)
* Synthetic test cases for edge conditions
* Error handling validation files