.. MolR documentation master file

Welcome to MolR's documentation!
================================

.. image:: https://img.shields.io/github/v/release/abhishektiwari/molr
   :alt: GitHub Release
   :target: https://github.com/abhishektiwari/molr/releases

.. image:: https://img.shields.io/github/actions/workflow/status/abhishektiwari/molr/test.yml?label=tests
   :alt: GitHub Actions Test Workflow Status
   :target: https://github.com/abhishektiwari/molr/actions/workflows/test.yml

.. pypi-shield::
   :project: molr
   :version:

.. pypi-shield::
   :wheel:

.. pypi-shield::
   :py-versions:
   
.. github-shield::
   :username: abhishektiwari
   :repository: molr
   :branch: main
   :last-commit:

.. image:: https://img.shields.io/pypi/status/molr
   :alt: PyPI - Status

.. image:: https://img.shields.io/conda/v/molr/molr
   :alt: Conda Version

.. github-shield::
   :username: abhishektiwari
   :repository: molr
   :license:

MolR (Molecular Realm for Spatial Indexed Structures) is a high-performance Python package that creates a spatial realm for molecular structures, enabling lightning-fast neighbor searches, geometric queries, and spatial operations through integrated KDTree indexing.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   quickstart
   user_guide
   api_reference
   examples

Key Features
------------

* **High-Performance Structure Representation**: NumPy-based Structure class with Structure of Arrays (SoA) design
* **Efficient Spatial Indexing**: Built-in scipy KDTree integration for O(log n) neighbor queries
* **Comprehensive Bond Detection**: Hierarchical system with multiple providers
* **Powerful Selection Language**: MDAnalysis/VMD-inspired syntax for complex atom queries
* **Multi-Format I/O Support**: PDB and mmCIF formats with automatic structure/trajectory detection

Quick Example
-------------

.. code-block:: python

   import molr

   # Load a structure
   structure = molr.Structure.from_pdb("protein.pdb")
   
   # Detect bonds
   bonds = structure.detect_bonds()
   
   # Use selection language
   active_site = structure.select("within 5.0 of (resname HIS)")
   
   # Fast spatial queries
   neighbors = structure.get_neighbors_within(atom_idx=100, radius=5.0)

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`