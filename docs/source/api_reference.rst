API Reference
=============

This section provides detailed documentation for all MolR classes and functions.

Core Classes
------------

Structure
~~~~~~~~~

.. autoclass:: molr.Structure
   :members:
   :undoc-members:
   :show-inheritance:

The main class for representing molecular structures with spatial indexing capabilities.

**Key Methods:**

* :meth:`~molr.Structure.from_pdb` - Load from PDB file
* :meth:`~molr.Structure.from_mmcif` - Load from mmCIF file
* :meth:`~molr.Structure.select` - Atom selection using query language
* :meth:`~molr.Structure.detect_bonds` - Automatic bond detection
* :meth:`~molr.Structure.get_neighbors_within` - Spatial neighbor queries

StructureEnsemble
~~~~~~~~~~~~~~~~~

.. autoclass:: molr.StructureEnsemble
   :members:
   :undoc-members:
   :show-inheritance:

Multi-model trajectory representation for handling structural ensembles.

**Key Methods:**

* :meth:`~molr.StructureEnsemble.from_pdb` - Load multi-model PDB
* :meth:`~molr.StructureEnsemble.__getitem__` - Access individual models
* :meth:`~molr.StructureEnsemble.__len__` - Number of models

BondList
~~~~~~~~

.. autoclass:: molr.BondList
   :members:
   :undoc-members:
   :show-inheritance:

Efficient storage and manipulation of molecular bonds.

**Key Methods:**

* :meth:`~molr.BondList.get_bond` - Get bond between atoms
* :meth:`~molr.BondList.get_neighbors` - Get bonded neighbors
* :meth:`~molr.BondList.to_connectivity_matrix` - Convert to adjacency matrix

Bond Detection
--------------

DefaultBondDetector
~~~~~~~~~~~~~~~~~~~

.. autoclass:: molr.bond_detection.DefaultBondDetector
   :members:
   :undoc-members:
   :show-inheritance:

Default bond detector that combines residue templates and distance-based detection.

Bond Detection Functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: molr.bond_detection.detect_bonds

Main function for bond detection in molecular structures.

I/O Parsers
-----------

PDB Parser
~~~~~~~~~~

.. autoclass:: molr.PDBParser
   :members:
   :undoc-members:
   :show-inheritance:

Parser for PDB format files with support for:

* Multi-model structures
* CONECT record parsing
* Alternate conformations
* Insertion codes
* Crystal information

mmCIF Parser
~~~~~~~~~~~~

.. autoclass:: molr.mmCIFParser
   :members:
   :undoc-members:
   :show-inheritance:

Parser for mmCIF format files with support for:

* Chemical bond information
* Large structure handling
* Complete metadata extraction

Selection System
----------------

Selection Engine
~~~~~~~~~~~~~~~~

.. autoclass:: molr.selection.SelectionEngine
   :members:
   :undoc-members:
   :show-inheritance:

Main engine for parsing and evaluating selection expressions.

Selection Functions
~~~~~~~~~~~~~~~~~~~

.. autofunction:: molr.selection.select

Main selection function for atom queries.

.. autofunction:: molr.selection.select_atoms

Alternative selection function.

Selection Parser
~~~~~~~~~~~~~~~~

.. autoclass:: molr.selection.SelectionParser
   :members:
   :undoc-members:
   :show-inheritance:

pyparsing-based parser for selection language syntax.

**Supported Expressions:**

* Atom properties: `name`, `element`, `resname`, `chain`
* Spatial queries: `within`, `around`, `cog`
* Boolean operations: `and`, `or`, `not`
* Residue modifiers: `byres`
* Predefined groups: `protein`, `backbone`, `sidechain`

Expression Classes
~~~~~~~~~~~~~~~~~~

Base Expression
^^^^^^^^^^^^^^^

.. autoclass:: molr.selection.SelectionExpression
   :members:
   :undoc-members:
   :show-inheritance:

Atom Property Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: molr.selection.ElementExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.AtomNameExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.ResidueNameExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.ResidueIdExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.ChainExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.IndexExpression
   :members:
   :undoc-members:
   :show-inheritance:

Structural Expressions
^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: molr.selection.BackboneExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.SidechainExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.ProteinExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.NucleicExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.DNAExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.RNAExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.LigandExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.AromaticExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.WaterExpression
   :members:
   :undoc-members:
   :show-inheritance:

Boolean Expressions
^^^^^^^^^^^^^^^^^^^

.. autoclass:: molr.selection.AndExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.OrExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.NotExpression
   :members:
   :undoc-members:
   :show-inheritance:

Special Expressions
^^^^^^^^^^^^^^^^^^^

.. autoclass:: molr.selection.AllExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.NoneExpression
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: molr.selection.ByResidueExpression
   :members:
   :undoc-members:
   :show-inheritance:

Utilities
---------

Atom Utilities
~~~~~~~~~~~~~~

.. automodule:: molr.utilities.atom_utils
   :members:
   :undoc-members:

**Key Functions:**

* :func:`~molr.utilities.atom_utils.classify_atom_types` - Classify backbone/sidechain
* :func:`~molr.utilities.atom_utils.calculate_center_of_geometry` - COG calculation
* :func:`~molr.utilities.atom_utils.get_vdw_radius` - Van der Waals radii lookup

Constants
---------

Atomic Data
~~~~~~~~~~~

.. automodule:: molr.constants.atomic_data
   :members:
   :undoc-members:

Standard atomic properties and constants.

**Available Data:**

* Element symbols and atomic numbers
* Van der Waals radii
* Covalent radii
* Atomic masses

Bond Parameters
~~~~~~~~~~~~~~~

.. automodule:: molr.constants.bond_parameters
   :members:
   :undoc-members:

Bond length and angle parameters for different atom types.

**Available Data:**

* Standard bond lengths
* Bond angle preferences
* Distance cutoffs for bond detection

PDB Constants
~~~~~~~~~~~~~

.. automodule:: molr.constants.pdb_constants
   :members:
   :undoc-members:

PDB format constants and mappings.

**Available Data:**

* Record type definitions
* Standard residue names
* Chain identifier mappings

Residue Bond Templates
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molr.constants.residue_bonds
   :members:
   :undoc-members:

Standard topology templates for common residues.

**Available Templates:**

* Amino acid topologies
* Nucleotide topologies
* Common ligand templates
* Metal coordination patterns

Configuration
-------------

.. automodule:: molr.config
   :members:
   :undoc-members:

Global configuration settings for MolR.

**Configuration Options:**

* Default bond detection parameters
* Spatial indexing settings
* I/O parser options
* Selection language settings


Type Hints
----------

MolR provides complete type hint coverage. Key type aliases:

.. code-block:: python

   from typing import Union, List, Tuple, Optional
   import numpy as np
   
   # Common type aliases used throughout MolR
   AtomIndex = int
   AtomMask = np.ndarray  # Boolean array for atom selection
   Coordinates = np.ndarray  # Shape (n_atoms, 3)
   BondPair = Tuple[AtomIndex, AtomIndex]
   SelectionString = str

Usage Examples
--------------

Here are some common usage patterns for the API:

.. code-block:: python

   import molr
   import numpy as np
   
   # Load and analyze structure
   structure = molr.Structure.from_pdb("protein.pdb")
   bonds = structure.detect_bonds()
   
   # Selection operations
   backbone = structure.select("backbone")
   active_site = structure.select("within 5.0 of (resname LIG)")
   
   # Spatial queries
   neighbors = structure.get_neighbors_within(100, 5.0)
   sphere_atoms = structure.get_atoms_within_sphere([0, 0, 0], 10.0)
   
   # Bond analysis
   atom_neighbors = bonds.get_neighbors(100)
   connectivity = bonds.to_connectivity_matrix(structure.n_atoms)
   
   # Custom bond detection
   from molr.bond_detection import DefaultBondDetector
   detector = DefaultBondDetector()
   custom_bonds = detector.detect_bonds(structure)

For more detailed examples, see the :doc:`examples` section.