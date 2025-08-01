User Guide
==========

This comprehensive guide covers all major features of MolR for molecular structure analysis.

Structure Representation
-------------------------

Understanding the Structure Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolR uses a Structure of Arrays (SoA) design for efficient memory usage and vectorized operations:

.. code-block:: python

   import molr
   import numpy as np

   structure = molr.Structure.from_pdb("protein.pdb")
   
   # All atom properties are NumPy arrays
   print(f"Coordinates shape: {structure.coord.shape}")
   print(f"Atom names: {structure.atom_name}")
   print(f"Elements: {structure.element}")
   print(f"Residue IDs: {structure.res_id}")

Key Properties
~~~~~~~~~~~~~~

The Structure class provides access to standard molecular properties:

.. code-block:: python

   # Atomic properties
   structure.atom_name     # Atom names (N, CA, C, O, etc.)
   structure.element       # Chemical elements (C, N, O, etc.)
   structure.coord         # 3D coordinates (n_atoms, 3)
   structure.charge        # Partial charges
   structure.occupancy     # Crystallographic occupancy
   structure.b_factor      # Temperature factors
   
   # Residue properties
   structure.res_name      # Residue names (ALA, GLY, etc.)
   structure.res_id        # Residue numbers
   structure.chain_id      # Chain identifiers
   structure.insertion_code  # PDB insertion codes
   
   # Structural hierarchy
   structure.n_atoms       # Number of atoms
   structure.get_center()  # Geometric center

File I/O Operations
-------------------

Loading Structures
~~~~~~~~~~~~~~~~~~

MolR supports multiple file formats and loading methods:

.. code-block:: python

   # From PDB file
   structure = molr.Structure.from_pdb("protein.pdb")
   
   # From mmCIF file
   structure = molr.Structure.from_mmcif("structure.cif")
   
   # From string data
   pdb_data = open("protein.pdb").read()
   structure = molr.Structure.from_pdb_string(pdb_data)
   
   # Multi-model files (trajectories)
   ensemble = molr.StructureEnsemble.from_pdb("trajectory.pdb")

Handling Multi-Model Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolR automatically detects whether a file contains single or multiple models:

.. code-block:: python

   # Single model returns Structure
   single_model = molr.Structure.from_pdb("1crn.pdb")
   
   # Multi-model returns StructureEnsemble
   ensemble = molr.StructureEnsemble.from_pdb("multi_model.pdb")
   
   # Access individual models
   first_model = ensemble[0]    # Returns Structure
   last_model = ensemble[-1]
   
   # Iterate over models
   for i, model in enumerate(ensemble):
       print(f"Model {i}: {model.n_atoms} atoms")

Bond Detection System
---------------------

Hierarchical Bond Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolR uses a sophisticated hierarchical system with multiple providers:

.. code-block:: python

   from molr.bond_detection import DefaultBondDetector
   
   # Default detector
   detector = DefaultBondDetector()
   
   # Detect bonds
   bonds = detector.detect_bonds(structure)
   print(f"Detected {len(bonds)} bonds")

Bond Provider Details
~~~~~~~~~~~~~~~~~~~~~

The DefaultBondDetector combines multiple detection methods:

1. **File-based bonds**: From PDB CONECT records and mmCIF chemical bonds
2. **Template-based**: Standard residue topologies (amino acids, nucleotides)
3. **Distance-based**: Fallback using Van der Waals radii

.. code-block:: python

   # Access different bond sources
   file_bonds = structure.file_bonds    # From original file
   
   # Complete bond set
   all_bonds = structure.bonds

Working with Bonds
~~~~~~~~~~~~~~~~~~

The BondList class provides efficient bond storage and analysis:

.. code-block:: python

   bonds = structure.detect_bonds()
   
   # Basic bond information
   print(f"Number of bonds: {len(bonds)}")
   
   # Iterate over bonds
   for i in range(len(bonds)):
       atom1_idx, atom2_idx = bonds.get_bond(i)
       atom1_name = structure.atom_name[atom1_idx]
       atom2_name = structure.atom_name[atom2_idx]
       print(f"Bond {i}: {atom1_name}-{atom2_name}")
   
   # Get neighbors of an atom
   neighbors = bonds.get_neighbors(atom_idx=100)
   
   # Bond connectivity matrix
   connectivity = bonds.to_connectivity_matrix(structure.n_atoms)

Selection Language
------------------

Basic Selections
~~~~~~~~~~~~~~~~

MolR provides a powerful selection language inspired by MDAnalysis and VMD:

.. code-block:: python

   # Select by atom name
   ca_atoms = structure.select("name CA")
   backbone = structure.select("backbone")
   
   # Select by residue
   his_residues = structure.select("resname HIS")
   protein = structure.select("protein")
   
   # Select by element
   carbons = structure.select("element C")
   nitrogens = structure.select("element N")
   
   # Select by chain
   chain_a = structure.select("chain A")

Boolean Operations
~~~~~~~~~~~~~~~~~~

Combine selections with logical operators:

.. code-block:: python

   # AND operation
   his_ca = structure.select("resname HIS and name CA")
   
   # OR operation
   aromatics = structure.select("resname HIS or resname PHE or resname TYR")
   
   # NOT operation
   non_hydrogen = structure.select("not element H")
   
   # Complex combinations
   active_site = structure.select("(resname HIS TYR CYS) and (chain A) and backbone")

Spatial Selections
~~~~~~~~~~~~~~~~~~

Select atoms based on spatial relationships:

.. code-block:: python

   # Within distance of selection
   near_ligand = structure.select("within 5.0 of (resname LIG)")
   
   # Around a point
   center_atoms = structure.select("around 8.0 of (10.0, 15.0, 20.0)")
   
   # Center of geometry based
   ligand = structure.select("resname LIG")
   binding_site = structure.select("protein and within 6.0 of cog (resname LIG)")

Residue-Based Selections
~~~~~~~~~~~~~~~~~~~~~~~~

Use the ``byres`` modifier to select entire residues:

.. code-block:: python

   # Select entire residues containing CA atoms within 5Å of ligand
   binding_residues = structure.select("byres (name CA and within 5.0 of (resname LIG))")
   
   # Select residues with any atom near the binding site
   contact_residues = structure.select("byres (protein and within 4.0 of (resname LIG))")

Spatial Indexing and Queries
-----------------------------

KDTree Integration
~~~~~~~~~~~~~~~~~~

MolR uses scipy's KDTree for efficient spatial queries:

.. code-block:: python

   # Built-in KDTree is automatically created and cached
   atom_idx = 100
   radius = 5.0
   
   # Find neighbors within radius (O(log n) complexity)
   neighbors = structure.get_neighbors_within(atom_idx, radius)
   print(f"Found {len(neighbors)} neighbors")

Advanced Spatial Queries
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Atoms within sphere around point
   center = [10.0, 15.0, 20.0]
   sphere_atoms = structure.get_atoms_within_sphere(center, radius=8.0)
   
   # Atoms within sphere around center of geometry
   ligand_mask = structure.select("resname LIG")
   nearby_atoms = structure.get_atoms_within_cog_sphere(ligand_mask, radius=10.0)
   
   # Contact analysis between selections
   protein_mask = structure.select("protein")
   ligand_mask = structure.select("resname LIG")
   contacts = structure.get_atoms_between_selections(
       protein_mask, ligand_mask, max_distance=4.0
   )

Structure Manipulation
----------------------

Creating Substructures
~~~~~~~~~~~~~~~~~~~~~~~

Extract parts of structures based on selections:

.. code-block:: python

   # Create new structure from selection
   protein_mask = structure.select("protein")
   protein_only = structure[protein_mask]
   
   # Chain extraction
   chain_a_mask = structure.chain_id == "A"
   chain_a = structure[chain_a_mask]
   
   # Combine multiple criteria
   ca_atoms_mask = structure.select("name CA")
   ca_structure = structure[ca_atoms_mask]

Coordinate Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Translation
   structure.translate([10.0, 0.0, 0.0])
   
   # Center at origin
   structure.center_at_origin()
   
   # Custom transformations
   import numpy as np
   
   # Rotate around z-axis
   angle = np.pi / 4
   rotation_matrix = np.array([
       [np.cos(angle), -np.sin(angle), 0],
       [np.sin(angle), np.cos(angle), 0],
       [0, 0, 1]
   ])
   structure.coord = structure.coord @ rotation_matrix.T

Adding Annotations
~~~~~~~~~~~~~~~~~~

Extend structures with custom properties:

.. code-block:: python

   # Add custom annotation
   structure.add_annotation("hydrophobicity", dtype=np.float32, default_value=0.0)
   
   # Set values
   hydrophobic_mask = structure.select("resname ALA VAL LEU ILE PHE TRP")
   structure.hydrophobicity[hydrophobic_mask] = 1.0
   
   # Use in selections
   hydrophobic_atoms = structure.select("hydrophobicity > 0.5")

Performance Optimization
------------------------

Memory Management
~~~~~~~~~~~~~~~~~

MolR uses lazy initialization for optional annotations:

.. code-block:: python

   # Annotations are loaded only when accessed
   structure = molr.Structure.from_pdb("large_protein.pdb")
   
   # These properties are always available (core data)
   coords = structure.coord  # Always loaded
   names = structure.atom_name  # Always loaded
   
   # These are loaded on first access (optional data)
   charges = structure.charge  # Loaded when first accessed
   bfactors = structure.b_factor  # Loaded when first accessed

Efficient Selections
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Use numpy operations for simple selections
   ca_mask = structure.atom_name == "CA"  # Fast numpy comparison
   
   # Selection language for complex queries
   complex_mask = structure.select("resname HIS and within 5.0 of (resname LIG)")
   
   # Cache selections for reuse
   protein_mask = structure.select("protein")
   # Reuse protein_mask multiple times instead of re-selecting

Working with Large Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For very large structures (>100k atoms), consider subsetting first
   structure = molr.Structure.from_pdb("large_complex.pdb")
   
   # Extract region of interest first
   roi_mask = structure.select("chain A and within 20.0 of (resname LIG)")
   roi_structure = structure[roi_mask]
   
   # Then perform expensive operations on smaller subset
   bonds = roi_structure.detect_bonds()
   neighbors = roi_structure.get_neighbors_within(atom_idx=10, radius=5.0)

Error Handling and Validation
------------------------------

Common Issues and Solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   try:
       structure = molr.Structure.from_pdb("protein.pdb")
   except FileNotFoundError:
       print("PDB file not found")
   except ValueError as e:
       print(f"Invalid PDB format: {e}")
   
   # Validate structure
   if structure.n_atoms == 0:
       print("Warning: No atoms loaded")
   
   # Check for missing coordinates
   if np.any(np.isnan(structure.coord)):
       print("Warning: NaN coordinates detected")

Best Practices
--------------

1. **Use appropriate file formats**: PDB for simple structures, mmCIF for large complexes
2. **Cache expensive operations**: Store bond detection and selection results
3. **Subset large structures**: Work with regions of interest when possible
4. **Validate input data**: Check for missing atoms, coordinates, or bonds
5. **Use spatial indexing**: Leverage built-in KDTree for neighbor searches
6. **Combine selections efficiently**: Use boolean operations instead of multiple passes