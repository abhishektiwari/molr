Quick Start
===========

This guide covers the basics of using MolR for molecular structure analysis.

Loading Structures
------------------

MolR can load molecular structures from PDB and mmCIF files:

.. code-block:: python

   import molr

   # Load from PDB file
   structure = molr.Structure.from_pdb("protein.pdb")
   print(f"Loaded {structure.n_atoms} atoms")

   # Load from mmCIF file
   structure = molr.Structure.from_mmcif("structure.cif")

   # Load from string
   pdb_string = """ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N"""
   structure = molr.Structure.from_pdb_string(pdb_string)

Basic Structure Information
---------------------------

Access basic properties of the structure:

.. code-block:: python

   # Number of atoms
   print(f"Number of atoms: {structure.n_atoms}")

   # Atom properties
   print(f"Atom names: {structure.atom_name[:5]}")
   print(f"Elements: {structure.element[:5]}")
   print(f"Coordinates shape: {structure.coord.shape}")

   # Residue information
   print(f"Residue names: {np.unique(structure.res_name)}")
   print(f"Chain IDs: {np.unique(structure.chain_id)}")

Bond Detection
--------------

MolR provides automatic bond detection with multiple methods:

.. code-block:: python

   # Detect bonds automatically
   bonds = structure.detect_bonds()
   print(f"Detected {len(bonds)} bonds")

   # Access bond information
   for i in range(min(5, len(bonds))):
       atom1_idx, atom2_idx = bonds.get_bond(i)
       atom1 = structure.atom_name[atom1_idx]
       atom2 = structure.atom_name[atom2_idx]
       print(f"Bond {i}: {atom1} - {atom2}")

Selection Language
------------------

MolR includes a powerful selection language inspired by MDAnalysis and VMD:

.. code-block:: python

   # Select backbone atoms
   backbone = structure.select("backbone")
   print(f"Selected {len(backbone)} backbone atoms")

   # Select by residue name
   his_residues = structure.select("resname HIS")
   
   # Select by element
   carbons = structure.select("element C")
   
   # Complex selections with boolean operators
   active_site = structure.select("(resname HIS TYR CYS) and (chain A)")
   
   # Spatial selections
   near_ligand = structure.select("within 5.0 of (resname LIG)")

Spatial Queries
---------------

Fast neighbor searches using built-in KDTree indexing:

.. code-block:: python

   # Find neighbors within radius
   atom_idx = 100
   radius = 5.0
   neighbors = structure.get_neighbors_within(atom_idx, radius)
   print(f"Found {len(neighbors)} neighbors within {radius} Å")

   # Find atoms in a sphere
   center = [10.0, 15.0, 20.0]
   atoms_in_sphere = structure.get_atoms_within_sphere(center, radius=8.0)

   # Center of geometry based queries
   protein = structure.select("protein")
   ligand = structure.select("resname LIG")
   nearby = structure.get_atoms_within_cog_sphere(ligand, radius=10.0)

Working with Subsets
--------------------

Create new structures from selections:

.. code-block:: python

   # Extract protein only
   protein_mask = structure.select("protein")
   protein_structure = structure[protein_mask]
   
   # Extract specific chain
   chain_a = structure[structure.chain_id == "A"]
   
   # Combine selections
   backbone_ca = structure.select("backbone and name CA")
   ca_structure = structure[backbone_ca]

Multi-Model Trajectories
------------------------

Handle structures with multiple models:

.. code-block:: python

   # Load trajectory
   ensemble = molr.StructureEnsemble.from_pdb("trajectory.pdb")
   print(f"Loaded {ensemble.n_models} models")

   # Access individual models
   first_model = ensemble[0]
   last_model = ensemble[-1]

   # Iterate over models
   for i, model in enumerate(ensemble):
       center = model.get_center()
       print(f"Model {i} center: {center}")

Simple Analysis Example
-----------------------

Here's a complete example analyzing a protein-ligand complex:

.. code-block:: python

   import molr
   import numpy as np

   # Load structure
   structure = molr.Structure.from_pdb("protein_ligand.pdb")
   
   # Detect bonds
   bonds = structure.detect_bonds()
   
   # Separate protein and ligand
   protein = structure[structure.select("protein")]
   ligand = structure[structure.select("resname LIG")]
   
   # Find binding site residues
   binding_site_mask = structure.select("protein and within 5.0 of (resname LIG)")
   binding_site = structure[binding_site_mask]
   
   # Get unique residues in binding site
   unique_residues = np.unique(binding_site.res_id)
   print(f"Binding site contains {len(unique_residues)} residues")
   
   # Calculate distances
   for res_id in unique_residues[:5]:  # First 5 residues
       res_mask = binding_site.res_id == res_id
       res_atoms = binding_site[res_mask]
       res_name = res_atoms.res_name[0]
       
       # Distance to ligand center
       ligand_center = ligand.get_center()
       res_center = res_atoms.get_center()
       distance = np.linalg.norm(ligand_center - res_center)
       
       print(f"{res_name} {res_id}: {distance:.2f} Å from ligand")

Next Steps
----------

- Explore the :doc:`user_guide` for detailed features
- Check :doc:`examples` for more complex use cases
- See the :doc:`api_reference` for complete API documentation