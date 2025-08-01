Examples
========

This section provides comprehensive examples demonstrating MolR's capabilities for various molecular analysis tasks.

Basic Structure Analysis
------------------------

Loading and Inspecting Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import molr
   import numpy as np

   # Load a protein structure
   structure = molr.Structure.from_pdb("1crn.pdb")
   
   print(f"Structure has {structure.n_atoms} atoms")
   print(f"Chains present: {np.unique(structure.chain_id)}")
   print(f"Residue types: {np.unique(structure.res_name)}")
   
   # Get basic structural information
   center = structure.get_center()
   print(f"Geometric center: {center}")
   
   # Coordinate statistics
   coord_range = np.ptp(structure.coord, axis=0)
   print(f"Structure dimensions: {coord_range}")

Atom and Residue Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Count atoms by element
   elements, counts = np.unique(structure.element, return_counts=True)
   for element, count in zip(elements, counts):
       print(f"{element}: {count} atoms")
   
   # Residue composition
   residues, res_counts = np.unique(structure.res_name, return_counts=True)
   for residue, count in zip(residues, res_counts):
       print(f"{residue}: {count} residues")
   
   # Chain analysis
   for chain in np.unique(structure.chain_id):
       chain_mask = structure.chain_id == chain
       n_residues = len(np.unique(structure.res_id[chain_mask]))
       n_atoms = np.sum(chain_mask)
       print(f"Chain {chain}: {n_residues} residues, {n_atoms} atoms")

Bond Detection and Analysis
---------------------------

Comprehensive Bond Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from molr.bond_detection import DefaultBondDetector

   # Load structure
   structure = molr.Structure.from_pdb("2ptc.pdb")
   
   # Configure bond detector
   detector = DefaultBondDetector()
   
   # Detect bonds
   bonds = detector.detect_bonds(structure)
   
   print(f"Total bonds detected: {len(bonds)}")

Bond Analysis by Residue Type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze bonds by residue type
   bond_stats = {}
   
   for i in range(len(bonds)):
       atom1_idx, atom2_idx = bonds.get_bond(i)
       
       res1 = structure.res_name[atom1_idx]
       res2 = structure.res_name[atom2_idx]
       
       # Intra-residue bonds
       if res1 == res2:
           key = f"intra-{res1}"
       else:
           # Inter-residue bonds
           key = f"inter-{res1}-{res2}"
       
       bond_stats[key] = bond_stats.get(key, 0) + 1
   
   # Display top bond types
   sorted_bonds = sorted(bond_stats.items(), key=lambda x: x[1], reverse=True)
   print("Top bond types:")
   for bond_type, count in sorted_bonds[:10]:
       print(f"  {bond_type}: {count}")

Protein-Ligand Interaction Analysis
-----------------------------------

Binding Site Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load protein-ligand complex
   structure = molr.Structure.from_pdb("2ptc.pdb")
   
   # Identify ligand and protein
   ligand_mask = structure.select("not protein and not water")
   protein_mask = structure.select("protein")
   
   if np.any(ligand_mask):
       ligand_residues = np.unique(structure.res_name[ligand_mask])
       print(f"Ligand residues: {ligand_residues}")
       
       # Find binding site residues
       binding_site_mask = structure.select(
           "protein and within 5.0 of (not protein and not water)"
       )
       
       # Get unique binding site residues
       binding_residues = []
       for res_id in np.unique(structure.res_id[binding_site_mask]):
           res_mask = structure.res_id == res_id
           res_name = structure.res_name[res_mask][0]
           chain = structure.chain_id[res_mask][0]
           binding_residues.append(f"{res_name}{res_id}{chain}")
       
       print(f"Binding site residues ({len(binding_residues)}):")
       for residue in binding_residues:
           print(f"  {residue}")

Contact Analysis
~~~~~~~~~~~~~~~~

.. code-block:: python

   # Detailed contact analysis
   contacts = structure.get_atoms_between_selections(
       protein_mask, ligand_mask, max_distance=4.0
   )
   
   print(f"Found {len(contacts)} protein-ligand contacts")
   
   # Analyze contact types
   contact_types = {}
   
   for protein_idx, ligand_idx in contacts:
       protein_atom = structure.atom_name[protein_idx]
       protein_res = structure.res_name[protein_idx]
       ligand_atom = structure.atom_name[ligand_idx]
       ligand_res = structure.res_name[ligand_idx]
       
       contact_key = f"{protein_res}:{protein_atom} - {ligand_res}:{ligand_atom}"
       contact_types[contact_key] = contact_types.get(contact_key, 0) + 1
   
   # Display most common contacts
   sorted_contacts = sorted(contact_types.items(), key=lambda x: x[1], reverse=True)
   print("Most common contacts:")
   for contact, count in sorted_contacts[:10]:
       print(f"  {contact}: {count}")

Advanced Selection Examples
---------------------------

Complex Spatial Selections
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load structure
   structure = molr.Structure.from_pdb("4hhb.pdb")  # Hemoglobin
   
   # Find heme groups
   heme_mask = structure.select("resname HEM")
   if np.any(heme_mask):
       print(f"Found {np.sum(heme_mask)} heme atoms")
       
       # Find residues coordinating iron
       iron_mask = structure.select("resname HEM and name FE")
       if np.any(iron_mask):
           # Get iron coordinates
           iron_indices = np.where(iron_mask)[0]
           
           for iron_idx in iron_indices:
               # Find coordinating residues
               coord_mask = structure.select(
                   f"protein and within 3.0 of index {iron_idx}"
               )
               
               coord_residues = []
               for res_id in np.unique(structure.res_id[coord_mask]):
                   res_mask = structure.res_id == res_id
                   res_name = structure.res_name[res_mask][0]
                   chain = structure.chain_id[res_mask][0]
                   coord_residues.append(f"{res_name}{res_id}{chain}")
               
               print(f"Iron {iron_idx} coordinated by: {coord_residues}")

Multi-Chain Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze each chain separately
   chains = np.unique(structure.chain_id)
   
   for chain in chains:
       chain_mask = structure.select(f"chain {chain}")
       chain_structure = structure[chain_mask]
       
       # Chain composition
       residue_count = len(np.unique(chain_structure.res_id))
       atom_count = chain_structure.n_atoms
       
       print(f"Chain {chain}: {residue_count} residues, {atom_count} atoms")
       
       # Secondary structure elements (simplified)
       backbone_mask = chain_structure.select("backbone")
       if np.any(backbone_mask):
           backbone_atoms = chain_structure[backbone_mask]
           # Could add secondary structure analysis here

Interface Analysis
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Find inter-chain contacts
   chain_pairs = []
   chains = np.unique(structure.chain_id)
   
   for i, chain1 in enumerate(chains):
       for chain2 in chains[i+1:]:
           chain1_mask = structure.chain_id == chain1
           chain2_mask = structure.chain_id == chain2
           
           # Find contacts between chains
           contacts = structure.get_atoms_between_selections(
               chain1_mask, chain2_mask, max_distance=4.0
           )
           
           if len(contacts) > 0:
               print(f"Chain {chain1} - Chain {chain2}: {len(contacts)} contacts")
               chain_pairs.append((chain1, chain2, contacts))
   
   # Analyze interface residues
   for chain1, chain2, contacts in chain_pairs:
       interface_residues = set()
       
       for atom1_idx, atom2_idx in contacts:
           res1_id = structure.res_id[atom1_idx]
           res2_id = structure.res_id[atom2_idx]
           
           interface_residues.add((chain1, res1_id))
           interface_residues.add((chain2, res2_id))
       
       print(f"Interface residues ({chain1}-{chain2}): {len(interface_residues)}")

Trajectory Analysis
-------------------

Multi-Model Structure Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load multi-model structure
   ensemble = molr.StructureEnsemble.from_pdb("1bq0.pdb")
   
   print(f"Loaded {ensemble.n_models} models")
   print(f"Each model has {ensemble.n_atoms} atoms")
   
   # Analyze structural variation
   centers = []
   ca_atoms_list = []
   
   for i, model in enumerate(ensemble):
       center = model.get_center()
       centers.append(center)
       
       # Get CA atoms for each model
       ca_mask = model.select("name CA")
       ca_coords = model.coord[ca_mask]
       ca_atoms_list.append(ca_coords)
   
   # Calculate center variations
   centers = np.array(centers)
   center_variation = np.std(centers, axis=0)
   print(f"Center variation (x,y,z): {center_variation}")

RMSD Calculation
~~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate RMSD between models
   if len(ca_atoms_list) > 1:
       reference = ca_atoms_list[0]  # First model as reference
       
       rmsds = []
       for i, coords in enumerate(ca_atoms_list[1:], 1):
           if coords.shape == reference.shape:
               # Simple RMSD calculation (without alignment)
               diff = coords - reference
               rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
               rmsds.append(rmsd)
               print(f"Model {i} RMSD to reference: {rmsd:.2f} Å")
       
       if rmsds:
           print(f"Average RMSD: {np.mean(rmsds):.2f} Å")
           print(f"RMSD range: {np.min(rmsds):.2f} - {np.max(rmsds):.2f} Å")

Structural Flexibility Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze per-residue flexibility
   if len(ca_atoms_list) > 2:
       ca_array = np.array(ca_atoms_list)  # Shape: (n_models, n_residues, 3)
       
       # Calculate per-residue RMSF (Root Mean Square Fluctuation)
       mean_positions = np.mean(ca_array, axis=0)
       fluctuations = []
       
       for i in range(ca_array.shape[1]):  # For each residue
           residue_coords = ca_array[:, i, :]  # All models for this residue
           deviations = residue_coords - mean_positions[i]
           rmsf = np.sqrt(np.mean(np.sum(deviations**2, axis=1)))
           fluctuations.append(rmsf)
       
       fluctuations = np.array(fluctuations)
       
       # Find most and least flexible residues
       max_flex_idx = np.argmax(fluctuations)
       min_flex_idx = np.argmin(fluctuations)
       
       # Get residue information from first model
       first_model = ensemble[0]
       ca_mask = first_model.select("name CA")
       res_ids = first_model.res_id[ca_mask]
       res_names = first_model.res_name[ca_mask]
       
       print(f"Most flexible residue: {res_names[max_flex_idx]}{res_ids[max_flex_idx]} "
             f"(RMSF: {fluctuations[max_flex_idx]:.2f} Å)")
       print(f"Least flexible residue: {res_names[min_flex_idx]}{res_ids[min_flex_idx]} "
             f"(RMSF: {fluctuations[min_flex_idx]:.2f} Å)")

Hydrogen Bond Analysis
----------------------

Simple Hydrogen Bond Detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load structure with hydrogens
   structure = molr.Structure.from_pdb("6rsa.pdb")
   
   # Detect all bonds first
   bonds = structure.detect_bonds()
   
   # Find potential hydrogen bonds (simplified criteria)
   h_atoms = structure.select("element H")
   donors = []
   
   # Find hydrogen bond donors (N-H, O-H)
   for h_idx in np.where(h_atoms)[0]:
       h_neighbors = bonds.get_neighbors(h_idx)
       
       for neighbor_idx in h_neighbors:
           neighbor_element = structure.element[neighbor_idx]
           if neighbor_element in ['N', 'O']:
               donors.append((neighbor_idx, h_idx))
   
   print(f"Found {len(donors)} potential hydrogen bond donors")
   
   # Find acceptors (N, O atoms not bonded to H)
   acceptors = []
   for atom_idx in range(structure.n_atoms):
       element = structure.element[atom_idx]
       if element in ['N', 'O']:
           neighbors = bonds.get_neighbors(atom_idx)
           neighbor_elements = [structure.element[n] for n in neighbors]
           
           # Simple criteria: N/O not saturated with H
           h_count = neighbor_elements.count('H')
           if (element == 'N' and h_count < 3) or (element == 'O' and h_count < 2):
               acceptors.append(atom_idx)
   
   print(f"Found {len(acceptors)} potential hydrogen bond acceptors")

Hydrogen Bond Geometry Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze hydrogen bond geometry
   potential_hbonds = []
   
   for donor_heavy, h_atom in donors[:10]:  # Limit for example
       donor_coord = structure.coord[donor_heavy]
       h_coord = structure.coord[h_atom]
       
       # Find nearby acceptors
       for acceptor in acceptors:
           acceptor_coord = structure.coord[acceptor]
           
           # Distance criteria
           da_distance = np.linalg.norm(donor_coord - acceptor_coord)
           ha_distance = np.linalg.norm(h_coord - acceptor_coord)
           
           if 2.5 <= da_distance <= 3.5 and ha_distance <= 2.5:
               # Angle criteria (simplified)
               v1 = h_coord - donor_coord
               v2 = acceptor_coord - donor_coord
               
               cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
               angle = np.arccos(np.clip(cos_angle, -1, 1)) * 180 / np.pi
               
               if angle > 120:  # Reasonable hydrogen bond angle
                   donor_res = f"{structure.res_name[donor_heavy]}{structure.res_id[donor_heavy]}"
                   acceptor_res = f"{structure.res_name[acceptor]}{structure.res_id[acceptor]}"
                   
                   potential_hbonds.append({
                       'donor': donor_res,
                       'acceptor': acceptor_res,
                       'distance': da_distance,
                       'angle': angle
                   })
   
   print(f"Found {len(potential_hbonds)} potential hydrogen bonds")
   for hb in potential_hbonds[:5]:  # Show first 5
       print(f"  {hb['donor']} -> {hb['acceptor']}: "
             f"{hb['distance']:.2f} Å, {hb['angle']:.1f}°")

Performance Benchmarking
------------------------

Spatial Query Performance
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import time

   # Load large structure
   structure = molr.Structure.from_pdb("3j3q.cif")  # Large ribosome structure
   
   print(f"Benchmarking with {structure.n_atoms} atoms")
   
   # Benchmark neighbor searches
   atom_idx = structure.n_atoms // 2
   radii = [3.0, 5.0, 8.0, 10.0]
   
   for radius in radii:
       start_time = time.time()
       neighbors = structure.get_neighbors_within(atom_idx, radius)
       end_time = time.time()
       
       print(f"Radius {radius} Å: {len(neighbors)} neighbors in {end_time - start_time:.4f}s")
   
   # Benchmark selections
   selections = [
       "protein",
       "backbone",
       "element C",
       "within 5.0 of (chain A)"
   ]
   
   for selection in selections:
       start_time = time.time()
       mask = structure.select(selection)
       end_time = time.time()
       
       print(f"Selection '{selection}': {np.sum(mask)} atoms in {end_time - start_time:.4f}s")

Memory Usage Analysis
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import sys

   # Analyze memory usage of different operations
   def get_size_mb(obj):
       return sys.getsizeof(obj) / 1024 / 1024

   print(f"Structure object: {get_size_mb(structure):.2f} MB")
   print(f"Coordinates: {get_size_mb(structure.coord):.2f} MB")
   print(f"Atom names: {get_size_mb(structure.atom_name):.2f} MB")
   
   # Bond detection memory usage
   bonds = structure.detect_bonds()
   print(f"Bonds object: {get_size_mb(bonds):.2f} MB")
   
   # Selection memory usage
   protein_mask = structure.select("protein")
   print(f"Selection mask: {get_size_mb(protein_mask):.2f} MB")

This comprehensive set of examples demonstrates MolR's capabilities across different use cases. Each example can be adapted and extended for specific research needs.