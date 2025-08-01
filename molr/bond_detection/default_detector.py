"""
Simplified default bond detector for molecular structures.

This module provides a straightforward bond detection approach:
1. Use file-based bonds if available (from PDB CONECT, mmCIF, etc.)
2. Apply residue templates for standard residues
3. Apply distance-based detection for remaining unbonded atoms
"""

from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..constants.atomic_data import AtomicData
from ..constants.residue_bonds import RESIDUE_BONDS
from ..core.bond_list import BondList, BondOrder
from ..core.structure import Structure


class DefaultBondDetector:
    """
    Simplified bond detector using templates and distance criteria.

    This replaces the complex hierarchical system with a straightforward approach:
    1. Apply residue templates (from residue_bonds.py or CCD)
    2. Apply distance-based detection as fallback
    """

    def __init__(self, vdw_factor: float = 0.75):
        """
        Initialize the default bond detector.

        Args:
            vdw_factor: Factor for Van der Waals radii in distance detection
                       (0.0 < factor <= 1.0). Default 0.75 works well for most cases.
        """
        if not 0.0 < vdw_factor <= 1.0:
            raise ValueError("vdw_factor must be between 0 and 1")
        self.vdw_factor = vdw_factor

    def detect_bonds(
        self, structure: Structure, use_file_bonds: bool = True
    ) -> BondList:
        """
        Detect bonds in a molecular structure.

        Args:
            structure: Structure to analyze
            use_file_bonds: Whether to include file-based bonds (CONECT, etc.)

        Returns:
            BondList containing all detected bonds
        """
        bonds = BondList()
        bonded_atoms: Set[int] = set()

        # Step 1: Include file-based bonds if available and requested
        if use_file_bonds and structure.file_bonds is not None:
            for i in range(len(structure.file_bonds)):
                bond_pair = structure.file_bonds[i]
                if isinstance(bond_pair, tuple):
                    atom1, atom2 = bond_pair[0], bond_pair[1]
                else:
                    continue  # Skip if not a tuple
                # Validate indices
                if 0 <= atom1 < structure.n_atoms and 0 <= atom2 < structure.n_atoms:
                    bonds.add_bond(
                        atom1,
                        atom2,
                        bond_order=structure.file_bonds.bond_order[i],
                        bond_type=structure.file_bonds.bond_type[i],
                        detection_method="file",
                    )
                    bonded_atoms.add(atom1)
                    bonded_atoms.add(atom2)

        # Step 2: Apply residue templates for standard residues
        template_bonds = self._apply_residue_templates(structure, bonded_atoms)
        for bond_info in template_bonds:
            bonds.add_bond(
                bond_info["atom1"],
                bond_info["atom2"],
                bond_order=bond_info["order"],
                bond_type="covalent",
                detection_method="template",
            )
            bonded_atoms.add(bond_info["atom1"])
            bonded_atoms.add(bond_info["atom2"])

        # Step 3: Apply distance-based detection for remaining unbonded atoms
        unbonded_atoms = set(range(structure.n_atoms)) - bonded_atoms
        if unbonded_atoms:
            distance_bonds = self._apply_distance_detection(
                structure, unbonded_atoms, bonded_atoms
            )
            for bond_info in distance_bonds:
                bonds.add_bond(
                    bond_info["atom1"],
                    bond_info["atom2"],
                    bond_order=1.0,
                    bond_type="covalent",
                    detection_method="distance",
                )

        return bonds

    def _apply_residue_templates(
        self, structure: Structure, existing_bonded: Set[int]
    ) -> List[Dict]:
        """
        Apply residue templates to detect bonds in standard residues.

        Args:
            structure: Structure to analyze
            existing_bonded: Set of atoms that already have bonds

        Returns:
            List of bond dictionaries with 'atom1', 'atom2', 'order' keys
        """
        bonds = []

        # Group atoms by residue
        residue_groups = self._group_atoms_by_residue(structure)

        for (res_name, res_id, chain_id), atom_indices in residue_groups.items():
            # Skip if residue not in templates
            if res_name not in RESIDUE_BONDS:
                continue

            template = RESIDUE_BONDS[res_name]

            # Create atom name to index mapping for this residue
            atom_name_to_idx = {}
            for idx in atom_indices:
                atom_name = structure.atom_name[idx]
                atom_name_to_idx[atom_name] = idx

            # Apply template bonds
            for bond in template.get("bonds", []):
                atom1_name = bond["atom1"]
                atom2_name = bond["atom2"]

                # Check if both atoms exist in this residue
                if atom1_name in atom_name_to_idx and atom2_name in atom_name_to_idx:
                    atom1_idx = atom_name_to_idx[atom1_name]
                    atom2_idx = atom_name_to_idx[atom2_name]

                    # Skip if either atom already has bonds from files
                    if atom1_idx in existing_bonded and atom2_idx in existing_bonded:
                        continue

                    # Map bond order
                    order_str = bond.get("order", "sing")
                    bond_order = self._map_bond_order(order_str)

                    bonds.append(
                        {"atom1": atom1_idx, "atom2": atom2_idx, "order": bond_order}
                    )

        return bonds

    def _apply_distance_detection(
        self, structure: Structure, unbonded_atoms: Set[int], existing_bonded: Set[int]
    ) -> List[Dict]:
        """
        Apply distance-based bond detection for unbonded atoms.

        Args:
            structure: Structure to analyze
            unbonded_atoms: Set of atoms without bonds
            existing_bonded: Set of atoms that already have bonds

        Returns:
            List of bond dictionaries with 'atom1', 'atom2' keys
        """
        bonds = []

        # Ensure spatial index exists
        structure._ensure_spatial_index()

        # Process each unbonded atom
        for atom_idx in unbonded_atoms:
            element = structure.element[atom_idx]

            # Get VdW radius
            vdw_radius = AtomicData.VDW_RADII.get(
                element, 1.7
            )  # Default to carbon-like

            # Search radius based on VdW radii
            search_radius = 2.0 * vdw_radius * self.vdw_factor

            # Find neighbors
            neighbors = structure.get_neighbors_within(atom_idx, search_radius)

            for neighbor_idx in neighbors:
                if neighbor_idx <= atom_idx:  # Avoid duplicates
                    continue

                # Get neighbor element and radius
                neighbor_element = structure.element[neighbor_idx]
                neighbor_vdw = AtomicData.VDW_RADII.get(neighbor_element, 1.7)

                # Calculate distance threshold
                threshold = (vdw_radius + neighbor_vdw) * self.vdw_factor

                # Calculate actual distance
                distance = np.linalg.norm(
                    structure.coord[atom_idx] - structure.coord[neighbor_idx]
                )

                # Check if within bonding distance
                if distance <= threshold:
                    bonds.append({"atom1": atom_idx, "atom2": neighbor_idx})

        return bonds

    def _group_atoms_by_residue(self, structure: Structure) -> Dict[Tuple, List[int]]:
        """
        Group atom indices by residue.

        Args:
            structure: Structure to analyze

        Returns:
            Dictionary mapping (res_name, res_id, chain_id) to list of atom indices
        """
        residue_groups: Dict[Tuple[str, int, str], List[int]] = {}

        for i in range(structure.n_atoms):
            key = (structure.res_name[i], structure.res_id[i], structure.chain_id[i])

            if key not in residue_groups:
                residue_groups[key] = []
            residue_groups[key].append(i)

        return residue_groups

    def _map_bond_order(self, order_str: str) -> float:
        """
        Map bond order string to numeric value.

        Args:
            order_str: Bond order string (e.g., 'sing', 'doub', 'trip')

        Returns:
            Numeric bond order
        """
        order_mapping = {
            "sing": 1.0,
            "doub": 2.0,
            "trip": 3.0,
            "arom": 1.5,
            "delo": 1.5,
            "pi": 1.0,
            "quad": 4.0,
        }
        return order_mapping.get(order_str.lower(), 1.0)


def detect_bonds(
    structure: Structure, vdw_factor: float = 0.75, use_file_bonds: bool = True
) -> BondList:
    """
    Convenience function to detect bonds in a structure.

    Args:
        structure: Structure to analyze
        vdw_factor: Factor for VdW radii in distance detection
        use_file_bonds: Whether to include file-based bonds

    Returns:
        BondList with detected bonds
    """
    detector = DefaultBondDetector(vdw_factor=vdw_factor)
    return detector.detect_bonds(structure, use_file_bonds=use_file_bonds)
