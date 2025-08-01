"""
Selection expression classes for molecular structure queries.

This module provides the expression system for atom selection language,
inspired by MDAnalysis and VMD selection syntax. Each expression class
represents a specific selection criterion or operation.
"""

from abc import ABC, abstractmethod
from typing import Any, List, Optional, Set, Union, cast

import numpy as np

from ..core.structure import Structure


class SelectionExpression(ABC):
    """
    Abstract base class for all selection expressions.

    Selection expressions form a tree structure that can be evaluated
    against a Structure to produce a boolean mask indicating which atoms
    are selected.
    """

    @abstractmethod
    def evaluate(self, structure: Structure) -> np.ndarray:
        """
        Evaluate the expression against a structure.

        Args:
            structure: The molecular structure to evaluate against

        Returns:
            Boolean array with True for selected atoms
        """
        pass

    def __and__(self, other: "SelectionExpression") -> "SelectionExpression":
        """Create AND expression using & operator."""
        return AndExpression(self, other)

    def __or__(self, other: "SelectionExpression") -> "SelectionExpression":
        """Create OR expression using | operator."""
        return OrExpression(self, other)

    def __invert__(self) -> "SelectionExpression":
        """Create NOT expression using ~ operator."""
        return NotExpression(self)

    @abstractmethod
    def __repr__(self) -> str:
        """String representation of the expression."""
        pass


# Basic Selection Expressions


class AllExpression(SelectionExpression):
    """Select all atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Return True for all atoms."""
        return cast(np.ndarray, np.ones(structure.n_atoms, dtype=bool))

    def __repr__(self) -> str:
        return "all"


class NoneExpression(SelectionExpression):
    """Select no atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Return False for all atoms."""
        return cast(np.ndarray, np.zeros(structure.n_atoms, dtype=bool))

    def __repr__(self) -> str:
        return "none"


class ElementExpression(SelectionExpression):
    """Select atoms by element type."""

    def __init__(self, elements: Union[str, List[str]]):
        """
        Initialize element selection.

        Args:
            elements: Element symbol(s) to select
        """
        if isinstance(elements, str):
            elements = [elements]
        self.elements = [elem.upper() for elem in elements]

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms matching the specified elements."""
        return np.isin(structure.element, self.elements)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if len(self.elements) == 1:
            return f"element {self.elements[0]}"
        return f"element {' '.join(self.elements)}"


class AtomNameExpression(SelectionExpression):
    """Select atoms by atom name."""

    def __init__(self, names: Union[str, List[str]]):
        """
        Initialize atom name selection.

        Args:
            names: Atom name(s) to select
        """
        if isinstance(names, str):
            names = [names]
        self.names = names

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms matching the specified names."""
        return np.isin(structure.atom_name, self.names)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if len(self.names) == 1:
            return f"name {self.names[0]}"
        return f"name {' '.join(self.names)}"


class ResidueNameExpression(SelectionExpression):
    """Select atoms by residue name."""

    def __init__(self, resnames: Union[str, List[str]]):
        """
        Initialize residue name selection.

        Args:
            resnames: Residue name(s) to select
        """
        if isinstance(resnames, str):
            resnames = [resnames]
        self.resnames = resnames

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms in residues matching the specified names."""
        return np.isin(structure.res_name, self.resnames)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if len(self.resnames) == 1:
            return f"resname {self.resnames[0]}"
        return f"resname {' '.join(self.resnames)}"


class ResidueIdExpression(SelectionExpression):
    """Select atoms by residue ID."""

    def __init__(self, resids: Union[int, List[int], range]):
        """
        Initialize residue ID selection.

        Args:
            resids: Residue ID(s) to select
        """
        if isinstance(resids, int):
            resids = [resids]
        elif isinstance(resids, range):
            resids = list(resids)
        self.resids = resids

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms in residues matching the specified IDs."""
        return np.isin(structure.res_id, self.resids)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if len(self.resids) == 1:
            return f"resid {self.resids[0]}"
        # Check if it's a continuous range
        if len(self.resids) > 2:
            sorted_ids = sorted(self.resids)
            if sorted_ids == list(range(sorted_ids[0], sorted_ids[-1] + 1)):
                return f"resid {sorted_ids[0]}:{sorted_ids[-1]}"
        return f"resid {' '.join(map(str, self.resids))}"


class ChainExpression(SelectionExpression):
    """Select atoms by chain ID."""

    def __init__(self, chains: Union[str, List[str]]):
        """
        Initialize chain selection.

        Args:
            chains: Chain ID(s) to select
        """
        if isinstance(chains, str):
            chains = list(chains)  # Split single string into characters
        self.chains = chains

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms in the specified chains."""
        return np.isin(structure.chain_id, self.chains)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if len(self.chains) == 1:
            return f"chain {self.chains[0]}"
        return f"chain {''.join(self.chains)}"


# Property-based Selection Expressions


class BackboneExpression(SelectionExpression):
    """Select backbone atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of the backbone."""
        return structure.is_backbone

    def __repr__(self) -> str:
        return "backbone"


class SidechainExpression(SelectionExpression):
    """Select sidechain atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of sidechains."""
        return structure.is_sidechain

    def __repr__(self) -> str:
        return "sidechain"


class ProteinExpression(SelectionExpression):
    """Select protein atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of protein residues."""
        return structure.residue_type == "PROTEIN"  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return "protein"


class NucleicExpression(SelectionExpression):
    """Select nucleic acid atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of DNA or RNA."""
        return (structure.residue_type == "DNA") | (structure.residue_type == "RNA")  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return "nucleic"


class DNAExpression(SelectionExpression):
    """Select DNA atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of DNA."""
        return structure.residue_type == "DNA"  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return "dna"


class RNAExpression(SelectionExpression):
    """Select RNA atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of RNA."""
        return structure.residue_type == "RNA"  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return "rna"


class LigandExpression(SelectionExpression):
    """Select ligand atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of ligands."""
        return structure.is_ligand

    def __repr__(self) -> str:
        return "ligand"


class AromaticExpression(SelectionExpression):
    """Select aromatic atoms."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of aromatic systems."""
        return structure.is_aromatic

    def __repr__(self) -> str:
        return "aromatic"


class WaterExpression(SelectionExpression):
    """Select water molecules."""

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms that are part of water molecules."""
        from ..constants.pdb_constants import WATER_MOLECULES

        return np.isin(structure.res_name, WATER_MOLECULES)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return "water"


# Boolean Operation Expressions


class AndExpression(SelectionExpression):
    """Logical AND of two expressions."""

    def __init__(self, left: SelectionExpression, right: SelectionExpression):
        """
        Initialize AND expression.

        Args:
            left: Left operand
            right: Right operand
        """
        self.left = left
        self.right = right

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Return atoms selected by both expressions."""
        return self.left.evaluate(structure) & self.right.evaluate(structure)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return f"({self.left!r} and {self.right!r})"


class OrExpression(SelectionExpression):
    """Logical OR of two expressions."""

    def __init__(self, left: SelectionExpression, right: SelectionExpression):
        """
        Initialize OR expression.

        Args:
            left: Left operand
            right: Right operand
        """
        self.left = left
        self.right = right

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Return atoms selected by either expression."""
        return self.left.evaluate(structure) | self.right.evaluate(structure)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return f"({self.left!r} or {self.right!r})"


class NotExpression(SelectionExpression):
    """Logical NOT of an expression."""

    def __init__(self, operand: SelectionExpression):
        """
        Initialize NOT expression.

        Args:
            operand: Expression to negate
        """
        self.operand = operand

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Return atoms not selected by the expression."""
        return ~self.operand.evaluate(structure)

    def __repr__(self) -> str:
        return f"(not {self.operand!r})"


# Index-based Selection Expressions


class IndexExpression(SelectionExpression):
    """Select atoms by index."""

    def __init__(self, indices: Union[int, List[int], range, slice]):
        """
        Initialize index selection.

        Args:
            indices: Atom indices to select
        """
        if isinstance(indices, int):
            indices = [indices]
        elif isinstance(indices, (range, slice)):
            indices = indices
        self.indices = indices

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms at the specified indices."""
        mask = np.zeros(structure.n_atoms, dtype=bool)
        if isinstance(self.indices, slice):
            mask[self.indices] = True
        elif isinstance(self.indices, range):
            mask[list(self.indices)] = True
        else:
            mask[self.indices] = True
        return mask  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        if isinstance(self.indices, slice):
            start = self.indices.start or 0
            stop = self.indices.stop
            step = self.indices.step or 1
            if step == 1:
                return f"index {start}:{stop}"
            return f"index {start}:{stop}:{step}"
        elif isinstance(self.indices, range):
            if self.indices.step == 1:
                return f"index {self.indices.start}:{self.indices.stop}"
            return f"index {self.indices.start}:{self.indices.stop}:{self.indices.step}"
        elif len(self.indices) == 1:
            return f"index {self.indices[0]}"
        return f"index {' '.join(map(str, self.indices))}"


# Composite Selection Expression


class ByResidueExpression(SelectionExpression):
    """Select complete residues based on atom selection."""

    def __init__(self, atom_selection: SelectionExpression):
        """
        Initialize by-residue selection.

        Args:
            atom_selection: Expression to identify residues
        """
        self.atom_selection = atom_selection

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select all atoms in residues that have any selected atoms."""
        # Get atoms selected by the inner expression
        atom_mask = self.atom_selection.evaluate(structure)

        # Find unique residues containing selected atoms
        selected_residues = np.unique(structure.res_id[atom_mask])

        # Select all atoms in those residues
        return np.isin(structure.res_id, selected_residues)  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return f"byres ({self.atom_selection!r})"


# Spatial Selection Expressions


class WithinExpression(SelectionExpression):
    """Select atoms within a distance of another selection."""

    def __init__(self, distance: float, selection: SelectionExpression):
        """
        Initialize within selection.

        Args:
            distance: Maximum distance in Angstroms
            selection: Selection to measure distance from
        """
        self.distance = distance
        self.selection = selection

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms within distance of the selection."""
        # Get reference atoms
        ref_mask = self.selection.evaluate(structure)
        ref_indices = np.where(ref_mask)[0]

        if len(ref_indices) == 0:
            return cast(np.ndarray, np.zeros(structure.n_atoms, dtype=bool))

        # Use Structure's spatial indexing if available
        if structure.has_spatial_index():
            # Get all neighbors within distance for each reference atom
            result_mask = np.zeros(structure.n_atoms, dtype=bool)

            for ref_idx in ref_indices:
                neighbors = structure.get_neighbors_within(ref_idx, self.distance)
                result_mask[neighbors] = True
                result_mask[ref_idx] = True  # Include reference atoms themselves

            return result_mask  # type: ignore[no-any-return]
        else:
            # Fall back to brute force calculation
            result_mask = np.zeros(structure.n_atoms, dtype=bool)
            ref_coords = structure.coord[ref_indices]

            for i in range(structure.n_atoms):
                atom_coord = structure.coord[i]
                distances = np.linalg.norm(ref_coords - atom_coord, axis=1)
                if np.any(distances <= self.distance):
                    result_mask[i] = True

            return result_mask

    def __repr__(self) -> str:
        return f"within {self.distance} of ({self.selection!r})"


class AroundExpression(SelectionExpression):
    """Select atoms around a selection (alternative syntax for within)."""

    def __init__(self, selection: SelectionExpression, distance: float):
        """
        Initialize around selection.

        Args:
            selection: Selection to find atoms around
            distance: Maximum distance in Angstroms
        """
        self.selection = selection
        self.distance = distance

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms around the selection."""
        # Delegate to WithinExpression
        within_expr = WithinExpression(self.distance, self.selection)
        return within_expr.evaluate(structure)

    def __repr__(self) -> str:
        return f"around ({self.selection!r}) {self.distance}"


class SphericalExpression(SelectionExpression):
    """Select atoms within a spherical region."""

    def __init__(self, center: np.ndarray, radius: float):
        """
        Initialize spherical selection.

        Args:
            center: Center point (x, y, z) of sphere
            radius: Radius of sphere in Angstroms
        """
        self.center = np.array(center)
        self.radius = radius

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms within the spherical region."""
        if structure.has_spatial_index():
            # Use Structure's spatial indexing
            atom_indices = structure.get_atoms_within_sphere(self.center, self.radius)
            result_mask = np.zeros(structure.n_atoms, dtype=bool)
            result_mask[atom_indices] = True
            return result_mask  # type: ignore[no-any-return]
        else:
            # Fall back to direct calculation
            distances = np.linalg.norm(structure.coord - self.center, axis=1)
            return distances <= self.radius  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return f"sphere center {self.center} radius {self.radius}"


class CenterOfGeometryExpression(SelectionExpression):
    """Select atoms within distance of center of geometry of a selection."""

    def __init__(self, selection: SelectionExpression, distance: float):
        """
        Initialize COG-based selection.

        Args:
            selection: Selection to calculate COG from
            distance: Distance from COG to select atoms
        """
        self.selection = selection
        self.distance = distance

    def evaluate(self, structure: Structure) -> np.ndarray:
        """Select atoms within distance of selection's center of geometry."""
        # Get reference atoms
        ref_mask = self.selection.evaluate(structure)
        ref_indices = np.where(ref_mask)[0]

        if len(ref_indices) == 0:
            return cast(np.ndarray, np.zeros(structure.n_atoms, dtype=bool))

        if structure.has_spatial_index():
            # Use Structure's spatial indexing
            atom_indices = structure.get_atoms_within_cog_sphere(
                ref_indices, self.distance
            )
            result_mask = np.zeros(structure.n_atoms, dtype=bool)
            result_mask[atom_indices] = True
            return result_mask  # type: ignore[no-any-return]
        else:
            # Fall back to manual COG calculation
            ref_coords = structure.coord[ref_indices]
            cog = np.mean(ref_coords, axis=0)
            distances = np.linalg.norm(structure.coord - cog, axis=1)
            return distances <= self.distance  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        return f"cog ({self.selection!r}) {self.distance}"
