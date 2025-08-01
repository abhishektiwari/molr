"""
Structure class for NumPy-based molecular structure representation.

This module provides the Structure class,
for efficient storage and manipulation of molecular structures using NumPy arrays.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    from scipy.spatial import KDTree

    HAS_SCIPY = True
except ImportError:
    KDTree = None
    HAS_SCIPY = False

from ..constants.atomic_data import AtomicData
from ..constants.pdb_constants import (
    DNA_RESIDUES,
    DNA_RNA_BACKBONE_ATOMS,
    DNA_RNA_BASE_ATOMS,
    PROTEIN_BACKBONE_ATOMS,
    PROTEIN_RESIDUES,
    PROTEIN_SIDECHAIN_ATOMS,
    RESIDUES_WITH_AROMATIC_RINGS,
    RNA_RESIDUES,
    WATER_MOLECULES,
)

if TYPE_CHECKING:
    from .bond_list import BondList


class Structure:
    """
    A NumPy-based molecular structure representation using Structure of Arrays (SoA) design.

    This class stores molecular data in separate NumPy arrays
    for efficient vectorized operations and memory usage. Supports optional annotations
    with lazy initialization to save memory when not needed.

    Core annotations (always present):
        - coord: Atomic coordinates (x, y, z) as float64
        - atom_name: PDB atom names as U4 strings
        - element: Element symbols as U2 strings
        - res_name: Residue names as U3 strings
        - res_id: Residue sequence numbers as int32
        - chain_id: Chain identifiers as U1 strings

    Optional annotations (lazy initialization):
        - alt_loc: Alternate location indicators
        - occupancy: Occupancy values
        - b_factor: Temperature factors
        - charge: Formal charges
        - serial: Atom serial numbers
        - insertion_code: Residue insertion codes
        - segment_id: Segment identifiers

    Classification flags (computed on demand):
        - is_backbone: Boolean array for backbone atoms
        - is_sidechain: Boolean array for sidechain atoms
        - is_aromatic: Boolean array for aromatic atoms
        - is_ligand: Boolean array for ligand atoms
        - residue_type: Residue type classification

    Example:
        >>> structure = Structure(n_atoms=100)
        >>> structure.coord = np.random.rand(100, 3)
        >>> structure.atom_name[:] = "CA"
        >>> structure.add_annotation("custom_prop", dtype=np.float32, default_value=1.0)
    """

    def __init__(self, n_atoms: int):
        """
        Initialize Structure with core annotations only.

        Args:
            n_atoms: Number of atoms in the structure

        Raises:
            ValueError: If n_atoms <= 0
        """
        if n_atoms <= 0:
            raise ValueError("Number of atoms must be positive")

        self.n_atoms = n_atoms

        # Core annotations - always present
        self._coord = np.zeros((n_atoms, 3), dtype=np.float64)
        self.atom_name = np.empty(n_atoms, dtype="U4")
        self.element = np.empty(n_atoms, dtype="U2")
        self.res_name = np.empty(n_atoms, dtype="U3")
        self.res_id = np.zeros(n_atoms, dtype=np.int32)
        self.chain_id = np.empty(n_atoms, dtype="U1")

        # Optional annotations - lazy initialization
        self.alt_loc: Optional[np.ndarray] = None
        self.occupancy: Optional[np.ndarray] = None
        self.b_factor: Optional[np.ndarray] = None
        self.charge: Optional[np.ndarray] = None
        self.serial: Optional[np.ndarray] = None
        self.insertion_code: Optional[np.ndarray] = None
        self.segment_id: Optional[np.ndarray] = None

        # Classification flags - computed on demand
        self._is_backbone: Optional[np.ndarray] = None
        self._is_sidechain: Optional[np.ndarray] = None
        self._is_aromatic: Optional[np.ndarray] = None
        self._is_ligand: Optional[np.ndarray] = None
        self._residue_type: Optional[np.ndarray] = None

        # Bond storage
        self._bonds: Optional["BondList"] = None
        self._file_bonds: Optional["BondList"] = None  # Bonds from PDB/mmCIF files

        # Spatial indexing - lazy initialization
        self._spatial_index: Optional["KDTree"] = None
        self._coords_hash: Optional[int] = None

        # Track custom annotations
        self._custom_annotations: set = set()

    @property
    def coord(self) -> np.ndarray:
        """
        Atomic coordinates array (n_atoms, 3).

        Returns:
            NumPy array of atomic coordinates
        """
        return self._coord  # type: ignore[no-any-return]

    @coord.setter
    def coord(self, value: np.ndarray) -> None:
        """
        Set atomic coordinates and invalidate spatial index.

        Args:
            value: New coordinate array (n_atoms, 3)
        """
        if value.shape != (self.n_atoms, 3):
            raise ValueError(
                f"Coordinates must have shape ({self.n_atoms}, 3), got {value.shape}"
            )
        self._coord = value.astype(np.float64)
        self._invalidate_spatial_index()

    def add_annotation(
        self, name: str, dtype: Any = np.float32, default_value: Any = None
    ) -> None:
        """
        Add custom annotation to structure.

        Args:
            name: Name of the annotation
            dtype: NumPy data type for the annotation
            default_value: Default value to fill the array (optional)

        Raises:
            ValueError: If annotation name already exists
        """
        if hasattr(self, name):
            raise ValueError(f"Annotation '{name}' already exists")

        array = np.empty(self.n_atoms, dtype=dtype)
        if default_value is not None:
            array.fill(default_value)

        setattr(self, name, array)
        self._custom_annotations.add(name)

    def _ensure_annotation(
        self, annotation_name: str, dtype: Any, default_value: Any = None
    ) -> np.ndarray:
        """
        Ensure optional annotation exists, creating it if necessary.

        Args:
            annotation_name: Name of the annotation attribute
            dtype: NumPy data type
            default_value: Default value to fill array

        Returns:
            The annotation array
        """
        annotation = getattr(self, annotation_name)
        if annotation is None:
            annotation = np.empty(self.n_atoms, dtype=dtype)
            if default_value is not None:
                annotation.fill(default_value)
            setattr(self, annotation_name, annotation)
        return annotation  # type: ignore

    @property
    def is_backbone(self) -> np.ndarray:
        """Boolean array indicating backbone atoms."""
        if self._is_backbone is None:
            self._update_classification_flags()
        return self._is_backbone  # type: ignore

    @property
    def is_sidechain(self) -> np.ndarray:
        """Boolean array indicating sidechain atoms."""
        if self._is_sidechain is None:
            self._update_classification_flags()
        return self._is_sidechain  # type: ignore

    @property
    def is_aromatic(self) -> np.ndarray:
        """Boolean array indicating aromatic atoms."""
        if self._is_aromatic is None:
            self._update_classification_flags()
        return self._is_aromatic  # type: ignore

    @property
    def is_ligand(self) -> np.ndarray:
        """Boolean array indicating ligand atoms."""
        if self._is_ligand is None:
            self._update_classification_flags()
        return self._is_ligand  # type: ignore

    @property
    def residue_type(self) -> np.ndarray:
        """String array with residue type classification (PROTEIN/DNA/RNA/LIGAND)."""
        if self._residue_type is None:
            self._update_classification_flags()
        return self._residue_type  # type: ignore

    def _update_classification_flags(self) -> None:
        """Update all classification flags based on current structure data."""
        # Initialize arrays
        self._is_backbone = np.zeros(self.n_atoms, dtype=bool)
        self._is_sidechain = np.zeros(self.n_atoms, dtype=bool)
        self._is_aromatic = np.zeros(self.n_atoms, dtype=bool)
        self._is_ligand = np.zeros(self.n_atoms, dtype=bool)
        self._residue_type = np.full(self.n_atoms, "LIGAND", dtype="U7")

        # Classify residue types first
        protein_mask = np.isin(self.res_name, PROTEIN_RESIDUES)
        dna_mask = np.isin(self.res_name, DNA_RESIDUES)
        rna_mask = np.isin(self.res_name, RNA_RESIDUES)
        water_mask = np.isin(self.res_name, WATER_MOLECULES)

        self._residue_type[protein_mask] = "PROTEIN"
        self._residue_type[dna_mask] = "DNA"
        self._residue_type[rna_mask] = "RNA"

        # Water is not considered ligand for our purposes
        ligand_mask = ~(protein_mask | dna_mask | rna_mask | water_mask)
        self._is_ligand[ligand_mask] = True

        # Backbone classification
        protein_backbone_mask = protein_mask & np.isin(
            self.atom_name, PROTEIN_BACKBONE_ATOMS
        )
        nucleic_backbone_mask = (dna_mask | rna_mask) & np.isin(
            self.atom_name, DNA_RNA_BACKBONE_ATOMS
        )
        self._is_backbone[protein_backbone_mask | nucleic_backbone_mask] = True

        # Sidechain classification
        protein_sidechain_mask = protein_mask & np.isin(
            self.atom_name, PROTEIN_SIDECHAIN_ATOMS
        )
        nucleic_base_mask = (dna_mask | rna_mask) & np.isin(
            self.atom_name, DNA_RNA_BASE_ATOMS
        )
        self._is_sidechain[protein_sidechain_mask | nucleic_base_mask] = True

        # Aromatic classification - atoms in aromatic residues
        aromatic_residue_mask = np.isin(self.res_name, RESIDUES_WITH_AROMATIC_RINGS)
        self._is_aromatic[aromatic_residue_mask] = True

    def __getitem__(self, index: Union[int, slice, np.ndarray]) -> "Structure":
        """
        Get subset of structure by index.

        Args:
            index: Integer, slice, or boolean/integer array for indexing

        Returns:
            New Structure containing selected atoms

        Example:
            >>> subset = structure[structure.element == "C"]
            >>> single_atom = structure[0]
            >>> chain_a = structure[structure.chain_id == "A"]
        """
        if isinstance(index, int):
            index = slice(index, index + 1)
        elif isinstance(index, np.ndarray) and index.dtype == bool:
            index = np.where(index)[0]

        # Determine number of atoms in subset
        if isinstance(index, slice):
            n_subset = len(range(*index.indices(self.n_atoms)))
        else:
            n_subset = len(index)  # type: ignore

        # Create new structure
        subset = Structure(n_subset)

        # Copy core annotations
        subset.coord = self.coord[index].copy()
        subset.atom_name = self.atom_name[index].copy()
        subset.element = self.element[index].copy()
        subset.res_name = self.res_name[index].copy()
        subset.res_id = self.res_id[index].copy()
        subset.chain_id = self.chain_id[index].copy()

        # Copy optional annotations if they exist
        for attr_name in [
            "alt_loc",
            "occupancy",
            "b_factor",
            "charge",
            "serial",
            "insertion_code",
            "segment_id",
        ]:
            attr_value = getattr(self, attr_name)
            if attr_value is not None:
                setattr(subset, attr_name, attr_value[index].copy())

        # Copy custom annotations
        for annotation_name in self._custom_annotations:
            annotation = getattr(self, annotation_name)
            setattr(subset, annotation_name, annotation[index].copy())
            subset._custom_annotations.add(annotation_name)

        return subset

    def __len__(self) -> int:
        """Return number of atoms in structure."""
        return self.n_atoms

    def copy(self) -> "Structure":
        """
        Create deep copy of structure.

        Returns:
            New Structure with copied data
        """
        return self[:]

    def get_center(self, weights: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Calculate geometric or mass-weighted center of structure.

        Args:
            weights: Optional weights for each atom (e.g., atomic masses)

        Returns:
            3D coordinate of center as numpy array
        """
        if weights is None:
            return np.mean(self.coord, axis=0)  # type: ignore
        else:
            if len(weights) != self.n_atoms:
                raise ValueError(
                    "Weights array must have same length as number of atoms"
                )
            return np.average(self.coord, axis=0, weights=weights)  # type: ignore

    def get_masses(self) -> np.ndarray:
        """
        Get atomic masses for all atoms.

        Returns:
            Array of atomic masses in amu
        """
        masses = np.zeros(self.n_atoms)
        for i, element in enumerate(self.element):
            masses[i] = AtomicData.ATOMIC_MASSES.get(
                element.upper(), AtomicData.DEFAULT_ATOMIC_MASS
            )
        return masses  # type: ignore

    def translate(self, vector: np.ndarray) -> None:
        """
        Translate structure by given vector.

        Args:
            vector: 3D translation vector
        """
        if len(vector) != 3:
            raise ValueError("Translation vector must be 3D")
        self.coord += vector
        self._invalidate_spatial_index()

    def center_at_origin(self, weights: Optional[np.ndarray] = None) -> None:
        """
        Center structure at origin.

        Args:
            weights: Optional weights for center calculation
        """
        center = self.get_center(weights)
        self.translate(-center)

    def _ensure_spatial_index(self) -> "KDTree":
        """
        Ensure spatial index is built and up-to-date.

        Returns:
            KDTree spatial index for this structure

        Raises:
            ImportError: If scipy is not available
        """
        if not HAS_SCIPY:
            raise ImportError(
                "scipy is required for spatial indexing. Install with: pip install scipy"
            )

        # Calculate hash of current coordinates
        current_hash = hash(self.coord.tobytes())

        # Build/rebuild index if needed
        if self._spatial_index is None or self._coords_hash != current_hash:
            self._spatial_index = KDTree(self.coord)
            self._coords_hash = current_hash

        return self._spatial_index

    def _invalidate_spatial_index(self) -> None:
        """Invalidate spatial index after coordinate changes."""
        self._spatial_index = None
        self._coords_hash = None

    def get_neighbors_within(self, atom_idx: int, radius: float) -> np.ndarray:
        """
        Get atom indices within radius of specified atom.

        Args:
            atom_idx: Index of query atom
            radius: Search radius in Angstroms

        Returns:
            Array of neighbor atom indices (excluding query atom)

        Example:
            >>> neighbors = structure.get_neighbors_within(100, 5.0)
        """
        if atom_idx < 0 or atom_idx >= self.n_atoms:
            raise IndexError(f"Atom index {atom_idx} out of range")

        spatial_index = self._ensure_spatial_index()
        query_point = self.coord[atom_idx]
        neighbor_indices = spatial_index.query_ball_point(query_point, radius)

        # Remove query atom from results
        return np.array([idx for idx in neighbor_indices if idx != atom_idx])  # type: ignore[no-any-return]

    def get_atoms_within_sphere(self, center: np.ndarray, radius: float) -> np.ndarray:
        """
        Get atoms within spherical region.

        Args:
            center: Center point of sphere (x, y, z)
            radius: Radius of sphere in Angstroms

        Returns:
            Array of atom indices within the sphere

        Example:
            >>> center = np.array([10.0, 15.0, 20.0])
            >>> atoms = structure.get_atoms_within_sphere(center, 8.0)
        """
        if len(center) != 3:
            raise ValueError("Center must be 3D coordinate")

        spatial_index = self._ensure_spatial_index()
        neighbor_indices = spatial_index.query_ball_point(center, radius)
        return np.array(neighbor_indices)  # type: ignore[no-any-return]

    def get_atoms_within_cog_sphere(
        self, selection: np.ndarray, radius: float
    ) -> np.ndarray:
        """
        Get atoms within spherical zone centered at center of geometry of selection.

        Args:
            selection: Boolean mask or indices of atoms to define COG
            radius: Radius of spherical zone in Angstroms

        Returns:
            Array of atom indices within the COG sphere

        Example:
            >>> active_site = structure.select("resname HIS")
            >>> nearby = structure.get_atoms_within_cog_sphere(active_site, 10.0)
        """
        if len(selection) == 0:
            return np.array([])  # type: ignore[no-any-return]

        # Handle both boolean masks and index arrays
        if selection.dtype == bool:
            if len(selection) != self.n_atoms:
                raise ValueError("Boolean selection must match number of atoms")
            selected_coords = self.coord[selection]
        else:
            # Assume it's an array of indices
            selected_coords = self.coord[selection]

        if len(selected_coords) == 0:
            return np.array([])  # type: ignore[no-any-return]

        # Calculate center of geometry
        cog = np.mean(selected_coords, axis=0)

        # Find atoms within sphere around COG
        return self.get_atoms_within_sphere(cog, radius)

    def get_neighbors_for_atoms(
        self, atom_indices: np.ndarray, radius: float
    ) -> Dict[int, np.ndarray]:
        """
        Get neighbors for multiple atoms at once (batch operation).

        Args:
            atom_indices: Array of atom indices to query
            radius: Search radius in Angstroms

        Returns:
            Dictionary mapping atom_idx -> array of neighbor indices

        Example:
            >>> ca_atoms = structure.select("name CA")
            >>> neighbors = structure.get_neighbors_for_atoms(ca_atoms, 8.0)
        """
        spatial_index = self._ensure_spatial_index()
        query_points = self.coord[atom_indices]
        neighbor_lists = spatial_index.query_ball_point(query_points, radius)

        result = {}
        for i, atom_idx in enumerate(atom_indices):
            # Remove query atom from neighbors
            neighbors = [idx for idx in neighbor_lists[i] if idx != atom_idx]
            result[atom_idx] = np.array(neighbors)

        return result

    def get_closest_atoms(
        self, query_point: np.ndarray, k: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get k nearest atoms to a query point.

        Args:
            query_point: 3D coordinate to query
            k: Number of nearest neighbors to return

        Returns:
            Tuple of (distances, atom_indices) for k nearest atoms

        Example:
            >>> center = np.array([0.0, 0.0, 0.0])
            >>> distances, indices = structure.get_closest_atoms(center, k=5)
        """
        if len(query_point) != 3:
            raise ValueError("Query point must be 3D coordinate")

        spatial_index = self._ensure_spatial_index()
        distances, indices = spatial_index.query(query_point, k=min(k, self.n_atoms))

        # Ensure we return arrays even for k=1
        if k == 1:
            distances = np.array([distances])
            indices = np.array([indices])

        return distances, indices

    def get_atoms_between_selections(
        self, selection1: np.ndarray, selection2: np.ndarray, max_distance: float
    ) -> Dict[str, np.ndarray]:
        """
        Find atoms from two selections within max_distance of each other.

        Args:
            selection1: First selection (boolean mask or indices)
            selection2: Second selection (boolean mask or indices)
            max_distance: Maximum distance between selections

        Returns:
            Dictionary with 'selection1_atoms', 'selection2_atoms', 'distances'

        Example:
            >>> protein = structure.select("protein")
            >>> ligand = structure.select("resname LIG")
            >>> contacts = structure.get_atoms_between_selections(protein, ligand, 5.0)
        """
        # Convert selections to index arrays if needed
        if selection1.dtype == bool:
            sel1_indices = np.where(selection1)[0]
        else:
            sel1_indices = selection1

        if selection2.dtype == bool:
            sel2_indices = np.where(selection2)[0]
        else:
            sel2_indices = selection2

        spatial_index = self._ensure_spatial_index()

        # Find all pairs within max_distance
        sel1_contacts = []
        sel2_contacts = []
        distances = []

        for atom1_idx in sel1_indices:
            query_point = self.coord[atom1_idx]
            neighbors = spatial_index.query_ball_point(query_point, max_distance)

            for neighbor_idx in neighbors:
                if neighbor_idx in sel2_indices:
                    distance = np.linalg.norm(
                        self.coord[atom1_idx] - self.coord[neighbor_idx]
                    )
                    sel1_contacts.append(atom1_idx)
                    sel2_contacts.append(neighbor_idx)
                    distances.append(distance)

        return {
            "selection1_atoms": np.array(sel1_contacts),
            "selection2_atoms": np.array(sel2_contacts),
            "distances": np.array(distances),
        }

    def has_spatial_index(self) -> bool:
        """
        Check if spatial index is available.

        Returns:
            True if scipy is available and spatial indexing is possible
        """
        return HAS_SCIPY

    def get_bonds_to(
        self, other_atoms: np.ndarray, max_distance: float = 2.0
    ) -> np.ndarray:
        """
        Find potential bonds to other atoms based on distance.

        Args:
            other_atoms: Indices of other atoms to check
            max_distance: Maximum distance for bond consideration

        Returns:
            Boolean array indicating which atoms have potential bonds
        """
        # Simple distance-based bond detection
        bonds = np.zeros(self.n_atoms, dtype=bool)

        for other_idx in other_atoms:
            distances = np.linalg.norm(self.coord - self.coord[other_idx], axis=1)
            # Exclude self-bonds
            distances[other_idx] = np.inf
            bonds |= distances <= max_distance

        return bonds  # type: ignore

    def select(self, selection_string: str) -> np.ndarray:
        """
        Select atoms using selection language.

        Args:
            selection_string: Selection expression

        Returns:
            Boolean array of selected atoms

        Examples:
            >>> mask = structure.select("protein and backbone")
            >>> mask = structure.select("resname ALA GLY")
            >>> mask = structure.select("chain A and resid 1:50")

        Raises:
            NotImplementedError: For unsupported selection syntax
        """
        from pyparsing import ParseException

        from ..selection import select

        try:
            return select(self, selection_string)
        except ParseException as e:
            # Convert parser errors to NotImplementedError for unsupported selections
            raise NotImplementedError(f"Unsupported selection syntax: {e}") from e

    def get_annotation_info(self) -> Dict[str, Dict[str, Any]]:
        """
        Get information about all annotations in the structure.

        Returns:
            Dictionary with annotation info including dtype and whether it's initialized
        """
        info = {}

        # Core annotations
        core_annotations = [
            "coord",
            "atom_name",
            "element",
            "res_name",
            "res_id",
            "chain_id",
        ]
        for name in core_annotations:
            attr = getattr(self, name)
            info[name] = {
                "type": "core",
                "dtype": attr.dtype,
                "shape": attr.shape,
                "initialized": True,
            }

        # Optional annotations
        optional_annotations = [
            "alt_loc",
            "occupancy",
            "b_factor",
            "charge",
            "serial",
            "insertion_code",
            "segment_id",
        ]
        for name in optional_annotations:
            attr = getattr(self, name)
            info[name] = {
                "type": "optional",
                "dtype": attr.dtype if attr is not None else None,
                "shape": attr.shape if attr is not None else None,
                "initialized": attr is not None,
            }

        # Custom annotations
        for name in self._custom_annotations:
            attr = getattr(self, name)
            info[name] = {
                "type": "custom",
                "dtype": attr.dtype,
                "shape": attr.shape,
                "initialized": True,
            }

        return info

    def __repr__(self) -> str:
        """String representation of Structure."""
        return f"Structure(n_atoms={self.n_atoms})"

    def __str__(self) -> str:
        """Detailed string representation."""
        lines = [f"Structure with {self.n_atoms} atoms"]

        # Show residue type distribution
        if self._residue_type is not None:
            unique_types, counts = np.unique(self.residue_type, return_counts=True)
            type_info = ", ".join(
                [f"{count} {type_}" for type_, count in zip(unique_types, counts)]
            )
            lines.append(f"Residue types: {type_info}")

        # Show chain distribution
        unique_chains, counts = np.unique(self.chain_id, return_counts=True)
        if len(unique_chains) > 1:
            chain_info = ", ".join(
                [
                    f"{count} in chain {chain}"
                    for chain, count in zip(unique_chains, counts)
                ]
            )
            lines.append(f"Chains: {chain_info}")

        # Show annotation status
        annotation_info = self.get_annotation_info()
        optional_initialized = sum(
            1
            for info in annotation_info.values()
            if info["type"] == "optional" and info["initialized"]
        )
        custom_count = len(self._custom_annotations)

        lines.append(
            f"Annotations: {optional_initialized} optional, {custom_count} custom"
        )

        return "\n".join(lines)

    @classmethod
    def from_pdb(cls, filename: str) -> "Structure":
        """
        Create Structure from PDB file.

        Args:
            filename: Path to PDB file

        Returns:
            Structure object with all atoms and annotations

        Raises:
            ValueError: If PDB file contains multiple models

        Example:
            >>> structure = Structure.from_pdb("example.pdb")
            >>> print(f"Loaded {structure.n_atoms} atoms")
        """
        from ..io.pdb import PDBParser

        parser = PDBParser()
        result = parser.parse_file(filename)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "PDB file contains multiple models, use StructureEnsemble.from_pdb() instead"
            )

    @classmethod
    def from_mmcif(cls, filename: str) -> "Structure":
        """
        Create Structure from mmCIF file.

        Args:
            filename: Path to mmCIF file

        Returns:
            Structure object with all atoms and annotations

        Raises:
            ValueError: If mmCIF file contains multiple models

        Example:
            >>> structure = Structure.from_mmcif("example.cif")
            >>> print(f"Loaded {structure.n_atoms} atoms")
        """
        from ..io.mmcif import mmCIFParser

        parser = mmCIFParser()
        result = parser.parse_file(filename)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "mmCIF file contains multiple models, use StructureEnsemble.from_mmcif() instead"
            )

    @classmethod
    def from_pdb_string(cls, pdb_content: str) -> "Structure":
        """
        Create Structure from PDB content string.

        Args:
            pdb_content: PDB file content as string

        Returns:
            Structure object with all atoms and annotations

        Raises:
            ValueError: If PDB content contains multiple models

        Example:
            >>> pdb_data = "ATOM      1  N   ALA A   1      20.154  16.967  22.478  1.00 10.00           N"
            >>> structure = Structure.from_pdb_string(pdb_data)
        """
        from ..io.pdb import PDBParser

        parser = PDBParser()
        result = parser.parse_string(pdb_content)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "PDB content contains multiple models, use StructureEnsemble.from_pdb_string() instead"
            )

    @classmethod
    def from_mmcif_string(cls, mmcif_content: str) -> "Structure":
        """
        Create Structure from mmCIF content string.

        Args:
            mmcif_content: mmCIF file content as string

        Returns:
            Structure object with all atoms and annotations

        Raises:
            ValueError: If mmCIF content contains multiple models

        Example:
            >>> mmcif_data = "data_test\\nloop_\\n_atom_site.group_PDB\\n..."
            >>> structure = Structure.from_mmcif_string(mmcif_data)
        """
        from ..io.mmcif import mmCIFParser

        parser = mmCIFParser()
        result = parser.parse_string(mmcif_content)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "mmCIF content contains multiple models, use StructureEnsemble.from_mmcif_string() instead"
            )

    def detect_bonds(
        self,
        vdw_factor: float = 0.75,
        use_file_bonds: bool = True,
        store_bonds: bool = True,
    ) -> "BondList":
        """
        Detect bonds using the simplified default detector.

        Args:
            vdw_factor: Factor for VdW radii in distance detection (0.0 < factor <= 1.0)
            use_file_bonds: Whether to include file-based bonds (CONECT, mmCIF)
            store_bonds: Whether to store detected bonds on structure

        Returns:
            BondList with detected bonds

        Example:
            >>> structure = Structure.from_pdb("protein.pdb")
            >>> bonds = structure.detect_bonds()
            >>> print(f"Detected {len(bonds)} bonds")
        """
        from ..bond_detection import DefaultBondDetector

        detector = DefaultBondDetector(vdw_factor=vdw_factor)
        bonds = detector.detect_bonds(self, use_file_bonds=use_file_bonds)

        if store_bonds:
            self._bonds = bonds

        return bonds

    @property
    def bonds(self) -> Optional["BondList"]:
        """
        Get bonds associated with this structure.

        Returns:
            BondList if bonds have been detected/assigned, None otherwise
        """
        return getattr(self, "_bonds", None)

    @bonds.setter
    def bonds(self, bond_list: "BondList") -> None:
        """
        Set bonds for this structure.

        Args:
            bond_list: BondList to associate with structure
        """
        self._bonds = bond_list

    def has_bonds(self) -> bool:
        """
        Check if structure has bond information.

        Returns:
            True if bonds are available
        """
        bonds = self.bonds
        return bonds is not None and len(bonds) > 0

    @property
    def file_bonds(self) -> Optional["BondList"]:
        """
        Get bonds loaded from structure files (PDB CONECT, mmCIF bonds).

        Returns:
            BondList with file-based bonds or None if not available
        """
        return self._file_bonds

    @file_bonds.setter
    def file_bonds(self, bond_list: "BondList") -> None:
        """
        Set file-based bonds for this structure.

        Args:
            bond_list: BondList with file-based bonds
        """
        self._file_bonds = bond_list

    def has_file_bonds(self) -> bool:
        """
        Check if structure has file-based bond information.

        Returns:
            True if file bonds are available
        """
        file_bonds = self.file_bonds
        return file_bonds is not None and len(file_bonds) > 0
