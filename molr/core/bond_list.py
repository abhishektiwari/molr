"""
BondList class for efficient bond storage and manipulation.

This module provides the BondList class for storing molecular bonds with smart indexing
that automatically adjusts when the parent structure is modified.
"""

from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


class BondOrder(Enum):
    """Bond order enumeration."""

    SINGLE = 1.0
    DOUBLE = 2.0
    TRIPLE = 3.0
    AROMATIC = 1.5
    UNKNOWN = 0.0


class BondList:
    """
    Efficient storage and manipulation of molecular bonds with smart indexing.

    The BondList class stores bonds as pairs of atom indices with additional metadata
    such as bond order, detection method, and confidence scores. It supports smart
    indexing that automatically adjusts bond indices when the parent structure is
    sliced or modified.

    Bond storage uses Structure of Arrays (SoA) design:
        - bonds: (N, 2) array of atom index pairs
        - bond_order: Bond order (1=single, 2=double, 3=triple, 1.5=aromatic)
        - bond_type: Bond type classification
        - detection_method: How the bond was detected
        - confidence: Confidence score for bond existence

    Smart indexing features:
        - Automatic bond index adjustment when structure is sliced
        - Efficient bond filtering based on atom selections
        - Bond validation against structure changes

    Example:
        >>> bond_list = BondList()
        >>> bond_list.add_bond(0, 1, bond_order=1.0, bond_type="covalent")
        >>> bond_list.add_bonds([(2, 3), (3, 4)], bond_orders=[1.0, 2.0])
        >>> subset_bonds = bond_list.filter_by_atoms([0, 1, 2])
    """

    def __init__(self, n_bonds: int = 0):
        """
        Initialize BondList.

        Args:
            n_bonds: Initial number of bonds (default: 0 for dynamic growth)
        """
        self.n_bonds = n_bonds
        self._capacity = max(n_bonds, 10)  # Minimum capacity for dynamic growth

        # Core bond data - always present
        self.bonds = np.full((self._capacity, 2), -1, dtype=np.int32)
        self.bond_order = np.ones(self._capacity, dtype=np.float32)
        self.bond_type = np.full(self._capacity, "covalent", dtype="U16")

        # Optional bond metadata - lazy initialization
        self.detection_method: Optional[np.ndarray] = None
        self.confidence: Optional[np.ndarray] = None
        self.distance: Optional[np.ndarray] = None
        self.is_hydrogen_bond: Optional[np.ndarray] = None
        self.is_halogen_bond: Optional[np.ndarray] = None

        # Track custom bond properties
        self._custom_properties: set = set()

        # Index mapping for smart indexing (atom_index -> new_index)
        self._index_mapping: Optional[Dict[int, int]] = None

    def _ensure_capacity(self, required_capacity: int) -> None:
        """
        Ensure arrays have sufficient capacity, growing if necessary.

        Args:
            required_capacity: Required minimum capacity
        """
        if required_capacity <= self._capacity:
            return

        # Grow by 50% or to required capacity, whichever is larger
        new_capacity = max(int(self._capacity * 1.5), required_capacity)

        # Resize core arrays
        new_bonds = np.full((new_capacity, 2), -1, dtype=np.int32)
        new_bonds[: self._capacity] = self.bonds
        self.bonds = new_bonds

        new_bond_order = np.ones(new_capacity, dtype=np.float32)
        new_bond_order[: self._capacity] = self.bond_order
        self.bond_order = new_bond_order

        new_bond_type = np.full(new_capacity, "covalent", dtype="U16")
        new_bond_type[: self._capacity] = self.bond_type
        self.bond_type = new_bond_type

        # Resize optional arrays if they exist
        for attr_name in [
            "detection_method",
            "confidence",
            "distance",
            "is_hydrogen_bond",
            "is_halogen_bond",
        ]:
            attr_value = getattr(self, attr_name)
            if attr_value is not None:
                if attr_name == "detection_method":
                    new_array = np.full(new_capacity, "unknown", dtype="U16")
                elif attr_name in ["is_hydrogen_bond", "is_halogen_bond"]:
                    new_array = np.zeros(new_capacity, dtype=bool)
                else:
                    new_array = np.zeros(new_capacity, dtype=np.float32)
                new_array[: self._capacity] = attr_value
                setattr(self, attr_name, new_array)

        # Resize custom properties
        for prop_name in self._custom_properties:
            prop_array = getattr(self, prop_name)
            new_array = np.zeros(new_capacity, dtype=prop_array.dtype)
            new_array[: self._capacity] = prop_array
            setattr(self, prop_name, new_array)

        self._capacity = new_capacity

    def _ensure_property(
        self, property_name: str, dtype: Any, default_value: Any = None
    ) -> np.ndarray:
        """
        Ensure optional property exists, creating it if necessary.

        Args:
            property_name: Name of the property attribute
            dtype: NumPy data type
            default_value: Default value to fill array

        Returns:
            The property array
        """
        prop_array = getattr(self, property_name)
        if prop_array is None:
            prop_array = np.zeros(self._capacity, dtype=dtype)
            if default_value is not None:
                prop_array[: self.n_bonds].fill(default_value)
            setattr(self, property_name, prop_array)
        return prop_array  # type: ignore

    def add_property(
        self, name: str, dtype: Any = np.float32, default_value: Any = None
    ) -> None:
        """
        Add custom property to bonds.

        Args:
            name: Name of the property
            dtype: NumPy data type for the property
            default_value: Default value to fill existing bonds

        Raises:
            ValueError: If property name already exists
        """
        if hasattr(self, name):
            raise ValueError(f"Property '{name}' already exists")

        prop_array = np.zeros(self._capacity, dtype=dtype)
        if default_value is not None:
            prop_array[: self.n_bonds].fill(default_value)

        setattr(self, name, prop_array)
        self._custom_properties.add(name)

    def add_bond(
        self,
        atom1: int,
        atom2: int,
        bond_order: float = 1.0,
        bond_type: str = "covalent",
        **kwargs: Any,
    ) -> int:
        """
        Add a single bond.

        Args:
            atom1: Index of first atom
            atom2: Index of second atom
            bond_order: Bond order (1.0=single, 2.0=double, etc.)
            bond_type: Type of bond
            **kwargs: Additional bond properties

        Returns:
            Index of the added bond

        Raises:
            ValueError: If atoms are the same or invalid
        """
        if atom1 == atom2:
            raise ValueError("Cannot create bond between same atom")
        if atom1 < 0 or atom2 < 0:
            raise ValueError("Atom indices must be non-negative")

        # Ensure capacity
        self._ensure_capacity(self.n_bonds + 1)

        # Store bond with consistent ordering (smaller index first)
        if atom1 > atom2:
            atom1, atom2 = atom2, atom1

        bond_idx = self.n_bonds
        self.bonds[bond_idx] = [atom1, atom2]
        self.bond_order[bond_idx] = bond_order
        self.bond_type[bond_idx] = bond_type

        # Set optional properties if provided
        for prop_name, value in kwargs.items():
            if hasattr(self, prop_name):
                # Use _ensure_property for lazy initialization of optional properties
                if prop_name == "detection_method":
                    prop_array = self._ensure_property(
                        "detection_method", dtype="U10", default_value="unknown"
                    )
                    prop_array[bond_idx] = value
                else:
                    prop_array = getattr(self, prop_name)
                    if prop_array is not None:
                        prop_array[bond_idx] = value

        self.n_bonds += 1
        return bond_idx

    def add_bonds(
        self,
        bond_pairs: List[Tuple[int, int]],
        bond_orders: Optional[List[float]] = None,
        bond_types: Optional[List[str]] = None,
        **kwargs: Any,
    ) -> np.ndarray:
        """
        Add multiple bonds efficiently.

        Args:
            bond_pairs: List of (atom1, atom2) tuples
            bond_orders: Optional list of bond orders (default: all 1.0)
            bond_types: Optional list of bond types (default: all "covalent")
            **kwargs: Additional properties as lists

        Returns:
            Array of bond indices for added bonds

        Raises:
            ValueError: If list lengths don't match
        """
        n_new_bonds = len(bond_pairs)
        if n_new_bonds == 0:
            return np.array([], dtype=np.int32)  # type: ignore

        # Validate input lengths
        if bond_orders is not None and len(bond_orders) != n_new_bonds:
            raise ValueError("bond_orders length must match bond_pairs length")
        if bond_types is not None and len(bond_types) != n_new_bonds:
            raise ValueError("bond_types length must match bond_pairs length")

        for prop_name, values in kwargs.items():
            if len(values) != n_new_bonds:
                raise ValueError(f"{prop_name} length must match bond_pairs length")

        # Ensure capacity
        self._ensure_capacity(self.n_bonds + n_new_bonds)

        # Prepare bond indices
        start_idx = self.n_bonds
        end_idx = start_idx + n_new_bonds
        bond_indices = np.arange(start_idx, end_idx)

        # Process bond pairs with consistent ordering
        bonds_array = np.array(bond_pairs, dtype=np.int32)
        # Sort each bond pair so smaller index comes first
        bonds_array.sort(axis=1)

        # Validate bonds
        if np.any(bonds_array[:, 0] == bonds_array[:, 1]):
            raise ValueError("Cannot create bonds between same atoms")
        if np.any(bonds_array < 0):
            raise ValueError("Atom indices must be non-negative")

        # Store bonds
        self.bonds[start_idx:end_idx] = bonds_array

        # Store bond orders
        if bond_orders is not None:
            self.bond_order[start_idx:end_idx] = bond_orders
        else:
            self.bond_order[start_idx:end_idx] = 1.0

        # Store bond types
        if bond_types is not None:
            self.bond_type[start_idx:end_idx] = bond_types
        else:
            self.bond_type[start_idx:end_idx] = "covalent"

        # Store optional properties
        for prop_name, values in kwargs.items():
            if hasattr(self, prop_name):
                prop_array = getattr(self, prop_name)
                if prop_array is not None:
                    prop_array[start_idx:end_idx] = values

        self.n_bonds += n_new_bonds
        return bond_indices  # type: ignore

    def remove_bonds(
        self, bond_indices: Union[int, List[int], np.ndarray[Any, Any]]
    ) -> None:
        """
        Remove bonds by index.

        Args:
            bond_indices: Bond index or array of bond indices to remove
        """
        if isinstance(bond_indices, int):
            bond_indices = [bond_indices]
        elif isinstance(bond_indices, np.ndarray):
            bond_indices = bond_indices.tolist()

        # Validate indices
        bond_indices = [idx for idx in bond_indices if 0 <= idx < self.n_bonds]  # type: ignore
        if not bond_indices:
            return

        # Sort in descending order for safe removal
        bond_indices = sorted(set(bond_indices), reverse=True)

        # Remove bonds by shifting arrays
        for idx in bond_indices:
            # Shift remaining bonds down
            self.bonds[idx : self.n_bonds - 1] = self.bonds[idx + 1 : self.n_bonds]
            self.bond_order[idx : self.n_bonds - 1] = self.bond_order[
                idx + 1 : self.n_bonds
            ]
            self.bond_type[idx : self.n_bonds - 1] = self.bond_type[
                idx + 1 : self.n_bonds
            ]

            # Shift optional properties
            for attr_name in [
                "detection_method",
                "confidence",
                "distance",
                "is_hydrogen_bond",
                "is_halogen_bond",
            ]:
                attr_value = getattr(self, attr_name)
                if attr_value is not None:
                    attr_value[idx : self.n_bonds - 1] = attr_value[
                        idx + 1 : self.n_bonds
                    ]

            # Shift custom properties
            for prop_name in self._custom_properties:
                prop_array = getattr(self, prop_name)
                prop_array[idx : self.n_bonds - 1] = prop_array[idx + 1 : self.n_bonds]

            self.n_bonds -= 1

    def get_bonds_for_atom(self, atom_index: int) -> np.ndarray:
        """
        Get all bonds involving a specific atom.

        Args:
            atom_index: Index of the atom

        Returns:
            Array of bond indices involving the atom
        """
        if self.n_bonds == 0:
            return np.array([], dtype=np.int32)  # type: ignore

        active_bonds = self.bonds[: self.n_bonds]
        mask = (active_bonds[:, 0] == atom_index) | (active_bonds[:, 1] == atom_index)
        return np.where(mask)[0]  # type: ignore

    def get_neighbors(self, atom_index: int) -> np.ndarray:
        """
        Get neighbor atoms for a specific atom.

        Args:
            atom_index: Index of the atom

        Returns:
            Array of neighbor atom indices
        """
        bond_indices = self.get_bonds_for_atom(atom_index)
        if len(bond_indices) == 0:
            return np.array([], dtype=np.int32)  # type: ignore

        bonds = self.bonds[bond_indices]
        # Get the other atom in each bond
        neighbors = np.where(bonds[:, 0] == atom_index, bonds[:, 1], bonds[:, 0])
        return neighbors  # type: ignore

    def filter_by_atoms(self, atom_indices: Union[List[int], np.ndarray]) -> "BondList":
        """
        Create new BondList containing only bonds between specified atoms.

        Args:
            atom_indices: List or array of atom indices to keep

        Returns:
            New BondList with filtered bonds and remapped indices
        """
        atom_set = set(atom_indices)

        # Find bonds where both atoms are in the selection
        if self.n_bonds == 0:
            return BondList(0)

        active_bonds = self.bonds[: self.n_bonds]
        mask = np.array([both_atoms_in_set(bond, atom_set) for bond in active_bonds])

        if not np.any(mask):
            return BondList(0)

        # Create index mapping
        index_mapping = {
            old_idx: new_idx for new_idx, old_idx in enumerate(sorted(atom_indices))
        }

        # Create new BondList
        selected_indices = np.where(mask)[0]
        new_bond_list = BondList(len(selected_indices))

        # Copy bonds with remapped indices
        for new_idx, old_idx in enumerate(selected_indices):
            old_bond = active_bonds[old_idx]
            new_bond = [index_mapping[old_bond[0]], index_mapping[old_bond[1]]]
            new_bond_list.bonds[new_idx] = new_bond
            new_bond_list.bond_order[new_idx] = self.bond_order[old_idx]
            new_bond_list.bond_type[new_idx] = self.bond_type[old_idx]

            # Copy optional properties
            for attr_name in [
                "detection_method",
                "confidence",
                "distance",
                "is_hydrogen_bond",
                "is_halogen_bond",
            ]:
                old_attr = getattr(self, attr_name)
                if old_attr is not None:
                    new_attr = new_bond_list._ensure_property(attr_name, old_attr.dtype)
                    new_attr[new_idx] = old_attr[old_idx]

            # Copy custom properties
            for prop_name in self._custom_properties:
                old_prop = getattr(self, prop_name)
                if not hasattr(new_bond_list, prop_name):
                    new_bond_list.add_property(prop_name, old_prop.dtype)
                new_prop = getattr(new_bond_list, prop_name)
                new_prop[new_idx] = old_prop[old_idx]

        new_bond_list._custom_properties = self._custom_properties.copy()
        return new_bond_list

    def get_bond_matrix(self, n_atoms: int) -> np.ndarray:
        """
        Create bond adjacency matrix.

        Args:
            n_atoms: Total number of atoms in structure

        Returns:
            (n_atoms, n_atoms) boolean adjacency matrix
        """
        matrix = np.zeros((n_atoms, n_atoms), dtype=bool)

        if self.n_bonds > 0:
            active_bonds = self.bonds[: self.n_bonds]
            valid_mask = np.all(active_bonds < n_atoms, axis=1)
            valid_bonds = active_bonds[valid_mask]

            matrix[valid_bonds[:, 0], valid_bonds[:, 1]] = True
            matrix[valid_bonds[:, 1], valid_bonds[:, 0]] = True

        return matrix  # type: ignore

    def validate_bonds(self, n_atoms: int) -> Tuple[bool, List[int]]:
        """
        Validate that all bonds reference valid atom indices.

        Args:
            n_atoms: Number of atoms in the structure

        Returns:
            Tuple of (all_valid, list_of_invalid_bond_indices)
        """
        if self.n_bonds == 0:
            return True, []

        active_bonds = self.bonds[: self.n_bonds]
        invalid_mask = np.any((active_bonds < 0) | (active_bonds >= n_atoms), axis=1)
        invalid_indices = np.where(invalid_mask)[0].tolist()

        return len(invalid_indices) == 0, invalid_indices

    def __len__(self) -> int:
        """Return number of bonds."""
        return self.n_bonds

    def __getitem__(
        self, index: Union[int, slice, np.ndarray]
    ) -> Union[Tuple[int, int], "BondList"]:
        """
        Get bond(s) by index.

        Args:
            index: Integer, slice, or array for indexing

        Returns:
            Single bond tuple or new BondList with selected bonds
        """
        if isinstance(index, int):
            if index < 0 or index >= self.n_bonds:
                raise IndexError("Bond index out of range")
            return tuple(self.bonds[index])

        # Handle slice or array indexing
        if isinstance(index, slice):
            indices = list(range(*index.indices(self.n_bonds)))
        else:
            indices = list(index)

        # Create new BondList with selected bonds
        new_bond_list = BondList(len(indices))

        for new_idx, old_idx in enumerate(indices):
            new_bond_list.bonds[new_idx] = self.bonds[old_idx]
            new_bond_list.bond_order[new_idx] = self.bond_order[old_idx]
            new_bond_list.bond_type[new_idx] = self.bond_type[old_idx]

            # Copy optional properties
            for attr_name in [
                "detection_method",
                "confidence",
                "distance",
                "is_hydrogen_bond",
                "is_halogen_bond",
            ]:
                old_attr = getattr(self, attr_name)
                if old_attr is not None:
                    new_attr = new_bond_list._ensure_property(attr_name, old_attr.dtype)
                    new_attr[new_idx] = old_attr[old_idx]

            # Copy custom properties
            for prop_name in self._custom_properties:
                old_prop = getattr(self, prop_name)
                if not hasattr(new_bond_list, prop_name):
                    new_bond_list.add_property(prop_name, old_prop.dtype)
                new_prop = getattr(new_bond_list, prop_name)
                new_prop[new_idx] = old_prop[old_idx]

        new_bond_list._custom_properties = self._custom_properties.copy()
        return new_bond_list

    def __repr__(self) -> str:
        """String representation of BondList."""
        return f"BondList(n_bonds={self.n_bonds})"

    def __str__(self) -> str:
        """Detailed string representation."""
        lines = [f"BondList with {self.n_bonds} bonds"]

        if self.n_bonds > 0:
            # Show bond type distribution
            active_types = self.bond_type[: self.n_bonds]
            unique_types, counts = np.unique(active_types, return_counts=True)
            type_info = ", ".join(
                [
                    f"{count} {bond_type}"
                    for bond_type, count in zip(unique_types, counts)
                ]
            )
            lines.append(f"Bond types: {type_info}")

            # Show bond order distribution
            active_orders = self.bond_order[: self.n_bonds]
            unique_orders, counts = np.unique(active_orders, return_counts=True)
            order_info = ", ".join(
                [
                    f"{count} order-{order:.1f}"
                    for order, count in zip(unique_orders, counts)
                ]
            )
            lines.append(f"Bond orders: {order_info}")

            # Show optional properties status
            optional_count = sum(
                1
                for attr in [
                    "detection_method",
                    "confidence",
                    "distance",
                    "is_hydrogen_bond",
                    "is_halogen_bond",
                ]
                if getattr(self, attr) is not None
            )
            custom_count = len(self._custom_properties)
            lines.append(
                f"Properties: {optional_count} optional, {custom_count} custom"
            )

        return "\n".join(lines)


def both_atoms_in_set(bond: np.ndarray, atom_set: set) -> bool:
    """
    Check if both atoms in a bond are in the given set.

    Args:
        bond: Array with two atom indices
        atom_set: Set of atom indices

    Returns:
        True if both atoms are in the set
    """
    return bond[0] in atom_set and bond[1] in atom_set
