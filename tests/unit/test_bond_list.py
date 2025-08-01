"""
Unit tests for BondList class.

Tests cover core functionality based on acceptance criteria from requirements.md:
- AC-016 through AC-030 (BondList functionality)
"""

import numpy as np
import pytest

from molr.core.bond_list import BondList, both_atoms_in_set


class TestBondListCreation:
    """Test BondList creation and initialization (AC-016 to AC-018)."""

    def test_bond_list_creation_empty(self):
        """Test creating empty BondList."""
        bond_list = BondList()
        assert bond_list.n_bonds == 0
        assert len(bond_list) == 0
        assert bond_list._capacity >= 10  # Minimum capacity

    def test_bond_list_creation_with_capacity(self):
        """Test creating BondList with initial capacity."""
        bond_list = BondList(n_bonds=5)
        assert bond_list.n_bonds == 5
        assert len(bond_list) == 5
        assert bond_list._capacity >= 5

    def test_core_bond_arrays_initialized(self):
        """Test that core bond arrays are properly initialized."""
        bond_list = BondList(n_bonds=3)

        # Check core arrays exist and have correct types
        assert bond_list.bonds.shape[1] == 2  # Pairs of atoms
        assert bond_list.bonds.dtype == np.int32
        assert bond_list.bond_order.dtype == np.float32
        assert bond_list.bond_type.dtype.kind == "U"

        # Check optional arrays are None initially
        assert bond_list.detection_method is None
        assert bond_list.confidence is None
        assert bond_list.distance is None
        assert bond_list.is_hydrogen_bond is None
        assert bond_list.is_halogen_bond is None

    def test_custom_properties_initialized_empty(self):
        """Test that custom properties set is empty initially."""
        bond_list = BondList()
        assert len(bond_list._custom_properties) == 0


class TestBondListProperties:
    """Test BondList property management (AC-019 to AC-021)."""

    def test_add_property_basic(self):
        """Test adding basic custom property."""
        bond_list = BondList(n_bonds=3)
        bond_list.add_property("custom_energy", dtype=np.float32, default_value=42.0)

        assert hasattr(bond_list, "custom_energy")
        assert "custom_energy" in bond_list._custom_properties
        assert bond_list.custom_energy.shape[0] >= 3
        assert bond_list.custom_energy.dtype == np.float32
        assert np.all(bond_list.custom_energy[:3] == 42.0)

    def test_add_property_no_default(self):
        """Test adding property without default value."""
        bond_list = BondList(n_bonds=2)
        bond_list.add_property("test_prop", dtype=np.int32)

        assert hasattr(bond_list, "test_prop")
        assert bond_list.test_prop.dtype == np.int32

    def test_add_property_duplicate_name_fails(self):
        """Test that adding property with existing name fails."""
        bond_list = BondList()
        bond_list.add_property("test_prop", dtype=np.float32)

        with pytest.raises(ValueError, match="Property 'test_prop' already exists"):
            bond_list.add_property("test_prop", dtype=np.float64)

    def test_add_property_conflicts_with_existing_attribute(self):
        """Test that adding property conflicts with existing attributes."""
        bond_list = BondList()

        with pytest.raises(ValueError, match="Property 'bonds' already exists"):
            bond_list.add_property("bonds", dtype=np.float32)

    def test_ensure_property_lazy_initialization(self):
        """Test lazy initialization of optional properties."""
        bond_list = BondList(n_bonds=2)

        # Initially None
        assert bond_list.confidence is None

        # Ensure it exists
        confidence = bond_list._ensure_property("confidence", np.float32, 0.5)

        assert bond_list.confidence is not None
        assert bond_list.confidence is confidence
        assert bond_list.confidence.dtype == np.float32
        assert np.all(bond_list.confidence[:2] == 0.5)


class TestBondListAddBonds:
    """Test adding bonds to BondList (AC-022 to AC-024)."""

    def test_add_single_bond_basic(self):
        """Test adding single bond with basic parameters."""
        bond_list = BondList()

        bond_idx = bond_list.add_bond(0, 1, bond_order=1.0, bond_type="covalent")

        assert bond_idx == 0
        assert bond_list.n_bonds == 1
        assert np.array_equal(bond_list.bonds[0], [0, 1])
        assert bond_list.bond_order[0] == 1.0
        assert bond_list.bond_type[0] == "covalent"

    def test_add_single_bond_atom_ordering(self):
        """Test that atoms are ordered consistently (smaller index first)."""
        bond_list = BondList()

        bond_idx = bond_list.add_bond(3, 1)  # Reversed order

        assert np.array_equal(bond_list.bonds[0], [1, 3])  # Should be sorted

    def test_add_single_bond_with_properties(self):
        """Test adding bond with optional properties."""
        bond_list = BondList()
        bond_list._ensure_property("confidence", np.float32)
        bond_list._ensure_property("distance", np.float32)

        bond_idx = bond_list.add_bond(0, 1, confidence=0.95, distance=1.5)

        assert bond_list.confidence[0] == 0.95
        assert bond_list.distance[0] == 1.5

    def test_add_single_bond_same_atom_fails(self):
        """Test that adding bond between same atom fails."""
        bond_list = BondList()

        with pytest.raises(ValueError, match="Cannot create bond between same atom"):
            bond_list.add_bond(2, 2)

    def test_add_single_bond_negative_indices_fail(self):
        """Test that negative atom indices fail."""
        bond_list = BondList()

        with pytest.raises(ValueError, match="Atom indices must be non-negative"):
            bond_list.add_bond(-1, 2)

        with pytest.raises(ValueError, match="Atom indices must be non-negative"):
            bond_list.add_bond(1, -2)

    def test_add_multiple_bonds_basic(self):
        """Test adding multiple bonds at once."""
        bond_list = BondList()

        bond_pairs = [(0, 1), (2, 3), (1, 4)]
        bond_orders = [1.0, 2.0, 1.5]
        bond_types = ["single", "double", "aromatic"]

        indices = bond_list.add_bonds(bond_pairs, bond_orders, bond_types)

        assert len(indices) == 3
        assert bond_list.n_bonds == 3
        assert np.array_equal(indices, [0, 1, 2])

        # Check bonds are stored correctly
        assert np.array_equal(bond_list.bonds[:3], [[0, 1], [2, 3], [1, 4]])
        assert np.array_equal(bond_list.bond_order[:3], [1.0, 2.0, 1.5])
        assert np.array_equal(bond_list.bond_type[:3], ["single", "double", "aromatic"])

    def test_add_multiple_bonds_default_parameters(self):
        """Test adding multiple bonds with default parameters."""
        bond_list = BondList()

        bond_pairs = [(0, 1), (2, 3)]
        indices = bond_list.add_bonds(bond_pairs)

        assert bond_list.n_bonds == 2
        assert np.all(bond_list.bond_order[:2] == 1.0)
        assert np.all(bond_list.bond_type[:2] == "covalent")

    def test_add_multiple_bonds_with_properties(self):
        """Test adding multiple bonds with custom properties."""
        bond_list = BondList()
        bond_list._ensure_property("confidence", np.float32)

        bond_pairs = [(0, 1), (2, 3)]
        confidences = [0.9, 0.8]

        bond_list.add_bonds(bond_pairs, confidence=confidences)

        assert np.allclose(bond_list.confidence[:2], [0.9, 0.8])

    def test_add_multiple_bonds_mismatched_lengths_fail(self):
        """Test that mismatched list lengths fail."""
        bond_list = BondList()

        bond_pairs = [(0, 1), (2, 3)]
        bond_orders = [1.0]  # Wrong length

        with pytest.raises(ValueError, match="bond_orders length must match"):
            bond_list.add_bonds(bond_pairs, bond_orders)

    def test_add_bonds_empty_list(self):
        """Test adding empty list of bonds."""
        bond_list = BondList()

        indices = bond_list.add_bonds([])

        assert len(indices) == 0
        assert bond_list.n_bonds == 0

    def test_capacity_growth(self):
        """Test that capacity grows automatically when needed."""
        bond_list = BondList()
        initial_capacity = bond_list._capacity

        # Add enough bonds to exceed initial capacity
        bond_pairs = [(i, i + 1) for i in range(initial_capacity + 5)]
        bond_list.add_bonds(bond_pairs)

        assert bond_list.n_bonds == initial_capacity + 5
        assert bond_list._capacity > initial_capacity


class TestBondListRemoveBonds:
    """Test removing bonds from BondList (AC-025)."""

    def test_remove_single_bond(self):
        """Test removing single bond by index."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5)])

        bond_list.remove_bonds(1)  # Remove middle bond

        assert bond_list.n_bonds == 2
        assert np.array_equal(bond_list.bonds[:2], [[0, 1], [4, 5]])

    def test_remove_multiple_bonds(self):
        """Test removing multiple bonds by indices."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5), (6, 7)])

        bond_list.remove_bonds([0, 2])  # Remove first and third bonds

        assert bond_list.n_bonds == 2
        assert np.array_equal(bond_list.bonds[:2], [[2, 3], [6, 7]])

    def test_remove_bonds_with_properties(self):
        """Test that removing bonds also removes properties."""
        bond_list = BondList()
        bond_list.add_property("energy", np.float32)
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5)], energy=[1.0, 2.0, 3.0])

        bond_list.remove_bonds(1)  # Remove middle bond

        assert bond_list.n_bonds == 2
        assert np.array_equal(bond_list.energy[:2], [1.0, 3.0])

    def test_remove_bonds_invalid_indices(self):
        """Test removing bonds with invalid indices."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        # Should not crash with invalid indices
        bond_list.remove_bonds([1, 5, -1])  # Only index 1 is valid

        assert bond_list.n_bonds == 1
        assert np.array_equal(bond_list.bonds[0], [0, 1])

    def test_remove_bonds_empty_list(self):
        """Test removing empty list of bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        bond_list.remove_bonds([])

        assert bond_list.n_bonds == 2  # Should be unchanged


class TestBondListQuerying:
    """Test querying bonds in BondList (AC-026 to AC-028)."""

    def test_get_bonds_for_atom(self):
        """Test getting all bonds for a specific atom."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2), (2, 3), (1, 4)])

        bonds_for_1 = bond_list.get_bonds_for_atom(1)

        # Atom 1 is in bonds 0, 1, and 3
        expected = np.array([0, 1, 3])
        assert np.array_equal(bonds_for_1, expected)

    def test_get_bonds_for_atom_no_bonds(self):
        """Test getting bonds for atom with no bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        bonds_for_5 = bond_list.get_bonds_for_atom(5)

        assert len(bonds_for_5) == 0

    def test_get_bonds_for_atom_empty_bond_list(self):
        """Test getting bonds for atom in empty bond list."""
        bond_list = BondList()

        bonds_for_0 = bond_list.get_bonds_for_atom(0)

        assert len(bonds_for_0) == 0

    def test_get_neighbors(self):
        """Test getting neighbor atoms for a specific atom."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2), (1, 3), (2, 4)])

        neighbors_of_1 = bond_list.get_neighbors(1)

        # Atom 1 is connected to atoms 0, 2, and 3
        expected = np.array([0, 2, 3])
        assert np.array_equal(np.sort(neighbors_of_1), np.sort(expected))

    def test_get_neighbors_no_bonds(self):
        """Test getting neighbors for atom with no bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        neighbors_of_5 = bond_list.get_neighbors(5)

        assert len(neighbors_of_5) == 0

    def test_get_bond_matrix(self):
        """Test creating bond adjacency matrix."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2), (0, 2)])  # Triangle

        matrix = bond_list.get_bond_matrix(n_atoms=4)

        assert matrix.shape == (4, 4)
        assert matrix.dtype == bool

        # Check symmetry
        assert matrix[0, 1] and matrix[1, 0]  # Bond 0-1
        assert matrix[1, 2] and matrix[2, 1]  # Bond 1-2
        assert matrix[0, 2] and matrix[2, 0]  # Bond 0-2
        assert not matrix[3, 0]  # No bond to atom 3

        # Diagonal should be False
        assert not np.any(np.diag(matrix))

    def test_get_bond_matrix_invalid_bonds(self):
        """Test bond matrix with bonds referencing invalid atoms."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 5)])  # Bond 2-5 exceeds n_atoms=4

        matrix = bond_list.get_bond_matrix(n_atoms=4)

        assert matrix[0, 1] and matrix[1, 0]  # Valid bond
        assert not matrix[2, 3]  # Invalid bond should be ignored

    def test_validate_bonds_all_valid(self):
        """Test bond validation with all valid bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2), (2, 3)])

        is_valid, invalid_indices = bond_list.validate_bonds(n_atoms=4)

        assert is_valid
        assert len(invalid_indices) == 0

    def test_validate_bonds_some_invalid(self):
        """Test bond validation with some invalid bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 5), (1, 3)])  # Bond 1 is invalid

        is_valid, invalid_indices = bond_list.validate_bonds(n_atoms=4)

        assert not is_valid
        assert invalid_indices == [1]  # Bond at index 1 is invalid

    def test_validate_bonds_empty(self):
        """Test bond validation with empty bond list."""
        bond_list = BondList()

        is_valid, invalid_indices = bond_list.validate_bonds(n_atoms=10)

        assert is_valid
        assert len(invalid_indices) == 0


class TestBondListFiltering:
    """Test filtering bonds by atoms (AC-029 to AC-030)."""

    def test_filter_by_atoms_basic(self):
        """Test basic filtering by atom selection."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2), (2, 3), (3, 4)])

        # Keep only atoms 0, 1, 2
        filtered = bond_list.filter_by_atoms([0, 1, 2])

        assert filtered.n_bonds == 2  # Bonds (0,1) and (1,2)
        # Indices should be remapped: 0->0, 1->1, 2->2
        assert np.array_equal(filtered.bonds[:2], [[0, 1], [1, 2]])

    def test_filter_by_atoms_with_remapping(self):
        """Test filtering with index remapping."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 2), (2, 4), (1, 3)])

        # Keep only atoms 0, 2, 4 (non-consecutive)
        filtered = bond_list.filter_by_atoms([0, 2, 4])

        assert filtered.n_bonds == 2  # Bonds (0,2) and (2,4)
        # Remapping: 0->0, 2->1, 4->2
        assert np.array_equal(filtered.bonds[:2], [[0, 1], [1, 2]])

    def test_filter_by_atoms_preserves_properties(self):
        """Test that filtering preserves bond properties."""
        bond_list = BondList()
        bond_list.add_property("energy", np.float32)
        bond_list.add_bonds(
            [(0, 1), (1, 2), (2, 3)],
            bond_orders=[1.0, 2.0, 1.5],
            energy=[10.0, 20.0, 15.0],
        )

        filtered = bond_list.filter_by_atoms([0, 1, 2])

        assert filtered.n_bonds == 2
        assert np.array_equal(filtered.bond_order[:2], [1.0, 2.0])
        assert np.array_equal(filtered.energy[:2], [10.0, 20.0])
        assert "energy" in filtered._custom_properties

    def test_filter_by_atoms_no_matching_bonds(self):
        """Test filtering with no matching bonds."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        # Keep atoms that don't form any bonds together
        filtered = bond_list.filter_by_atoms([0, 2])

        assert filtered.n_bonds == 0

    def test_filter_by_atoms_single_atom(self):
        """Test filtering with single atom (no bonds possible)."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2)])

        filtered = bond_list.filter_by_atoms([1])

        assert filtered.n_bonds == 0  # Can't have bonds with only one atom

    def test_filter_by_atoms_empty_selection(self):
        """Test filtering with empty atom selection."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (1, 2)])

        filtered = bond_list.filter_by_atoms([])

        assert filtered.n_bonds == 0


class TestBondListIndexing:
    """Test BondList indexing operations."""

    def test_integer_indexing(self):
        """Test single bond indexing."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5)])

        bond = bond_list[1]
        assert bond == (2, 3)

    def test_integer_indexing_out_of_range(self):
        """Test indexing with out of range index."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1)])

        with pytest.raises(IndexError, match="Bond index out of range"):
            bond_list[5]

    def test_slice_indexing(self):
        """Test slice indexing."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5), (6, 7)])
        bond_list.add_property("energy", np.float32)
        bond_list.energy[:4] = [1.0, 2.0, 3.0, 4.0]

        subset = bond_list[1:3]

        assert subset.n_bonds == 2
        assert np.array_equal(subset.bonds[:2], [[2, 3], [4, 5]])
        assert np.array_equal(subset.energy[:2], [2.0, 3.0])

    def test_array_indexing(self):
        """Test array indexing."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3), (4, 5), (6, 7)])

        indices = [0, 2]
        subset = bond_list[indices]

        assert subset.n_bonds == 2
        assert np.array_equal(subset.bonds[:2], [[0, 1], [4, 5]])


class TestBondListStringRepresentation:
    """Test BondList string representation methods."""

    def test_repr(self):
        """Test BondList __repr__ method."""
        bond_list = BondList()
        bond_list.add_bonds([(0, 1), (2, 3)])

        repr_str = repr(bond_list)
        assert "BondList(n_bonds=2)" == repr_str

    def test_str_empty(self):
        """Test BondList __str__ method with empty list."""
        bond_list = BondList()

        str_repr = str(bond_list)
        assert "BondList with 0 bonds" in str_repr

    def test_str_with_bonds(self):
        """Test BondList __str__ method with bonds."""
        bond_list = BondList()
        bond_list.add_bonds(
            [(0, 1), (2, 3)], bond_orders=[1.0, 2.0], bond_types=["single", "double"]
        )

        str_repr = str(bond_list)
        assert "BondList with 2 bonds" in str_repr
        assert "single" in str_repr
        assert "double" in str_repr
        assert "order-1.0" in str_repr
        assert "order-2.0" in str_repr


class TestBothAtomsInSetHelper:
    """Test both_atoms_in_set helper function."""

    def test_both_atoms_in_set_true(self):
        """Test when both atoms are in set."""
        bond = np.array([1, 3])
        atom_set = {0, 1, 2, 3, 4}

        assert both_atoms_in_set(bond, atom_set)

    def test_both_atoms_in_set_false(self):
        """Test when one atom is not in set."""
        bond = np.array([1, 5])
        atom_set = {0, 1, 2, 3, 4}

        assert not both_atoms_in_set(bond, atom_set)

    def test_both_atoms_in_set_neither(self):
        """Test when neither atom is in set."""
        bond = np.array([7, 8])
        atom_set = {0, 1, 2, 3, 4}

        assert not both_atoms_in_set(bond, atom_set)
