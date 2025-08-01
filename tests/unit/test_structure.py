"""
Unit tests for Structure class.

Tests cover core functionality based on acceptance criteria from requirements.md:
- AC-001 through AC-015 (Structure core functionality)
"""

import numpy as np
import pytest

from molr.core.structure import Structure


class TestStructureCreation:
    """Test Structure creation and initialization (AC-001 to AC-003)."""

    def test_structure_creation_valid_size(self):
        """Test creating structure with valid number of atoms."""
        structure = Structure(n_atoms=10)
        assert structure.n_atoms == 10
        assert len(structure) == 10
        assert structure.coord.shape == (10, 3)
        assert structure.atom_name.shape == (10,)
        assert structure.element.shape == (10,)
        assert structure.res_name.shape == (10,)
        assert structure.res_id.shape == (10,)
        assert structure.chain_id.shape == (10,)

    def test_structure_creation_zero_atoms_fails(self):
        """Test that creating structure with zero atoms fails."""
        with pytest.raises(ValueError, match="Number of atoms must be positive"):
            Structure(n_atoms=0)

    def test_structure_creation_negative_atoms_fails(self):
        """Test that creating structure with negative atoms fails."""
        with pytest.raises(ValueError, match="Number of atoms must be positive"):
            Structure(n_atoms=-5)

    def test_core_annotations_initialized(self):
        """Test that core annotations are properly initialized."""
        structure = Structure(n_atoms=5)

        # Check core annotations exist and have correct types
        assert structure.coord.dtype == np.float64
        assert structure.atom_name.dtype.kind == "U"
        assert structure.element.dtype.kind == "U"
        assert structure.res_name.dtype.kind == "U"
        assert structure.res_id.dtype == np.int32
        assert structure.chain_id.dtype.kind == "U"

        # Check optional annotations are None initially
        assert structure.alt_loc is None
        assert structure.occupancy is None
        assert structure.b_factor is None
        assert structure.charge is None
        assert structure.serial is None

    def test_custom_annotations_initialized_empty(self):
        """Test that custom annotations set is empty initially."""
        structure = Structure(n_atoms=3)
        assert len(structure._custom_annotations) == 0


class TestStructureAnnotations:
    """Test Structure annotation system (AC-004 to AC-006)."""

    def test_add_annotation_basic(self):
        """Test adding basic custom annotation."""
        structure = Structure(n_atoms=5)
        structure.add_annotation("custom_prop", dtype=np.float32, default_value=1.5)

        assert hasattr(structure, "custom_prop")
        assert "custom_prop" in structure._custom_annotations
        assert structure.custom_prop.shape == (5,)
        assert structure.custom_prop.dtype == np.float32
        assert np.all(structure.custom_prop == 1.5)

    def test_add_annotation_no_default(self):
        """Test adding annotation without default value."""
        structure = Structure(n_atoms=3)
        structure.add_annotation("test_prop", dtype=np.int32)

        assert hasattr(structure, "test_prop")
        assert structure.test_prop.shape == (3,)
        assert structure.test_prop.dtype == np.int32

    def test_add_annotation_duplicate_name_fails(self):
        """Test that adding annotation with existing name fails."""
        structure = Structure(n_atoms=3)
        structure.add_annotation("test_prop", dtype=np.float32)

        with pytest.raises(ValueError, match="Annotation 'test_prop' already exists"):
            structure.add_annotation("test_prop", dtype=np.float64)

    def test_add_annotation_conflicts_with_existing_attribute(self):
        """Test that adding annotation conflicts with existing attributes."""
        structure = Structure(n_atoms=3)

        with pytest.raises(ValueError, match="Annotation 'coord' already exists"):
            structure.add_annotation("coord", dtype=np.float32)

    def test_annotation_info(self):
        """Test getting annotation information."""
        structure = Structure(n_atoms=4)
        structure.add_annotation("custom1", dtype=np.float32)
        structure.add_annotation("custom2", dtype=np.int32, default_value=42)

        info = structure.get_annotation_info()

        # Check core annotations
        assert info["coord"]["type"] == "core"
        assert info["coord"]["initialized"] is True
        assert info["coord"]["dtype"] == np.float64

        # Check optional annotations
        assert info["alt_loc"]["type"] == "optional"
        assert info["alt_loc"]["initialized"] is False
        assert info["alt_loc"]["dtype"] is None

        # Check custom annotations
        assert info["custom1"]["type"] == "custom"
        assert info["custom1"]["initialized"] is True
        assert info["custom2"]["type"] == "custom"
        assert info["custom2"]["initialized"] is True


class TestStructureIndexing:
    """Test Structure indexing and slicing (AC-007 to AC-009)."""

    def test_integer_indexing(self):
        """Test single atom indexing."""
        structure = Structure(n_atoms=5)
        structure.coord = np.random.rand(5, 3)
        structure.atom_name[:] = ["CA", "CB", "CG", "CD", "CE"]
        structure.element[:] = ["C", "C", "C", "C", "C"]

        subset = structure[2]
        assert subset.n_atoms == 1
        assert np.array_equal(subset.coord, structure.coord[2:3])
        assert subset.atom_name[0] == "CG"

    def test_slice_indexing(self):
        """Test slice indexing."""
        structure = Structure(n_atoms=6)
        structure.coord = np.random.rand(6, 3)
        structure.atom_name[:] = ["CA", "CB", "CG", "CD", "CE", "CZ"]

        subset = structure[1:4]
        assert subset.n_atoms == 3
        assert np.array_equal(subset.coord, structure.coord[1:4])
        assert np.array_equal(subset.atom_name, ["CB", "CG", "CD"])

    def test_boolean_array_indexing(self):
        """Test boolean array indexing."""
        structure = Structure(n_atoms=5)
        structure.coord = np.random.rand(5, 3)
        structure.element[:] = ["C", "N", "C", "O", "C"]

        carbon_mask = structure.element == "C"
        carbon_subset = structure[carbon_mask]

        assert carbon_subset.n_atoms == 3
        assert np.all(carbon_subset.element == "C")

    def test_integer_array_indexing(self):
        """Test integer array indexing."""
        structure = Structure(n_atoms=6)
        structure.coord = np.random.rand(6, 3)
        structure.atom_name[:] = ["CA", "CB", "CG", "CD", "CE", "CZ"]

        indices = np.array([0, 2, 4])
        subset = structure[indices]

        assert subset.n_atoms == 3
        assert np.array_equal(subset.atom_name, ["CA", "CG", "CE"])

    def test_indexing_preserves_custom_annotations(self):
        """Test that indexing preserves custom annotations."""
        structure = Structure(n_atoms=4)
        structure.add_annotation("custom_prop", dtype=np.float32)
        structure.custom_prop[:] = [1.0, 2.0, 3.0, 4.0]

        subset = structure[1:3]
        assert hasattr(subset, "custom_prop")
        assert "custom_prop" in subset._custom_annotations
        assert np.array_equal(subset.custom_prop, [2.0, 3.0])

    def test_indexing_copies_optional_annotations(self):
        """Test that indexing copies optional annotations if they exist."""
        structure = Structure(n_atoms=4)
        structure._ensure_annotation("occupancy", np.float32, 1.0)
        structure.occupancy[:] = [0.8, 0.9, 1.0, 0.7]

        subset = structure[::2]  # Every other atom
        assert subset.occupancy is not None
        assert np.allclose(subset.occupancy, [0.8, 1.0])


class TestStructureOperations:
    """Test Structure operations and methods (AC-010 to AC-012)."""

    def test_copy_structure(self):
        """Test copying structure."""
        structure = Structure(n_atoms=3)
        structure.coord = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)
        structure.atom_name[:] = ["CA", "CB", "CG"]
        structure.add_annotation("custom", dtype=np.float32, default_value=42.0)

        copied = structure.copy()

        assert copied.n_atoms == structure.n_atoms
        assert np.array_equal(copied.coord, structure.coord)
        assert np.array_equal(copied.atom_name, structure.atom_name)
        assert hasattr(copied, "custom")
        assert np.array_equal(copied.custom, structure.custom)

        # Verify it's a deep copy
        copied.coord[0, 0] = 999
        assert structure.coord[0, 0] != 999

    def test_get_center_geometric(self):
        """Test geometric center calculation."""
        structure = Structure(n_atoms=3)
        structure.coord = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]], dtype=np.float64)

        center = structure.get_center()
        expected = np.array([2 / 3, 2 / 3, 0])
        assert np.allclose(center, expected)

    def test_get_center_weighted(self):
        """Test weighted center calculation."""
        structure = Structure(n_atoms=2)
        structure.coord = np.array([[0, 0, 0], [2, 0, 0]], dtype=np.float64)
        weights = np.array([1, 3])  # Second atom has 3x weight

        center = structure.get_center(weights=weights)
        expected = np.array([1.5, 0, 0])  # Weighted average
        assert np.allclose(center, expected)

    def test_translate(self):
        """Test structure translation."""
        structure = Structure(n_atoms=2)
        structure.coord = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64)

        translation = np.array([10, 20, 30])
        structure.translate(translation)

        expected = np.array([[11, 22, 33], [14, 25, 36]])
        assert np.allclose(structure.coord, expected)

    def test_translate_invalid_vector_fails(self):
        """Test that invalid translation vector fails."""
        structure = Structure(n_atoms=2)

        with pytest.raises(ValueError, match="Translation vector must be 3D"):
            structure.translate([1, 2])  # Only 2D

    def test_center_at_origin(self):
        """Test centering structure at origin."""
        structure = Structure(n_atoms=3)
        structure.coord = np.array([[1, 1, 1], [3, 3, 3], [5, 5, 5]], dtype=np.float64)

        structure.center_at_origin()

        # Center should now be at origin
        center = structure.get_center()
        assert np.allclose(center, [0, 0, 0], atol=1e-10)

    def test_get_masses(self):
        """Test atomic mass calculation."""
        structure = Structure(n_atoms=3)
        structure.element[:] = ["C", "N", "O"]

        masses = structure.get_masses()

        # Check that masses are reasonable (approximate values)
        assert masses[0] == pytest.approx(12.011, abs=0.1)  # Carbon
        assert masses[1] == pytest.approx(14.007, abs=0.1)  # Nitrogen
        assert masses[2] == pytest.approx(15.999, abs=0.1)  # Oxygen


class TestStructureClassification:
    """Test Structure classification flags (AC-013 to AC-015)."""

    def test_classification_flags_lazy_initialization(self):
        """Test that classification flags are computed on demand."""
        structure = Structure(n_atoms=3)

        # Flags should be None initially
        assert structure._is_backbone is None
        assert structure._is_sidechain is None
        assert structure._is_aromatic is None
        assert structure._is_ligand is None

        # Accessing property should trigger computation
        backbone = structure.is_backbone
        assert structure._is_backbone is not None
        assert isinstance(backbone, np.ndarray)
        assert backbone.dtype == bool

    def test_protein_backbone_classification(self):
        """Test protein backbone atom classification."""
        structure = Structure(n_atoms=4)
        structure.res_name[:] = "ALA"
        structure.atom_name[:] = ["N", "CA", "C", "O"]

        backbone = structure.is_backbone
        assert np.all(backbone)  # All should be backbone

        sidechain = structure.is_sidechain
        assert not np.any(sidechain)  # None should be sidechain

    def test_protein_sidechain_classification(self):
        """Test protein sidechain atom classification."""
        structure = Structure(n_atoms=2)
        structure.res_name[:] = "ALA"
        structure.atom_name[:] = ["CB", "HB1"]  # CB is sidechain, HB1 is hydrogen

        sidechain = structure.is_sidechain
        assert sidechain[0]  # CB should be sidechain
        assert not sidechain[1]  # HB1 is not in predefined sidechain atoms

    def test_nucleic_acid_classification(self):
        """Test nucleic acid backbone and base classification."""
        structure = Structure(n_atoms=3)
        structure.res_name[:] = ["DA", "DA", "DA"]
        structure.atom_name[:] = ["P", "C1'", "N1"]  # Phosphate, sugar, base

        backbone = structure.is_backbone
        sidechain = structure.is_sidechain

        assert backbone[0]  # P is backbone
        assert backbone[1]  # C1' is backbone
        assert not backbone[2]  # N1 is base

        assert not sidechain[0]  # P is not sidechain/base
        assert not sidechain[1]  # C1' is not sidechain/base
        assert sidechain[2]  # N1 is base (treated as "sidechain")

    def test_ligand_classification(self):
        """Test ligand atom classification."""
        structure = Structure(n_atoms=3)
        structure.res_name[:] = ["ATP", "ALA", "HOH"]  # Ligand, protein, water

        ligand = structure.is_ligand
        residue_type = structure.residue_type

        assert ligand[0]  # ATP is ligand
        assert not ligand[1]  # ALA is protein
        assert not ligand[2]  # HOH is water (not ligand)

        assert residue_type[0] == "LIGAND"
        assert residue_type[1] == "PROTEIN"
        assert residue_type[2] == "LIGAND"  # Water is classified as LIGAND by default

    def test_aromatic_classification(self):
        """Test aromatic atom classification."""
        structure = Structure(n_atoms=3)
        structure.res_name[:] = [
            "PHE",
            "ALA",
            "DA",
        ]  # Aromatic protein, non-aromatic protein, aromatic DNA

        aromatic = structure.is_aromatic

        assert aromatic[0]  # PHE is aromatic
        assert not aromatic[1]  # ALA is not aromatic
        assert aromatic[2]  # DA (adenine) is aromatic

    def test_residue_type_classification(self):
        """Test residue type classification."""
        structure = Structure(n_atoms=5)
        structure.res_name[:] = ["ALA", "DA", "A", "HOH", "HEM"]

        residue_type = structure.residue_type

        assert residue_type[0] == "PROTEIN"
        assert residue_type[1] == "DNA"
        assert residue_type[2] == "RNA"
        assert residue_type[3] == "LIGAND"  # Water classified as ligand
        assert residue_type[4] == "LIGAND"  # Heme is ligand


class TestStructureSelection:
    """Test Structure selection methods (placeholder for future implementation)."""

    def test_basic_element_selection(self):
        """Test basic element selection."""
        structure = Structure(n_atoms=4)
        structure.element[:] = ["C", "N", "C", "O"]

        carbon_mask = structure.select("element C")
        expected = np.array([True, False, True, False])
        assert np.array_equal(carbon_mask, expected)

    def test_basic_resname_selection(self):
        """Test basic residue name selection."""
        structure = Structure(n_atoms=3)
        structure.res_name[:] = ["ALA", "GLY", "ALA"]

        ala_mask = structure.select("resname ALA")
        expected = np.array([True, False, True])
        assert np.array_equal(ala_mask, expected)

    def test_basic_chain_selection(self):
        """Test basic chain selection."""
        structure = Structure(n_atoms=3)
        structure.chain_id[:] = ["A", "B", "A"]

        chain_a_mask = structure.select("chain A")
        expected = np.array([True, False, True])
        assert np.array_equal(chain_a_mask, expected)

    def test_backbone_selection(self):
        """Test backbone selection."""
        structure = Structure(n_atoms=2)
        structure.res_name[:] = "ALA"
        structure.atom_name[:] = ["CA", "CB"]

        backbone_mask = structure.select("backbone")
        expected = np.array([True, False])  # CA is backbone, CB is not
        assert np.array_equal(backbone_mask, expected)

    def test_ligand_selection(self):
        """Test ligand selection."""
        structure = Structure(n_atoms=2)
        structure.res_name[:] = ["ALA", "HEM"]

        ligand_mask = structure.select("ligand")
        expected = np.array([False, True])  # ALA is protein, HEM is ligand
        assert np.array_equal(ligand_mask, expected)

    def test_unsupported_selection_fails(self):
        """Test that unsupported selections fail appropriately."""
        structure = Structure(n_atoms=2)

        with pytest.raises(NotImplementedError):
            structure.select("within 5 of resname ALA")


class TestStructureStringRepresentation:
    """Test Structure string representation methods."""

    def test_repr(self):
        """Test Structure __repr__ method."""
        structure = Structure(n_atoms=42)
        repr_str = repr(structure)
        assert "Structure(n_atoms=42)" == repr_str

    def test_str_basic(self):
        """Test Structure __str__ method with basic structure."""
        structure = Structure(n_atoms=5)
        str_repr = str(structure)

        assert "Structure with 5 atoms" in str_repr
        assert "Annotations:" in str_repr

    def test_str_with_multiple_chains(self):
        """Test Structure __str__ method with multiple chains."""
        structure = Structure(n_atoms=4)
        structure.chain_id[:] = ["A", "A", "B", "B"]

        str_repr = str(structure)
        assert "2 in chain A" in str_repr
        assert "2 in chain B" in str_repr

    def test_str_with_residue_types(self):
        """Test Structure __str__ method showing residue types."""
        structure = Structure(n_atoms=3)
        structure.res_name[:] = ["ALA", "DA", "HEM"]

        str_repr = str(structure)
        # Force classification by accessing property
        _ = structure.residue_type
        str_repr = str(structure)

        assert "PROTEIN" in str_repr
        assert "DNA" in str_repr
        assert "LIGAND" in str_repr


class TestStructureBondOperations:
    """Test Structure bond-related operations (placeholder for integration with BondList)."""

    def test_get_bonds_to_basic(self):
        """Test basic bond detection by distance."""
        structure = Structure(n_atoms=3)
        structure.coord = np.array(
            [[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0]], dtype=np.float64
        )

        # Find atoms potentially bonded to atom 0
        bonds = structure.get_bonds_to([1], max_distance=2.0)

        # Atom 0 should be bonded to atom 1 (distance 1.5 < 2.0)
        # Atom 1 should be bonded to atom 0 (same bond)
        # Atom 2 should not be bonded to atom 1 (distance 1.5 < 2.0) - wait this should be bonded
        # Let's check atom 0 bonds
        bonds_to_1 = structure.get_bonds_to([1], max_distance=2.0)
        assert bonds_to_1[0]  # Atom 0 is within 2.0 of atom 1
        assert not bonds_to_1[1]  # Atom 1 doesn't bond to itself
        assert bonds_to_1[2]  # Atom 2 is within 2.0 of atom 1 (distance = 1.5)

    def test_get_bonds_to_no_self_bonds(self):
        """Test that atoms don't bond to themselves."""
        structure = Structure(n_atoms=2)
        structure.coord = np.array([[0, 0, 0], [0.5, 0, 0]], dtype=np.float64)

        bonds = structure.get_bonds_to([0], max_distance=2.0)
        assert not bonds[0]  # Atom 0 should not bond to itself
        assert bonds[1]  # Atom 1 should bond to atom 0
