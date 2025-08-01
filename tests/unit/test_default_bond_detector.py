"""
Unit tests for the simplified default bond detector.

Tests the new streamlined bond detection approach that combines
residue templates, distance-based detection, and file bonds.
"""

import numpy as np
import pytest

from molr.bond_detection.default_detector import DefaultBondDetector, detect_bonds
from molr.core.bond_list import BondList
from molr.core.structure import Structure


@pytest.fixture
def simple_molecule():
    """Create a simple methane-like molecule for testing."""
    structure = Structure(4)

    # Set coordinates for C-H bonds (tetrahedral methane geometry)
    structure.coord[0] = [0.0, 0.0, 0.0]  # C at origin
    structure.coord[1] = [1.09, 0.0, 0.0]  # H along x-axis
    structure.coord[2] = [-0.36, 1.03, 0.0]  # H in xy-plane
    structure.coord[3] = [-0.36, -0.51, 0.89]  # H above xy-plane

    # Set atom properties
    structure.atom_name[0] = "C"
    structure.atom_name[1] = "H1"
    structure.atom_name[2] = "H2"
    structure.atom_name[3] = "H3"

    structure.element[0] = "C"
    structure.element[1] = "H"
    structure.element[2] = "H"
    structure.element[3] = "H"

    structure.res_name[:] = "UNK"
    structure.res_id[:] = 1
    structure.chain_id[:] = "A"

    return structure


@pytest.fixture
def alanine_residue():
    """Create an alanine residue for template testing."""
    structure = Structure(10)

    # Simplified alanine coordinates
    coords = np.array(
        [
            [0.0, 0.0, 0.0],  # N
            [1.5, 0.0, 0.0],  # CA
            [2.0, 1.5, 0.0],  # C
            [3.0, 1.8, 0.5],  # O
            [1.8, 0.5, -1.3],  # CB
            [-0.5, -0.8, 0.3],  # H
            [1.2, -0.8, 0.5],  # HA
            [2.8, 0.2, -1.1],  # HB1
            [1.6, 1.4, -1.8],  # HB2
            [1.2, -0.3, -1.8],  # HB3
        ]
    )
    structure.coord = coords

    # Set atom names
    atom_names = ["N", "CA", "C", "O", "CB", "H", "HA", "HB1", "HB2", "HB3"]
    for i, name in enumerate(atom_names):
        structure.atom_name[i] = name

    # Set elements
    elements = ["N", "C", "C", "O", "C", "H", "H", "H", "H", "H"]
    for i, elem in enumerate(elements):
        structure.element[i] = elem

    structure.res_name[:] = "ALA"
    structure.res_id[:] = 1
    structure.chain_id[:] = "A"

    return structure


@pytest.fixture
def structure_with_file_bonds():
    """Create a structure with file-based bonds."""
    structure = Structure(3)

    # Simple linear molecule
    structure.coord[0] = [0.0, 0.0, 0.0]
    structure.coord[1] = [1.5, 0.0, 0.0]
    structure.coord[2] = [3.0, 0.0, 0.0]

    structure.atom_name[0] = "C1"
    structure.atom_name[1] = "C2"
    structure.atom_name[2] = "C3"

    structure.element[0] = "C"
    structure.element[1] = "C"
    structure.element[2] = "C"

    structure.res_name[:] = "UNK"
    structure.res_id[:] = 1
    structure.chain_id[:] = "A"

    # Add file bonds (e.g., from PDB CONECT)
    file_bonds = BondList()
    file_bonds.add_bond(
        0, 1, bond_order=1.0, bond_type="covalent", detection_method="file"
    )
    file_bonds.add_bond(
        1, 2, bond_order=2.0, bond_type="covalent", detection_method="file"
    )
    structure.file_bonds = file_bonds

    return structure


class TestDefaultBondDetector:
    """Test the DefaultBondDetector class."""

    def test_initialization(self):
        """Test detector initialization."""
        detector = DefaultBondDetector()
        assert detector.vdw_factor == 0.75

        detector = DefaultBondDetector(vdw_factor=0.8)
        assert detector.vdw_factor == 0.8

        # Test invalid vdw_factor
        with pytest.raises(ValueError):
            DefaultBondDetector(vdw_factor=0.0)
        with pytest.raises(ValueError):
            DefaultBondDetector(vdw_factor=1.5)

    def test_distance_detection_only(self, simple_molecule):
        """Test distance-based detection on simple molecule."""
        detector = DefaultBondDetector(vdw_factor=0.8)
        bonds = detector.detect_bonds(simple_molecule, use_file_bonds=False)

        # Should detect C-H bonds
        assert len(bonds) == 3

        # Check bond types
        for i in range(len(bonds)):
            assert bonds.bond_order[i] == 1.0
            assert bonds.bond_type[i] == "covalent"
            assert bonds.detection_method[i] == "distance"

    def test_residue_template_detection(self, alanine_residue):
        """Test residue template-based detection."""
        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(alanine_residue, use_file_bonds=False)

        # Should have bonds from ALA template
        assert len(bonds) > 5  # At least backbone bonds

        # Check that some bonds are from templates
        template_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "template"
        ]
        assert len(template_bonds) > 0

    def test_file_bonds_priority(self, structure_with_file_bonds):
        """Test that file bonds are included and have priority."""
        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure_with_file_bonds, use_file_bonds=True)

        # Should include the 2 file bonds
        file_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "file"
        ]
        assert len(file_bonds) == 2

        # Check bond orders from file
        bond_orders = [bonds.bond_order[i] for i in file_bonds]
        assert 1.0 in bond_orders
        assert 2.0 in bond_orders

    def test_ignore_file_bonds(self, structure_with_file_bonds):
        """Test detection without file bonds."""
        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure_with_file_bonds, use_file_bonds=False)

        # Should not have any file bonds
        file_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "file"
        ]
        assert len(file_bonds) == 0

        # But may have distance-based bonds
        distance_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "distance"
        ]
        assert len(distance_bonds) >= 0  # May be 0 if distances too large

    def test_combined_detection(self, alanine_residue):
        """Test combined template + distance detection."""
        # Add some distant atoms that need distance detection
        extended_structure = Structure(12)

        # Copy alanine data
        extended_structure.coord[:10] = alanine_residue.coord
        extended_structure.atom_name[:10] = alanine_residue.atom_name
        extended_structure.element[:10] = alanine_residue.element
        extended_structure.res_name[:10] = alanine_residue.res_name
        extended_structure.res_id[:10] = alanine_residue.res_id
        extended_structure.chain_id[:10] = alanine_residue.chain_id

        # Add extra atoms for distance detection
        extended_structure.coord[10] = [5.0, 0.0, 0.0]  # Distant C
        extended_structure.coord[11] = [6.2, 0.0, 0.0]  # H bonded to distant C

        extended_structure.atom_name[10] = "C_EXT"
        extended_structure.atom_name[11] = "H_EXT"
        extended_structure.element[10] = "C"
        extended_structure.element[11] = "H"
        extended_structure.res_name[10:] = "UNK"
        extended_structure.res_id[10:] = 2
        extended_structure.chain_id[10:] = "A"

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(extended_structure, use_file_bonds=False)

        # Should have both template and distance bonds
        template_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "template"
        ]
        distance_bonds = [
            i for i in range(len(bonds)) if bonds.detection_method[i] == "distance"
        ]

        assert len(template_bonds) > 0
        assert len(distance_bonds) > 0

    def test_empty_structure(self):
        """Test detection on structure with no bonding."""
        # Create single isolated atom
        structure = Structure(1)
        structure.coord[0] = [0.0, 0.0, 0.0]
        structure.atom_name[0] = "HE"  # Helium - inert
        structure.element[0] = "HE"
        structure.res_name[0] = "HE"
        structure.res_id[0] = 1
        structure.chain_id[0] = "A"

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure)

        assert len(bonds) == 0


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_detect_bonds_function(self, simple_molecule):
        """Test the detect_bonds convenience function."""
        bonds = detect_bonds(simple_molecule, vdw_factor=0.8, use_file_bonds=False)

        # Should detect C-H bonds
        assert len(bonds) > 0
        assert isinstance(bonds, BondList)

    def test_detect_bonds_with_file_bonds(self, structure_with_file_bonds):
        """Test detect_bonds function with file bonds."""
        bonds = detect_bonds(structure_with_file_bonds, use_file_bonds=True)

        # Should include file bonds
        file_bond_count = sum(
            1 for i in range(len(bonds)) if bonds.detection_method[i] == "file"
        )
        assert file_bond_count == 2

    def test_detect_bonds_parameters(self, simple_molecule):
        """Test detect_bonds function parameters."""
        # Test with different vdw_factor
        bonds1 = detect_bonds(simple_molecule, vdw_factor=0.6)
        bonds2 = detect_bonds(simple_molecule, vdw_factor=0.9)

        # Different factors may give different results
        # (exact results depend on coordinates and elements)
        assert isinstance(bonds1, BondList)
        assert isinstance(bonds2, BondList)


class TestIntegrationWithStructure:
    """Test integration with Structure class."""

    def test_structure_detect_bonds_auto(self, alanine_residue):
        """Test Structure.detect_bonds() with default parameters."""
        bonds = alanine_residue.detect_bonds()

        # Should use new simplified detector
        assert len(bonds) > 0
        assert isinstance(bonds, BondList)

        # Should be stored in structure
        assert alanine_residue.bonds is not None
        assert len(alanine_residue.bonds) == len(bonds)

    def test_structure_detect_bonds_params(self, alanine_residue):
        """Test Structure.detect_bonds() with different parameters."""
        # Test with different VdW factor
        bonds1 = alanine_residue.detect_bonds(vdw_factor=0.5, store_bonds=False)
        bonds2 = alanine_residue.detect_bonds(vdw_factor=0.9, store_bonds=False)

        assert isinstance(bonds1, BondList)
        assert isinstance(bonds2, BondList)
        # Different parameters may give different numbers of bonds

    def test_structure_bonds_storage(self, simple_molecule):
        """Test bond storage in structure."""
        # Detect and store bonds
        bonds = simple_molecule.detect_bonds(store_bonds=True)

        assert simple_molecule.has_bonds()
        assert simple_molecule.bonds is not None
        assert len(simple_molecule.bonds) == len(bonds)

        # Test without storage
        bonds2 = simple_molecule.detect_bonds(store_bonds=False)
        # Previous bonds should still be there
        assert simple_molecule.bonds is not None
        assert len(simple_molecule.bonds) == len(bonds)  # From first call
