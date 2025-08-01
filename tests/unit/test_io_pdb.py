"""
Unit tests for the PDB I/O module.

Tests cover PDB parsing and writing functionality for the space module.
"""

import os
import tempfile

import numpy as np
import pytest

from molr.core.bond_list import BondList
from molr.core.structure import Structure
from molr.core.structure_ensemble import StructureEnsemble
from molr.io.pdb import PDBParser


@pytest.fixture
def sample_pdb_content():
    """Simple PDB content for testing."""
    return """ATOM      1  N   ALA A   1      20.154  16.967  22.478  1.00 20.00           N  
ATOM      2  CA  ALA A   1      21.618  17.209  22.700  1.00 20.00           C  
ATOM      3  C   ALA A   1      22.096  16.445  23.928  1.00 20.00           C  
ATOM      4  O   ALA A   1      22.031  15.218  23.963  1.00 20.00           O  
ATOM      5  CB  ALA A   1      22.430  16.794  21.490  1.00 20.00           C  
HETATM    6  O   HOH A   2      17.245  16.894  21.730  1.00 30.00           O  
CONECT    1    2
CONECT    2    1    3    5
CONECT    3    2    4
END"""


@pytest.fixture
def sample_pdb_with_models():
    """PDB content with multiple models for trajectory testing."""
    return """MODEL        1
ATOM      1  N   ALA A   1      20.154  16.967  22.478  1.00 20.00           N  
ATOM      2  CA  ALA A   1      21.618  17.209  22.700  1.00 20.00           C  
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1      20.200  17.000  22.500  1.00 20.00           N  
ATOM      2  CA  ALA A   1      21.650  17.240  22.720  1.00 20.00           C  
ENDMDL
END"""


class TestPDBParser:
    """Test PDB parsing functionality."""

    def test_parse_string_basic(self, sample_pdb_content):
        """Test basic PDB string parsing."""
        parser = PDBParser()
        structure = parser.parse_string(sample_pdb_content)

        # Check basic structure properties
        assert structure.n_atoms == 6
        assert len(structure.atom_name) == 6
        assert len(structure.coord) == 6

        # Check atom names
        expected_names = ["N", "CA", "C", "O", "CB", "O"]
        assert list(structure.atom_name) == expected_names

        # Check elements
        expected_elements = ["N", "C", "C", "O", "C", "O"]
        assert list(structure.element) == expected_elements

        # Check residue names
        expected_res_names = ["ALA", "ALA", "ALA", "ALA", "ALA", "HOH"]
        assert list(structure.res_name) == expected_res_names

        # Check coordinates for first atom
        np.testing.assert_array_almost_equal(
            structure.coord[0], [20.154, 16.967, 22.478], decimal=3
        )

    def test_parse_file(self, sample_pdb_content):
        """Test PDB file parsing."""
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(sample_pdb_content)
            temp_filename = f.name

        try:
            parser = PDBParser()
            structure = parser.parse_file(temp_filename)

            assert structure.n_atoms == 6
            assert structure.atom_name[0] == "N"
            assert structure.element[0] == "N"

        finally:
            os.unlink(temp_filename)

    def test_conect_records_parsing(self, sample_pdb_content):
        """Test CONECT record parsing."""
        parser = PDBParser()
        structure = parser.parse_string(sample_pdb_content)

        # Check if bonds were parsed
        assert hasattr(structure, "_bonds")
        bonds = structure._bonds

        if bonds is not None:
            assert len(bonds) > 0
            # Should have bonds from CONECT records
            # Bond between atoms 0 and 1 (N-CA)
            # Bond between atoms 1 and 2 (CA-C)
            # etc.

    def test_optional_annotations(self, sample_pdb_content):
        """Test that optional annotations are properly initialized."""
        parser = PDBParser()
        structure = parser.parse_string(sample_pdb_content)

        # Should have serial numbers (needed for CONECT)
        assert structure.serial is not None
        assert len(structure.serial) == 6
        assert structure.serial[0] == 1

        # Should have occupancy and b_factor data
        assert structure.occupancy is not None
        assert structure.b_factor is not None

        # Check occupancy values
        assert structure.occupancy[0] == 1.0
        assert structure.b_factor[0] == 20.0

    def test_structure_from_pdb_classmethod(self, sample_pdb_content):
        """Test Structure.from_pdb_string class method."""
        structure = Structure.from_pdb_string(sample_pdb_content)

        assert structure.n_atoms == 6
        assert structure.atom_name[0] == "N"
        assert structure.chain_id[0] == "A"
        assert structure.res_id[0] == 1

    def test_empty_pdb_raises_error(self):
        """Test that empty PDB content raises appropriate error."""
        parser = PDBParser()

        with pytest.raises(ValueError, match="No atom records found"):
            parser.parse_string("END")

    def test_malformed_coordinates_handled(self):
        """Test handling of malformed coordinate data."""
        malformed_pdb = """ATOM      1  N   ALA A   1      XXX   16.967  22.478  1.00 20.00           N  
END"""
        parser = PDBParser()
        structure = parser.parse_string(malformed_pdb)

        # Should have default coordinate value
        assert structure.coord[0][0] == 0.0  # Default for malformed X coordinate
        assert structure.coord[0][1] == 16.967  # Valid Y coordinate

    def test_missing_element_inferred(self):
        """Test that missing element symbols are inferred from atom names."""
        pdb_no_element = """ATOM      1  CA  ALA A   1      20.154  16.967  22.478  1.00 20.00              
END"""
        parser = PDBParser()
        structure = parser.parse_string(pdb_no_element)

        # Element should be inferred as C from CA atom name
        assert structure.element[0] == "C"

    def test_multi_model_returns_ensemble(self, sample_pdb_with_models):
        """Test that multi-model PDB returns StructureEnsemble."""
        parser = PDBParser()
        result = parser.parse_string(sample_pdb_with_models)

        # Should return StructureEnsemble
        assert isinstance(result, StructureEnsemble)
        assert result.n_frames == 2
        assert result.n_atoms == 2

        # Check that frames have different coordinates
        frame0 = result[0]
        frame1 = result[1]

        assert frame0.n_atoms == 2
        assert frame1.n_atoms == 2

        # Coordinates should be different between frames
        assert not np.array_equal(frame0.coord, frame1.coord)

        # But atom names should be the same
        assert list(frame0.atom_name) == list(frame1.atom_name)

    def test_structure_from_pdb_rejects_multimodel(self, sample_pdb_with_models):
        """Test that Structure.from_pdb_string rejects multi-model PDB."""
        with pytest.raises(ValueError, match="multiple models"):
            Structure.from_pdb_string(sample_pdb_with_models)

    def test_ensemble_from_pdb_rejects_single_model(self, sample_pdb_content):
        """Test that StructureEnsemble.from_pdb_string rejects single-model PDB."""
        with pytest.raises(ValueError, match="single model"):
            StructureEnsemble.from_pdb_string(sample_pdb_content)
