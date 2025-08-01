"""
Unit tests for mmCIF parser.

Tests the mmCIF parser functionality including parsing files,
string content, multi-model support, and bond information extraction.
"""

import os
import tempfile
from pathlib import Path

import pytest

from molr.core.bond_list import BondList
from molr.core.structure import Structure
from molr.core.structure_ensemble import StructureEnsemble
from molr.io.mmcif import mmCIFParser


@pytest.fixture
def simple_mmcif_content():
    """Simple mmCIF content for testing."""
    return """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N  . ALA A 1 1   ? 20.154 16.967 22.478 1.00 8.65  ? 1   ALA A N  1
ATOM   2  C  CA . ALA A 1 1   ? 20.238 18.376 22.896 1.00 8.63  ? 1   ALA A CA 1
ATOM   3  C  C  . ALA A 1 1   ? 21.652 18.840 22.496 1.00 8.55  ? 1   ALA A C  1
ATOM   4  O  O  . ALA A 1 1   ? 22.041 18.819 21.327 1.00 8.64  ? 1   ALA A O  1
ATOM   5  C  CB . ALA A 1 1   ? 19.191 19.321 22.291 1.00 8.80  ? 1   ALA A CB 1
#
"""


@pytest.fixture
def multimodel_mmcif_content():
    """Multi-model mmCIF content for testing."""
    return """data_multimodel
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N  . ALA A 1 1   ? 20.154 16.967 22.478 1.00 8.65  ? 1   ALA A N  1
ATOM   2  C  CA . ALA A 1 1   ? 20.238 18.376 22.896 1.00 8.63  ? 1   ALA A CA 1
ATOM   3  N  N  . ALA A 1 1   ? 20.154 16.967 22.478 1.00 8.65  ? 1   ALA A N  2
ATOM   4  C  CA . ALA A 1 1   ? 20.238 18.376 22.896 1.00 8.63  ? 1   ALA A CA 2
#
"""


class TestmmCIFParser:
    """Test mmCIF parser functionality."""

    def test_parser_initialization(self):
        """Test parser initialization."""
        parser = mmCIFParser()
        assert parser._current_structure is None
        assert parser._current_bonds is None

    def test_parse_string_simple(self, simple_mmcif_content):
        """Test parsing simple mmCIF content from string."""
        parser = mmCIFParser()

        result = parser.parse_string(simple_mmcif_content)

        assert isinstance(result, Structure)
        assert result.n_atoms == 5

        # Check atom data
        assert result.atom_name[0] == "N"
        assert result.atom_name[1] == "CA"
        assert result.element[0] == "N"
        assert result.element[1] == "C"
        assert result.res_name[0] == "ALA"
        assert result.chain_id[0] == "A"
        assert result.res_id[0] == 1

        # Check coordinates
        assert abs(result.coord[0][0] - 20.154) < 0.001
        assert abs(result.coord[0][1] - 16.967) < 0.001
        assert abs(result.coord[0][2] - 22.478) < 0.001

        # Check optional annotations
        assert result.occupancy is not None
        assert result.b_factor is not None
        assert abs(result.occupancy[0] - 1.00) < 0.001
        assert abs(result.b_factor[0] - 8.65) < 0.001

    def test_parse_string_multimodel(self, multimodel_mmcif_content):
        """Test parsing multi-model mmCIF content."""
        parser = mmCIFParser()

        result = parser.parse_string(multimodel_mmcif_content)

        assert isinstance(result, StructureEnsemble)
        assert result.n_frames == 2
        assert result.template.n_atoms == 2

        # Check template structure
        assert result.template.atom_name[0] == "N"
        assert result.template.atom_name[1] == "CA"
        assert result.template.res_name[0] == "ALA"

        # Check coordinates for both models
        # Note: coords array has capacity >= n_frames, so check only valid frames
        assert result.coords.shape[1:] == (2, 3)  # (n_atoms, 3)
        assert result.coords.shape[0] >= 2  # capacity >= n_frames
        assert abs(result.coords[0][0][0] - 20.154) < 0.001  # Model 1, atom 1, x
        assert abs(result.coords[1][0][0] - 20.154) < 0.001  # Model 2, atom 1, x

    def test_parse_file_functionality(self, simple_mmcif_content):
        """Test file parsing functionality."""
        parser = mmCIFParser()

        # Create temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
            f.write(simple_mmcif_content)
            temp_filename = f.name

        try:
            result = parser.parse_file(temp_filename)

            assert isinstance(result, Structure)
            assert result.n_atoms == 5
            assert result.atom_name[0] == "N"

        finally:
            os.unlink(temp_filename)

    def test_safe_conversion_functions(self):
        """Test safe conversion utility functions."""
        from molr.io.mmcif import (
            _safe_convert_float,
            _safe_convert_int,
            _safe_convert_str,
        )

        # Test integer conversion
        assert _safe_convert_int("123") == 123
        assert _safe_convert_int("?") == 0
        assert _safe_convert_int(".") == 0
        assert _safe_convert_int(None) == 0
        assert _safe_convert_int("invalid", 99) == 99

        # Test float conversion
        assert abs(_safe_convert_float("12.5") - 12.5) < 0.001
        assert _safe_convert_float("?") == 0.0
        assert _safe_convert_float(".") == 0.0
        assert _safe_convert_float(None) == 0.0
        assert abs(_safe_convert_float("invalid", 99.9) - 99.9) < 0.001

        # Test string conversion
        assert _safe_convert_str("test") == "test"
        assert _safe_convert_str("?") == ""
        assert _safe_convert_str(".") == ""
        assert _safe_convert_str(None) == ""
        assert _safe_convert_str("  test  ") == "test"  # Should strip

    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        parser = mmCIFParser()

        # Test invalid file
        with pytest.raises(IOError):
            parser.parse_file("nonexistent_file.cif")

        # Test empty content
        with pytest.raises(ValueError):
            parser.parse_string("")

        # Test invalid mmCIF format
        with pytest.raises(ValueError):
            parser.parse_string("invalid mmcif content")

    def test_hetatm_handling(self):
        """Test handling of HETATM records."""
        hetatm_content = """data_hetatm_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1  O  O  . HOH W 2 .   ? 25.000 25.000 25.000 1.00 10.00 ? .   HOH W O  1
HETATM 2  O  O1 . LIG L 3 1   ? 30.000 30.000 30.000 1.00 15.00 ? 1   LIG L O1 1
#
"""

        parser = mmCIFParser()
        result = parser.parse_string(hetatm_content)

        assert isinstance(result, Structure)
        assert result.n_atoms == 2
        assert result.atom_name[0] == "O"
        assert result.atom_name[1] == "O1"
        assert result.res_name[0] == "HOH"
        assert result.res_name[1] == "LIG"

    def test_alternative_locations(self):
        """Test handling of alternative locations."""
        alt_loc_content = """data_alt_loc_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N  . ALA A 1 1   ? 20.154 16.967 22.478 1.00 8.65  ? 1   ALA A N  1
ATOM   2  C  CA A ALA A 1 1   ? 20.238 18.376 22.896 1.00 8.63  ? 1   ALA A CA 1
ATOM   3  C  CA B ALA A 1 1   ? 20.300 18.400 22.900 1.00 8.63  ? 1   ALA A CA 1
#
"""

        parser = mmCIFParser()
        result = parser.parse_string(alt_loc_content)

        assert isinstance(result, Structure)
        assert result.n_atoms == 3
        assert result.alt_loc is not None
        assert result.alt_loc[0] == ""  # No alt loc
        assert result.alt_loc[1] == "A"
        assert result.alt_loc[2] == "B"

    def test_charge_handling(self):
        """Test handling of formal charges."""
        charged_content = """data_charge_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N  . LYS A 1 1   ? 20.154 16.967 22.478 1.00 8.65  +1 1   LYS A N  1
ATOM   2  O  O  . GLU A 1 2   ? 20.238 18.376 22.896 1.00 8.63  -1 2   GLU A O  1
#
"""

        parser = mmCIFParser()
        result = parser.parse_string(charged_content)

        assert isinstance(result, Structure)
        assert result.n_atoms == 2
        assert result.charge is not None
        assert abs(result.charge[0] - 1.0) < 0.001  # +1 charge
        assert abs(result.charge[1] - (-1.0)) < 0.001  # -1 charge

    def test_missing_data_handling(self):
        """Test handling of missing data fields."""
        minimal_content = """data_minimal_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM   1  C  CA ALA A 1   20.000 20.000 20.000 1
#
"""

        parser = mmCIFParser()
        result = parser.parse_string(minimal_content)

        assert isinstance(result, Structure)
        assert result.n_atoms == 1
        assert result.atom_name[0] == "CA"
        assert result.element[0] == "C"
        assert result.res_name[0] == "ALA"

        # Should have default values for missing fields
        assert result.occupancy[0] == 1.0  # Default occupancy
        assert result.b_factor[0] == 0.0  # Default B-factor


class TestBondParsing:
    """Test bond parsing from mmCIF data."""

    def test_bond_order_mapping(self):
        """Test bond order string mapping."""
        from molr.io.mmcif import mmCIFParser

        parser = mmCIFParser()

        assert parser._map_bond_order("SING") == 1.0
        assert parser._map_bond_order("DOUB") == 2.0
        assert parser._map_bond_order("TRIP") == 3.0
        assert parser._map_bond_order("AROM") == 1.5
        assert parser._map_bond_order("UNKNOWN") == 1.0  # Default

    def test_atom_mapping_creation(self, simple_mmcif_content):
        """Test atom mapping creation."""
        parser = mmCIFParser()
        structure = parser.parse_string(simple_mmcif_content)

        atom_mapping = parser._create_atom_mapping(structure)

        assert "N" in atom_mapping
        assert "CA" in atom_mapping
        assert "C" in atom_mapping
        assert "O" in atom_mapping
        assert "CB" in atom_mapping

        # Check that indices are correct
        assert atom_mapping["N"] == 0
        assert atom_mapping["CA"] == 1


if __name__ == "__main__":
    pytest.main([__file__])
