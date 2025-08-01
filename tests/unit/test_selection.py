"""
Unit tests for the atom selection system.

Tests cover functionality based on acceptance criteria from requirements.md:
- AC-046 through AC-100 (Selection language basic functionality)
- AC-101 through AC-150 (Selection parsing)
"""

import numpy as np
import pytest
from pyparsing import ParseException

from molr.core.structure import Structure
from molr.selection import (  # Expression classes
    AllExpression,
    AndExpression,
    AromaticExpression,
    AtomNameExpression,
    BackboneExpression,
    ByResidueExpression,
    ChainExpression,
    DNAExpression,
    ElementExpression,
    IndexExpression,
    LigandExpression,
    NoneExpression,
    NotExpression,
    NucleicExpression,
    OrExpression,
    ProteinExpression,
    ResidueIdExpression,
    ResidueNameExpression,
    RNAExpression,
    SelectionEngine,
    SelectionParser,
    SidechainExpression,
    WaterExpression,
    select,
    select_atoms,
)


@pytest.fixture
def sample_structure():
    """Create a sample structure for testing."""
    structure = Structure(n_atoms=10)

    # Set up a mixed system: protein, DNA, water, ligand
    structure.atom_name[:] = ["N", "CA", "C", "O", "CB", "P", "C1'", "OH2", "C1", "C2"]
    structure.element[:] = ["N", "C", "C", "O", "C", "P", "C", "O", "C", "C"]
    structure.res_name[:] = [
        "ALA",
        "ALA",
        "ALA",
        "ALA",
        "ALA",
        "DA",
        "DA",
        "HOH",
        "LIG",
        "LIG",
    ]
    structure.res_id[:] = [1, 1, 1, 1, 1, 2, 2, 3, 4, 4]
    structure.chain_id[:] = ["A", "A", "A", "A", "A", "B", "B", "W", "L", "L"]

    return structure


class TestBasicSelectionExpressions:
    """Test basic selection expressions (AC-046 to AC-055)."""

    def test_all_expression(self, sample_structure):
        """Test selecting all atoms."""
        expr = AllExpression()
        mask = expr.evaluate(sample_structure)
        assert np.all(mask)
        assert len(mask) == sample_structure.n_atoms

    def test_none_expression(self, sample_structure):
        """Test selecting no atoms."""
        expr = NoneExpression()
        mask = expr.evaluate(sample_structure)
        assert not np.any(mask)
        assert len(mask) == sample_structure.n_atoms

    def test_element_expression_single(self, sample_structure):
        """Test selecting by single element."""
        expr = ElementExpression("C")
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, True, True, False, True, False, True, False, True, True]
        )
        assert np.array_equal(mask, expected)

    def test_element_expression_multiple(self, sample_structure):
        """Test selecting by multiple elements."""
        expr = ElementExpression(["C", "N"])
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, True, True, False, True, False, True, False, True, True]
        )
        assert np.array_equal(mask, expected)

    def test_atom_name_expression(self, sample_structure):
        """Test selecting by atom name."""
        expr = AtomNameExpression("CA")
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, True, False, False, False, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_residue_name_expression(self, sample_structure):
        """Test selecting by residue name."""
        expr = ResidueNameExpression("ALA")
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_residue_id_expression_single(self, sample_structure):
        """Test selecting by single residue ID."""
        expr = ResidueIdExpression(1)
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_residue_id_expression_range(self, sample_structure):
        """Test selecting by residue ID range."""
        expr = ResidueIdExpression(range(2, 4))
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, True, True, True, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_chain_expression(self, sample_structure):
        """Test selecting by chain."""
        expr = ChainExpression("A")
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_index_expression(self, sample_structure):
        """Test selecting by index."""
        expr = IndexExpression([0, 2, 4])
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, False, True, False, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)


class TestPropertySelectionExpressions:
    """Test property-based selection expressions (AC-056 to AC-065)."""

    def test_backbone_expression(self, sample_structure):
        """Test selecting backbone atoms."""
        expr = BackboneExpression()
        mask = expr.evaluate(sample_structure)
        # N, CA, C, O are backbone atoms in protein
        # P, C1' are backbone atoms in DNA
        expected_indices = [0, 1, 2, 3, 5, 6]
        expected = np.zeros(sample_structure.n_atoms, dtype=bool)
        expected[expected_indices] = True
        assert np.array_equal(mask, expected)

    def test_sidechain_expression(self, sample_structure):
        """Test selecting sidechain atoms."""
        expr = SidechainExpression()
        mask = expr.evaluate(sample_structure)
        # CB is sidechain in protein
        expected_indices = [4]
        expected = np.zeros(sample_structure.n_atoms, dtype=bool)
        expected[expected_indices] = True
        assert np.array_equal(mask, expected)

    def test_protein_expression(self, sample_structure):
        """Test selecting protein atoms."""
        expr = ProteinExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_nucleic_expression(self, sample_structure):
        """Test selecting nucleic acid atoms."""
        expr = NucleicExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_dna_expression(self, sample_structure):
        """Test selecting DNA atoms."""
        expr = DNAExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_ligand_expression(self, sample_structure):
        """Test selecting ligand atoms."""
        expr = LigandExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, False, False, False, True, True]
        )
        assert np.array_equal(mask, expected)

    def test_water_expression(self, sample_structure):
        """Test selecting water atoms."""
        expr = WaterExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, False, False, True, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_aromatic_expression(self, sample_structure):
        """Test selecting aromatic atoms."""
        # Modify structure to have aromatic residue
        sample_structure.res_name[5:7] = "PHE"
        expr = AromaticExpression()
        mask = expr.evaluate(sample_structure)
        expected = np.array(
            [False, False, False, False, False, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)


class TestBooleanOperations:
    """Test boolean operations on selections (AC-066 to AC-075)."""

    def test_and_expression(self, sample_structure):
        """Test AND operation."""
        expr1 = ElementExpression("C")
        expr2 = ResidueNameExpression("ALA")
        and_expr = AndExpression(expr1, expr2)
        mask = and_expr.evaluate(sample_structure)
        # C atoms in ALA residue
        expected = np.array(
            [False, True, True, False, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_or_expression(self, sample_structure):
        """Test OR operation."""
        expr1 = ElementExpression("N")
        expr2 = ElementExpression("O")
        or_expr = OrExpression(expr1, expr2)
        mask = or_expr.evaluate(sample_structure)
        # N or O atoms
        expected = np.array(
            [True, False, False, True, False, False, False, True, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_not_expression(self, sample_structure):
        """Test NOT operation."""
        expr = ElementExpression("C")
        not_expr = NotExpression(expr)
        mask = not_expr.evaluate(sample_structure)
        # Not C atoms
        expected = np.array(
            [True, False, False, True, False, True, False, True, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_operator_overloads(self, sample_structure):
        """Test operator overloads for boolean operations."""
        expr1 = ElementExpression("C")
        expr2 = ResidueNameExpression("ALA")
        expr3 = ElementExpression("N")

        # Test & operator
        and_expr = expr1 & expr2
        assert isinstance(and_expr, AndExpression)

        # Test | operator
        or_expr = expr1 | expr3
        assert isinstance(or_expr, OrExpression)

        # Test ~ operator
        not_expr = ~expr1
        assert isinstance(not_expr, NotExpression)

    def test_complex_boolean_expression(self, sample_structure):
        """Test complex boolean combinations."""
        # (protein and backbone) or water
        protein = ProteinExpression()
        backbone = BackboneExpression()
        water = WaterExpression()

        expr = (protein & backbone) | water
        mask = expr.evaluate(sample_structure)

        # Protein backbone atoms (indices 0-3) or water (index 7)
        expected = np.array(
            [True, True, True, True, False, False, False, True, False, False]
        )
        assert np.array_equal(mask, expected)


class TestByResidueSelection:
    """Test by-residue selection (AC-076 to AC-080)."""

    def test_byres_basic(self, sample_structure):
        """Test basic by-residue selection."""
        # Select residues containing CA atoms
        atom_expr = AtomNameExpression("CA")
        byres_expr = ByResidueExpression(atom_expr)
        mask = byres_expr.evaluate(sample_structure)

        # All atoms in residue 1 (which contains CA)
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_byres_multiple_residues(self, sample_structure):
        """Test by-residue selection across multiple residues."""
        # Select residues containing P atoms
        atom_expr = ElementExpression("P")
        byres_expr = ByResidueExpression(atom_expr)
        mask = byres_expr.evaluate(sample_structure)

        # All atoms in residue 2 (which contains P)
        expected = np.array(
            [False, False, False, False, False, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)


class TestSelectionParser:
    """Test selection language parser (AC-101 to AC-120)."""

    def test_parse_simple_selections(self, sample_structure):
        """Test parsing simple selection strings."""
        parser = SelectionParser()

        # Test various simple selections
        expr = parser.parse("all")
        assert isinstance(expr, AllExpression)

        expr = parser.parse("protein")
        assert isinstance(expr, ProteinExpression)

        expr = parser.parse("element C")
        assert isinstance(expr, ElementExpression)
        assert expr.elements == ["C"]

        expr = parser.parse("resname ALA")
        assert isinstance(expr, ResidueNameExpression)
        assert expr.resnames == ["ALA"]

    def test_parse_multiple_values(self, sample_structure):
        """Test parsing selections with multiple values."""
        parser = SelectionParser()

        expr = parser.parse("element C N O")
        assert isinstance(expr, ElementExpression)
        assert sorted(expr.elements) == ["C", "N", "O"]

        expr = parser.parse("resname ALA GLY VAL")
        assert isinstance(expr, ResidueNameExpression)
        assert expr.resnames == ["ALA", "GLY", "VAL"]

    def test_parse_residue_ranges(self, sample_structure):
        """Test parsing residue ID ranges."""
        parser = SelectionParser()

        expr = parser.parse("resid 1:10")
        assert isinstance(expr, ResidueIdExpression)
        assert expr.resids == list(range(1, 11))

        expr = parser.parse("resid 1 5 10")
        assert isinstance(expr, ResidueIdExpression)
        assert expr.resids == [1, 5, 10]

    def test_parse_boolean_operations(self, sample_structure):
        """Test parsing boolean operations."""
        parser = SelectionParser()

        # AND operation
        expr = parser.parse("protein and backbone")
        assert isinstance(expr, AndExpression)
        assert isinstance(expr.left, ProteinExpression)
        assert isinstance(expr.right, BackboneExpression)

        # OR operation
        expr = parser.parse("protein or ligand")
        assert isinstance(expr, OrExpression)

        # NOT operation
        expr = parser.parse("not water")
        assert isinstance(expr, NotExpression)
        assert isinstance(expr.operand, WaterExpression)

    def test_parse_implicit_and(self, sample_structure):
        """Test parsing implicit AND operations."""
        parser = SelectionParser()

        # Space implies AND
        expr = parser.parse("protein backbone")
        assert isinstance(expr, AndExpression)

    def test_parse_parentheses(self, sample_structure):
        """Test parsing with parentheses."""
        parser = SelectionParser()

        expr = parser.parse("(protein and chain A) or ligand")
        assert isinstance(expr, OrExpression)
        assert isinstance(expr.left, AndExpression)
        assert isinstance(expr.right, LigandExpression)

    def test_parse_byres(self, sample_structure):
        """Test parsing by-residue selections."""
        parser = SelectionParser()

        expr = parser.parse("byres (ligand)")
        assert isinstance(expr, ByResidueExpression)
        assert isinstance(expr.atom_selection, LigandExpression)

    def test_parse_aliases(self, sample_structure):
        """Test parsing with keyword aliases."""
        parser = SelectionParser()

        # Test element aliases
        expr1 = parser.parse("element C")
        expr2 = parser.parse("elem C")
        assert type(expr1) == type(expr2)

        # Test name aliases
        expr1 = parser.parse("name CA")
        expr2 = parser.parse("atomname CA")
        assert type(expr1) == type(expr2)

    def test_parse_errors(self):
        """Test parser error handling."""
        parser = SelectionParser()

        with pytest.raises(ParseException):
            parser.parse("invalid_selection_keyword")

        with pytest.raises(ParseException):
            parser.parse("element")  # Missing value

        with pytest.raises(ParseException):
            parser.parse("and protein")  # Missing left operand


class TestSelectionEngine:
    """Test selection engine functionality (AC-121 to AC-130)."""

    def test_engine_select(self, sample_structure):
        """Test basic selection through engine."""
        engine = SelectionEngine()

        mask = engine.select(sample_structure, "protein")
        expected = np.array(
            [True, True, True, True, True, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_engine_select_atoms(self, sample_structure):
        """Test selecting atoms to create new structure."""
        engine = SelectionEngine()

        subset = engine.select_atoms(sample_structure, "protein")
        assert subset.n_atoms == 5
        assert np.all(subset.res_name == "ALA")

    def test_engine_count(self, sample_structure):
        """Test counting selected atoms."""
        engine = SelectionEngine()

        count = engine.count(sample_structure, "element C")
        assert count == 6

    def test_engine_get_indices(self, sample_structure):
        """Test getting indices of selected atoms."""
        engine = SelectionEngine()

        indices = engine.get_indices(sample_structure, "backbone")
        expected = np.array([0, 1, 2, 3, 5, 6])
        assert np.array_equal(indices, expected)

    def test_engine_caching(self, sample_structure):
        """Test that engine caches parsed selections."""
        engine = SelectionEngine(cache_size=10)

        # First call should parse
        mask1 = engine.select(sample_structure, "protein and backbone")
        assert len(engine.cache) == 1

        # Second call should use cache
        mask2 = engine.select(sample_structure, "protein and backbone")
        assert len(engine.cache) == 1
        assert np.array_equal(mask1, mask2)

    def test_engine_cache_limit(self, sample_structure):
        """Test cache size limiting."""
        engine = SelectionEngine(cache_size=2)

        engine.select(sample_structure, "protein")
        engine.select(sample_structure, "ligand")
        assert len(engine.cache) == 2

        # This should evict the oldest entry
        engine.select(sample_structure, "water")
        assert len(engine.cache) == 2
        assert "protein" not in engine.cache

    def test_convenience_functions(self, sample_structure):
        """Test module-level convenience functions."""
        # Test select function
        mask = select(sample_structure, "element C")
        assert isinstance(mask, np.ndarray)
        assert mask.dtype == bool

        # Test select_atoms function
        subset = select_atoms(sample_structure, "chain A")
        assert isinstance(subset, Structure)
        assert subset.n_atoms == 5


class TestStructureIntegration:
    """Test selection integration with Structure class (AC-131 to AC-140)."""

    def test_structure_select_method(self, sample_structure):
        """Test Structure.select() method."""
        mask = sample_structure.select("protein and backbone")
        expected = np.array(
            [True, True, True, True, False, False, False, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_structure_select_complex(self, sample_structure):
        """Test complex selection through Structure."""
        mask = sample_structure.select("(element C or element N) and not ligand")
        # C or N atoms that are not in ligand
        expected = np.array(
            [True, True, True, False, True, False, True, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_structure_indexing_with_selection(self, sample_structure):
        """Test using selection result for structure indexing."""
        protein_structure = sample_structure[sample_structure.select("protein")]
        assert protein_structure.n_atoms == 5
        assert np.all(protein_structure.res_name == "ALA")


class TestSelectionExamples:
    """Test real-world selection examples (AC-141 to AC-150)."""

    def test_active_site_selection(self, sample_structure):
        """Test selecting active site residues."""
        # Example: protein residues 1-5
        mask = sample_structure.select("protein and resid 1:5")
        assert np.sum(mask) == 5  # All protein atoms are in residue 1

    def test_interface_selection(self, sample_structure):
        """Test selecting interface atoms."""
        # Example: protein chain A or DNA chain B
        mask = sample_structure.select("(protein and chain A) or (dna and chain B)")
        expected = np.array(
            [True, True, True, True, True, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)

    def test_hydrophobic_selection(self, sample_structure):
        """Test selecting hydrophobic atoms."""
        # Add hydrophobic residue for testing
        sample_structure.res_name[0:5] = "VAL"
        mask = sample_structure.select("resname VAL ALA LEU ILE")
        assert np.sum(mask) == 5  # All VAL atoms

    def test_metal_coordination_selection(self, sample_structure):
        """Test selecting metal coordination sphere."""
        # This would require distance-based selection (future feature)
        # For now, test selecting by residue containing specific atoms
        mask = sample_structure.select("byres (element P)")
        expected = np.array(
            [False, False, False, False, False, True, True, False, False, False]
        )
        assert np.array_equal(mask, expected)
