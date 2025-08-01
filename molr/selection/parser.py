"""
Selection language parser using pyparsing.

This module provides a parser for atom selection strings, converting them
into SelectionExpression objects that can be evaluated against structures.
The syntax is inspired by MDAnalysis and VMD selection languages.
"""

from typing import Any, List, Optional, Union

import pyparsing as pp
from pyparsing import (
    CaselessKeyword,
    CaselessLiteral,
    Forward,
    Group,
    OneOrMore,
)
from pyparsing import Optional as PPOptional
from pyparsing import (
    ParseException,
    Regex,
    Suppress,
    Word,
    ZeroOrMore,
    alphanums,
    alphas,
    nums,
    oneOf,
    pyparsing_common,
)

from .expressions import (
    AllExpression,
    AndExpression,
    AromaticExpression,
    AroundExpression,
    AtomNameExpression,
    BackboneExpression,
    ByResidueExpression,
    CenterOfGeometryExpression,
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
    SelectionExpression,
    SidechainExpression,
    SphericalExpression,
    WaterExpression,
    WithinExpression,
)


class SelectionParser:
    """
    Parser for atom selection language.

    Supports syntax like:
        - "protein and backbone"
        - "resname ALA GLY"
        - "chain A and resid 1:100"
        - "element C N O"
        - "not water"
        - "(protein and chain A) or ligand"
        - "byres (ligand and within 5 of protein)"
    """

    def __init__(self) -> None:
        """Initialize the parser with grammar rules."""
        self._build_grammar()

    def _build_grammar(self) -> None:
        """Build the pyparsing grammar for selection language."""
        # Enable packrat parsing for better performance
        pp.ParserElement.enablePackrat()

        # Basic tokens
        integer = pyparsing_common.integer
        real = pyparsing_common.real
        identifier = Word(alphas, alphanums + "_")

        # Keywords (case-insensitive)
        AND = CaselessKeyword("and")
        OR = CaselessKeyword("or")
        NOT = CaselessKeyword("not")

        # Selection keywords
        ALL = CaselessKeyword("all")
        NONE = CaselessKeyword("none")
        PROTEIN = CaselessKeyword("protein")
        NUCLEIC = CaselessKeyword("nucleic")
        DNA = CaselessKeyword("dna")
        RNA = CaselessKeyword("rna")
        BACKBONE = CaselessKeyword("backbone")
        SIDECHAIN = CaselessKeyword("sidechain")
        WATER = CaselessKeyword("water")
        LIGAND = CaselessKeyword("ligand")
        AROMATIC = CaselessKeyword("aromatic")

        # Property keywords
        ELEMENT = CaselessKeyword("element") | CaselessKeyword("elem")
        NAME = CaselessKeyword("name") | CaselessKeyword("atomname")
        RESNAME = CaselessKeyword("resname") | CaselessKeyword("resn")
        RESID = CaselessKeyword("resid") | CaselessKeyword("resi")
        CHAIN = CaselessKeyword("chain") | CaselessKeyword("segid")
        INDEX = CaselessKeyword("index") | CaselessKeyword("idx")
        BYRES = CaselessKeyword("byres") | CaselessKeyword("byresidue")

        # Spatial keywords
        WITHIN = CaselessKeyword("within")
        AROUND = CaselessKeyword("around")
        OF = CaselessKeyword("of")
        COG = CaselessKeyword("cog") | CaselessKeyword("centerofgeometry")
        SPHERE = CaselessKeyword("sphere")
        CENTER = CaselessKeyword("center")
        RADIUS = CaselessKeyword("radius")

        # Operators
        LPAREN = Suppress("(")
        RPAREN = Suppress(")")
        COLON = Suppress(":")

        # Forward declaration for recursive grammar
        selection_expr = Forward()

        # Simple selections
        all_selection = ALL.setParseAction(lambda: AllExpression())
        none_selection = NONE.setParseAction(lambda: NoneExpression())
        protein_selection = PROTEIN.setParseAction(lambda: ProteinExpression())
        nucleic_selection = NUCLEIC.setParseAction(lambda: NucleicExpression())
        dna_selection = DNA.setParseAction(lambda: DNAExpression())
        rna_selection = RNA.setParseAction(lambda: RNAExpression())
        backbone_selection = BACKBONE.setParseAction(lambda: BackboneExpression())
        sidechain_selection = SIDECHAIN.setParseAction(lambda: SidechainExpression())
        water_selection = WATER.setParseAction(lambda: WaterExpression())
        ligand_selection = LIGAND.setParseAction(lambda: LigandExpression())
        aromatic_selection = AROMATIC.setParseAction(lambda: AromaticExpression())

        # Element selection: "element C N O" or "elem C"
        element_list = Group(ELEMENT + OneOrMore(Word(alphas, max=2)))
        element_selection = element_list.setParseAction(
            lambda t: ElementExpression(list(t[0][1:]))
        )

        # Atom name selection: "name CA CB" or "atomname CA"
        name_list = Group(NAME + OneOrMore(Word(alphanums + "_*")))
        name_selection = name_list.setParseAction(
            lambda t: AtomNameExpression(list(t[0][1:]))
        )

        # Residue name selection: "resname ALA GLY" or "resn ALA"
        resname_list = Group(RESNAME + OneOrMore(Word(alphanums)))
        resname_selection = resname_list.setParseAction(
            lambda t: ResidueNameExpression(list(t[0][1:]))
        )

        # Residue ID selection: "resid 1 2 3" or "resid 1:100" or "resi 50"
        resid_range = Group(integer + COLON + integer)
        resid_value = resid_range | integer
        resid_list = Group(RESID + ZeroOrMore(resid_value))

        def parse_resid(tokens: Any) -> Any:
            """Parse residue ID selection."""
            values: List[int] = []
            i = 1
            while i < len(tokens[0]):
                current = tokens[0][i]
                if hasattr(current, "asList") and len(current.asList()) == 2:
                    # Range notation (from resid_range)
                    start, end = current.asList()
                    values.extend(range(start, end + 1))
                elif isinstance(current, list) and len(current) == 2:
                    # Range notation (backup)
                    start, end = current
                    values.extend(range(start, end + 1))
                else:
                    # Single value
                    values.append(current)
                i += 1
            return ResidueIdExpression(values)

        resid_selection = resid_list.setParseAction(parse_resid)

        # Chain selection: "chain A B" or "chain AB"
        chain_chars = Word(alphanums)
        chain_list = Group(CHAIN + ZeroOrMore(chain_chars))

        def parse_chain(tokens: Any) -> Any:
            """Parse chain selection."""
            chains = []
            for item in tokens[0][1:]:
                if len(item) == 1:
                    chains.append(item)
                else:
                    # Multiple chains in one string
                    chains.extend(list(item))
            return ChainExpression(chains)

        chain_selection = chain_list.setParseAction(parse_chain)

        # Index selection: "index 0 1 2" or "index 0:100"
        index_range = integer + COLON + integer + PPOptional(COLON + integer)
        index_value = index_range | integer
        index_list = Group(INDEX + ZeroOrMore(index_value))

        def parse_index(tokens: Any) -> Any:
            """Parse index selection."""
            values: List[Union[int, slice]] = []
            i = 1
            while i < len(tokens[0]):
                if i + 2 < len(tokens[0]) and isinstance(tokens[0][i + 1], str):
                    # Range notation
                    start = tokens[0][i]
                    end = tokens[0][i + 2]
                    if i + 4 < len(tokens[0]) and isinstance(tokens[0][i + 3], str):
                        # Step notation
                        step = tokens[0][i + 4]
                        values.append(slice(start, end, step))
                        i += 5
                    else:
                        values.extend(list(range(start, end)))
                        i += 3
                else:
                    # Single value
                    values.append(tokens[0][i])
                    i += 1

            if len(values) == 1 and isinstance(values[0], slice):
                return IndexExpression(values[0])
            # Convert to list of ints, expanding slices
            int_values = []
            for v in values:
                if isinstance(v, slice):
                    # Expand slice to list of ints
                    start = v.start if v.start is not None else 0
                    stop = v.stop if v.stop is not None else 0
                    step = v.step if v.step is not None else 1
                    int_values.extend(list(range(start, stop, step)))
                else:
                    int_values.append(v)
            return IndexExpression(int_values)

        index_selection = index_list.setParseAction(parse_index)

        # Spatial selections
        # "within 5.0 of protein"
        within_selection = (
            WITHIN + real + OF + LPAREN + selection_expr + RPAREN
        ).setParseAction(lambda t: WithinExpression(t[1], t[4]))

        # "around protein 5.0"
        around_selection = (
            AROUND + LPAREN + selection_expr + RPAREN + real
        ).setParseAction(lambda t: AroundExpression(t[2], t[4]))

        # "cog protein 8.0" (atoms within distance of center of geometry)
        cog_selection = (COG + LPAREN + selection_expr + RPAREN + real).setParseAction(
            lambda t: CenterOfGeometryExpression(t[2], t[4])
        )

        # By residue selection: "byres (protein and chain A)"
        byres_selection = (BYRES + LPAREN + selection_expr + RPAREN).setParseAction(
            lambda t: ByResidueExpression(t[1])
        )

        # Atomic selection (all non-boolean selections)
        atomic_selection = (
            all_selection
            | none_selection
            | protein_selection
            | nucleic_selection
            | dna_selection
            | rna_selection
            | backbone_selection
            | sidechain_selection
            | water_selection
            | ligand_selection
            | aromatic_selection
            | element_selection
            | name_selection
            | resname_selection
            | resid_selection
            | chain_selection
            | index_selection
            | within_selection
            | around_selection
            | cog_selection
            | byres_selection
            | (LPAREN + selection_expr + RPAREN)
        )

        # NOT expression
        not_expr = (NOT + atomic_selection).setParseAction(
            lambda t: NotExpression(t[1])
        )

        # AND expression (implicit when no operator)
        and_expr = (atomic_selection | not_expr) + ZeroOrMore(
            PPOptional(AND) + (atomic_selection | not_expr)
        )

        def parse_and(tokens: Any) -> Any:
            """Parse AND expressions."""
            if len(tokens) == 1:
                return tokens[0]
            result = tokens[0]
            for i in range(1, len(tokens)):
                if not isinstance(tokens[i], SelectionExpression):
                    continue
                result = AndExpression(result, tokens[i])
            return result

        and_expr.setParseAction(parse_and)

        # OR expression
        or_expr = and_expr + ZeroOrMore(OR + and_expr)

        def parse_or(tokens: Any) -> Any:
            """Parse OR expressions."""
            if len(tokens) == 1:
                return tokens[0]
            result = tokens[0]
            for i in range(2, len(tokens), 2):
                result = OrExpression(result, tokens[i])
            return result

        or_expr.setParseAction(parse_or)

        # Complete expression
        selection_expr <<= or_expr

        # Set the complete grammar
        self.grammar = selection_expr + pp.StringEnd()

    def parse(self, selection_string: str) -> SelectionExpression:
        """
        Parse a selection string into a SelectionExpression.

        Args:
            selection_string: The selection string to parse

        Returns:
            SelectionExpression object

        Raises:
            ParseException: If the string cannot be parsed
        """
        try:
            result = self.grammar.parseString(selection_string, parseAll=True)
            return result[0]  # type: ignore[no-any-return]
        except ParseException as e:
            # Enhance error message
            col = e.column
            line = e.line
            error_msg = f"Invalid selection syntax at position {col}: {e.msg}"
            if line:
                marker = " " * (col - 1) + "^"
                error_msg = f"{error_msg}\n{line}\n{marker}"
            raise ParseException(error_msg)

    @classmethod
    def parse_selection(cls, selection_string: str) -> SelectionExpression:
        """
        Convenience class method to parse a selection string.

        Args:
            selection_string: The selection string to parse

        Returns:
            SelectionExpression object
        """
        parser = cls()
        return parser.parse(selection_string)
