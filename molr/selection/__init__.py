"""
Atom selection language implementation.

This module provides a powerful atom selection language for molecular structures,
inspired by MDAnalysis and VMD selection syntax.
"""

from .engine import SelectionEngine, select, select_atoms
from .expressions import (
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
    SelectionExpression,
    SidechainExpression,
    WaterExpression,
)
from .parser import SelectionParser

__all__ = [
    # Engine and convenience functions
    "SelectionEngine",
    "select",
    "select_atoms",
    # Parser
    "SelectionParser",
    # Base expression
    "SelectionExpression",
    # Specific expressions
    "AllExpression",
    "NoneExpression",
    "ElementExpression",
    "AtomNameExpression",
    "ResidueNameExpression",
    "ResidueIdExpression",
    "ChainExpression",
    "BackboneExpression",
    "SidechainExpression",
    "ProteinExpression",
    "NucleicExpression",
    "DNAExpression",
    "RNAExpression",
    "LigandExpression",
    "AromaticExpression",
    "WaterExpression",
    "AndExpression",
    "OrExpression",
    "NotExpression",
    "IndexExpression",
    "ByResidueExpression",
]
