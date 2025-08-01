"""
Selection engine for evaluating selection expressions.

This module provides the main interface for atom selection, integrating
the expression system and parser to provide a simple API for structure
selection operations.
"""

from typing import Dict, List, Optional, Union

import numpy as np

from ..core.structure import Structure
from .expressions import SelectionExpression
from .parser import SelectionParser


class SelectionEngine:
    """
    Engine for evaluating atom selections on structures.

    Provides caching and optimization for repeated selections.
    """

    def __init__(self, cache_size: int = 100):
        """
        Initialize selection engine.

        Args:
            cache_size: Maximum number of cached selections
        """
        self.parser = SelectionParser()
        self.cache: Dict[str, SelectionExpression] = {}
        self.cache_size = cache_size

    def select(
        self, structure: Structure, selection: Union[str, SelectionExpression]
    ) -> np.ndarray:
        """
        Select atoms from a structure.

        Args:
            structure: The structure to select from
            selection: Selection string or expression

        Returns:
            Boolean array indicating selected atoms

        Raises:
            ParseException: If selection string is invalid
        """
        if isinstance(selection, str):
            # Parse string to expression
            expr = self._get_or_parse(selection)
        else:
            expr = selection

        # Evaluate expression
        return expr.evaluate(structure)

    def select_atoms(
        self, structure: Structure, selection: Union[str, SelectionExpression]
    ) -> Structure:
        """
        Return a new Structure containing only selected atoms.

        Args:
            structure: The structure to select from
            selection: Selection string or expression

        Returns:
            New Structure with selected atoms
        """
        mask = self.select(structure, selection)
        return structure[mask]

    def count(
        self, structure: Structure, selection: Union[str, SelectionExpression]
    ) -> int:
        """
        Count atoms matching selection.

        Args:
            structure: The structure to select from
            selection: Selection string or expression

        Returns:
            Number of selected atoms
        """
        mask = self.select(structure, selection)
        return int(np.sum(mask))

    def get_indices(
        self, structure: Structure, selection: Union[str, SelectionExpression]
    ) -> np.ndarray:
        """
        Get indices of atoms matching selection.

        Args:
            structure: The structure to select from
            selection: Selection string or expression

        Returns:
            Array of atom indices
        """
        mask = self.select(structure, selection)
        return np.where(mask)[0]  # type: ignore[no-any-return]

    def _get_or_parse(self, selection_string: str) -> SelectionExpression:
        """
        Get expression from cache or parse it.

        Args:
            selection_string: Selection string to parse

        Returns:
            Parsed SelectionExpression
        """
        if selection_string in self.cache:
            return self.cache[selection_string]

        # Parse expression
        expr = self.parser.parse(selection_string)

        # Add to cache (with size limit)
        if len(self.cache) >= self.cache_size:
            # Remove oldest entry (simple FIFO)
            oldest = next(iter(self.cache))
            del self.cache[oldest]

        self.cache[selection_string] = expr
        return expr

    def clear_cache(self) -> None:
        """Clear the selection cache."""
        self.cache.clear()


# Convenience functions


def select(
    structure: Structure, selection: Union[str, SelectionExpression]
) -> np.ndarray:
    """
    Select atoms from a structure.

    Args:
        structure: The structure to select from
        selection: Selection string or expression

    Returns:
        Boolean array indicating selected atoms
    """
    engine = SelectionEngine()
    return engine.select(structure, selection)


def select_atoms(
    structure: Structure, selection: Union[str, SelectionExpression]
) -> Structure:
    """
    Return a new Structure containing only selected atoms.

    Args:
        structure: The structure to select from
        selection: Selection string or expression

    Returns:
        New Structure with selected atoms
    """
    engine = SelectionEngine()
    return engine.select_atoms(structure, selection)
