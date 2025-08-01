"""
Bond detection package for molecular structures.

This package provides a simplified bond detection approach:
- Default detection using residue templates and distance criteria
- File-based bond information from PDB CONECT and mmCIF
"""

from .default_detector import DefaultBondDetector, detect_bonds

__all__ = [
    "DefaultBondDetector",
    "detect_bonds",
]
