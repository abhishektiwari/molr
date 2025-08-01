"""
Bond detection method enumeration.

This module defines the bond detection methods used in MolR.
The complex hierarchical system has been replaced with a simpler approach.
"""

from enum import Enum


class BondDetectionMethod(Enum):
    """Bond detection method enumeration."""

    FILE = "file"  # From file-based bonds (PDB CONECT, mmCIF, etc.)
    TEMPLATE = "template"  # From residue_bonds.py templates
    DISTANCE = "distance"  # Distance-based detection
    MANUAL = "manual"  # Manually added bonds
