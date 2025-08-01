"""
MolR - Molecular Realm for Spatial Indexed Structures.

A Python package for efficient molecular structure representation with built-in
spatial indexing, providing fast neighbor queries and geometric operations.

Key Features:
- NumPy-based Structure class with Structure of Arrays design
- Efficient spatial indexing with KDTree integration
- Comprehensive bond detection with multiple providers
- Selection language for complex atom queries
- Support for PDB and mmCIF file formats
- Memory-efficient trajectory handling
"""

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.0.0+unknown"
__author__ = "Abhishek Tiwari"

# Import bond detection
from .bond_detection import detect_bonds
from .core.bond_list import BondList, BondOrder
from .core.connectivity import BondDetectionMethod
from .core.realm import Realm

# Import core classes
from .core.structure import Structure
from .core.structure_ensemble import StructureEnsemble
from .io.mmcif import mmCIFParser

# Import I/O parsers
from .io.pdb import PDBParser

# Import selection functionality
from .selection import select

__all__ = [
    # Core classes
    "Structure",
    "StructureEnsemble",
    "BondList",
    "BondOrder",
    "BondDetectionMethod",
    "Realm",
    # Selection
    "select",
    # I/O
    "PDBParser",
    "mmCIFParser",
    # Bond detection
    "detect_bonds",
]
