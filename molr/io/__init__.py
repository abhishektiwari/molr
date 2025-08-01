"""
I/O package for reading and writing molecular structure files.

This package provides file format support for the space package, including
PDB and mmCIF reading, CONECT record parsing, and connectivity interfaces.
"""

from .pdb import PDBParser
from .mmcif import mmCIFParser

__all__ = [
    "PDBParser",
    "mmCIFParser",
]
