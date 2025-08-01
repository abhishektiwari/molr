"""
Configuration for MolR package.

This module provides configuration paths and settings for the MolR package,
including paths for CCD data storage.
"""

import os
from pathlib import Path


def get_molr_data_dir() -> Path:
    """
    Get the MolR data directory path.

    Returns:
        Path to MolR data directory (~/.molr)
    """
    home = Path.home()
    molr_dir = home / ".molr"
    molr_dir.mkdir(exist_ok=True)
    return molr_dir


def get_ccd_data_path() -> Path:
    """
    Get the path for CCD data storage.

    Returns:
        Path to CCD data directory (~/.molr/ccd-data)
    """
    molr_dir = get_molr_data_dir()
    ccd_dir = molr_dir / "ccd-data"
    ccd_dir.mkdir(exist_ok=True)
    return ccd_dir
