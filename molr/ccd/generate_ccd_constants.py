#!/usr/bin/env python3
"""
Main script for generating CCD bond constants.

This script is designed to be used by the Makefile and provides a simple
interface for generating residue bond constants from CCD BinaryCIF files.

Can be run as:
1. python -m molr.ccd.generate_ccd_constants (recommended)
2. python molr/ccd/generate_ccd_constants.py (with path handling)
"""

import os
import sys

# Handle both direct execution and module execution
try:
    # Try relative imports first (when run as module)
    from ..constants.pdb_constants import RESIDUES
    from .ccd_analyzer import CCDDataManager
    from .constants_generator import CCDConstantsGenerator
except ImportError:
    # Fall back to absolute imports (when run as script)
    project_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    sys.path.insert(0, project_root)

    from molr.ccd.ccd_analyzer import CCDDataManager
    from molr.ccd.constants_generator import CCDConstantsGenerator
    from molr.constants.pdb_constants import RESIDUES


def main() -> int:
    """
    Main function for generating CCD bond constants.

    This function:
    1. Initializes the CCD data manager
    2. Downloads CCD files if needed
    3. Generates bond constants for all standard residues
    4. Provides a comprehensive analysis report
    """
    print("=== CCD Bond Constants Generator ===")
    print("Initializing CCD data manager...")

    # Initialize CCD data manager (uses ~/.molr/ccd-data by default)
    ccd_manager = CCDDataManager()

    # Ensure files are available (download if needed)
    if not ccd_manager.ensure_files_exist():
        print("Error: Failed to ensure CCD files are available")
        return 1

    # Initialize constants generator
    generator = CCDConstantsGenerator(ccd_manager)

    # Generate analysis report
    print("\nGenerating analysis report...")
    report = generator.generate_analysis_report(RESIDUES)
    generator.print_analysis_report(report)

    # Generate constants file
    print("\n=== Generating Bond Constants ===")
    success = generator.write_residue_bonds_constants(RESIDUES)

    if success:
        print("\n✅ Bond constants generated successfully!")
        print("📁 Output: molr/constants/residue_bonds.py")
        print(
            f"📊 Processed {report['found_residues']} residues with {report['total_bonds']:,} bonds"
        )
        return 0
    else:
        print("\n❌ Failed to generate bond constants")
        return 1


if __name__ == "__main__":
    sys.exit(main())
