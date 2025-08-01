"""
End-to-end integration tests for bond detection using real PDB structures.

This module tests the complete bond detection pipeline using actual PDB files
to validate that FILE, TEMPLATE, and DISTANCE-based bond detection methods
work correctly on real molecular structures. Tests are data-driven with
expected value ranges for validation.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pytest

from molr.bond_detection.default_detector import DefaultBondDetector
from molr.core.bond_list import BondList
from molr.core.structure import Structure
from molr.io.mmcif import mmCIFParser
from molr.io.pdb import PDBParser

# Get path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "example_pdb_files"


# Test data with expected ranges - using dictionaries for pytest parametrization
TEST_DATA = [
    {
        "filename": "6rsa.pdb",
        "file_bonds": (65, 75),  # Around 70 expected
        "template_bonds": (1700, 1800),  # Around 1752 expected
        "distance_bonds": (200, 280),  # Around 237 expected
        "total_bonds": (2000, 2100),  # Around 2059 expected
        "hydrogen_bonds_min": 650,  # At least 650 H bonds
        "hydrogen_fraction_min": 0.3,  # At least 30% H bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Ribonuclease A with hydrogen atoms and CONECT records",
        "has_conect_records": True,
        "has_hydrogen_atoms": True,
    },
    {
        "filename": "2izf.pdb",
        "file_bonds": (80, 90),  # Around 84 expected
        "template_bonds": (2400, 2550),  # Around 2494 expected
        "distance_bonds": (100, 150),  # Around 122 expected
        "total_bonds": (2650, 2750),  # Around 2700 expected
        "hydrogen_bonds_min": 850,  # At least 850 H bonds
        "hydrogen_fraction_min": 0.3,  # At least 30% H bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Streptavidin-biotin complex with hydrogen atoms and CONECT records",
        "has_conect_records": True,
        "has_hydrogen_atoms": True,
    },
    {
        "filename": "1crn.pdb",
        # Small structure for quick tests - crambin (no hydrogens)
        "file_bonds": (0, 10),  # May have a few CONECT records
        "template_bonds": (280, 320),  # Around 289 expected
        "distance_bonds": (0, 50),  # Not many distance bonds expected
        "total_bonds": (280, 350),  # Around 295 total
        "hydrogen_bonds_min": None,  # No hydrogen bonds expected
        "hydrogen_fraction_min": None,  # No hydrogen bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Small crambin structure for quick testing (no hydrogens)",
        "has_conect_records": True,  # Actually has some CONECT records
        "has_hydrogen_atoms": False,  # No hydrogen atoms in this file
    },
    # CIF file versions - have more connectivity info from chemical bond tables
    {
        "filename": "6rsa.cif",
        "file_bonds": (320, 340),  # CIF has ~329 bonds from connectivity tables
        "template_bonds": (1480, 1540),  # Around 1500-1530 template bonds
        "distance_bonds": (620, 670),  # Around 648 distance bonds
        "total_bonds": (2450, 2500),  # Around 2477 total bonds
        "hydrogen_bonds_min": 650,  # At least 650 H bonds
        "hydrogen_fraction_min": 0.3,  # At least 30% H bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Ribonuclease A CIF format with chemical bond connectivity",
        "has_conect_records": True,  # CIF has extensive bond connectivity
        "has_hydrogen_atoms": True,
    },
    {
        "filename": "2izf.cif",
        "file_bonds": (240, 260),  # CIF has ~253 bonds from connectivity tables
        "template_bonds": (2000, 2070),  # Around 2026-2060 template bonds
        "distance_bonds": (800, 850),  # Around 822 distance bonds
        "total_bonds": (3080, 3120),  # Around 3101 total bonds
        "hydrogen_bonds_min": 850,  # At least 850 H bonds
        "hydrogen_fraction_min": 0.3,  # At least 30% H bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Streptavidin-biotin complex CIF format with chemical bond connectivity",
        "has_conect_records": True,  # CIF has extensive bond connectivity
        "has_hydrogen_atoms": True,
    },
    {
        "filename": "1crn.cif",
        # Small structure - crambin CIF with connectivity tables
        "file_bonds": (110, 130),  # CIF has ~120 bonds from connectivity tables
        "template_bonds": (270, 295),  # Around 276-289 template bonds
        "distance_bonds": (0, 10),  # No distance bonds needed (all covered)
        "total_bonds": (390, 410),  # Around 396 total bonds
        "hydrogen_bonds_min": None,  # No hydrogen atoms in this structure
        "hydrogen_fraction_min": None,  # No hydrogen bonds
        "min_atoms_per_second": 10,
        "max_execution_time": 30.0,
        "max_unreasonable_fraction": 0.1,
        "description": "Small crambin structure CIF format with chemical bond connectivity",
        "has_conect_records": True,  # CIF has bond connectivity tables
        "has_hydrogen_atoms": False,  # No hydrogen atoms in this file
    },
]


class TestEndToEndBondDetection:
    """Data-driven end-to-end tests for bond detection on real structures."""

    @pytest.fixture(params=TEST_DATA)
    def test_data(self, request):
        """Load test structure and expected values."""
        data = request.param
        filename = data["filename"]
        file_path = TEST_DATA_DIR / filename
        if not file_path.exists():
            pytest.skip(f"Test file {filename} not found")

        # Choose parser based on file extension
        if filename.endswith(".cif"):
            parser = mmCIFParser()
        else:
            parser = PDBParser()

        structure = parser.parse_file(str(file_path))
        return {"structure": structure, "expected": data, "filename": filename}

    def test_bond_detection_methods(self, test_data):
        """Test that different bond detection methods are being used."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure, use_file_bonds=True)

        # Check that detection_method attribute is populated
        assert hasattr(
            bonds, "detection_method"
        ), "BondList should have detection_method"
        assert (
            bonds.detection_method is not None
        ), "detection_method should be initialized"

        # Count bonds by detection method
        method_counts = {}
        for i in range(len(bonds)):
            method = bonds.detection_method[i]
            method_counts[method] = method_counts.get(method, 0) + 1

        print(f"{filename} bond detection methods: {method_counts}")

        # Validate file bonds if expected
        if expected["file_bonds"] is not None:
            file_count = method_counts.get("file", 0)
            min_file, max_file = expected["file_bonds"]
            assert (
                min_file <= file_count <= max_file
            ), f"{filename}: Expected {min_file}-{max_file} file bonds, got {file_count}"

        # Validate template bonds if expected
        if expected["template_bonds"] is not None:
            template_count = method_counts.get("template", 0)
            min_template, max_template = expected["template_bonds"]
            assert (
                min_template <= template_count <= max_template
            ), f"{filename}: Expected {min_template}-{max_template} template bonds, got {template_count}"

        # Validate distance bonds if expected
        if expected["distance_bonds"] is not None:
            distance_count = method_counts.get("distance", 0)
            min_distance, max_distance = expected["distance_bonds"]
            assert (
                min_distance <= distance_count <= max_distance
            ), f"{filename}: Expected {min_distance}-{max_distance} distance bonds, got {distance_count}"

        # Validate total bonds if expected
        if expected["total_bonds"] is not None:
            total_count = len(bonds)
            min_total, max_total = expected["total_bonds"]
            assert (
                min_total <= total_count <= max_total
            ), f"{filename}: Expected {min_total}-{max_total} total bonds, got {total_count}"

        # Should have multiple detection methods (unless it's a very small structure)
        if len(bonds) > 100:
            assert (
                len(method_counts) > 1
            ), f"{filename}: Should use multiple detection methods"

        # Should have file-based bonds if CONECT records expected
        if expected["has_conect_records"]:
            assert (
                "file" in method_counts
            ), f"{filename}: Should detect file-based bonds from CONECT records"

        # Should have template or distance bonds
        assert (
            "template" in method_counts or "distance" in method_counts
        ), f"{filename}: Should detect template or distance-based bonds"

    def test_conect_record_bonds(self, test_data):
        """Test that CONECT record bonds are properly parsed and used."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Check that structure has file bonds from CONECT records
        assert hasattr(
            structure, "file_bonds"
        ), f"{filename} should have file_bonds attribute"

        if expected["has_conect_records"]:
            if structure.file_bonds is not None:
                file_bond_count = len(structure.file_bonds)
                print(f"{filename}: {file_bond_count} file bonds from CONECT records")

                assert file_bond_count > 0, f"Should parse CONECT records in {filename}"

                # Test that file bonds are used in detection
                detector = DefaultBondDetector()
                bonds = detector.detect_bonds(structure, use_file_bonds=True)

                # Count file method bonds
                file_method_count = 0
                if (
                    hasattr(bonds, "detection_method")
                    and bonds.detection_method is not None
                ):
                    for i in range(len(bonds)):
                        if bonds.detection_method[i] == "file":
                            file_method_count += 1

                assert (
                    file_method_count > 0
                ), f"Should use file bonds in detection for {filename}"
                assert (
                    file_method_count <= file_bond_count
                ), "File method bonds <= parsed file bonds"

    def test_hydrogen_bond_detection(self, test_data):
        """Test that hydrogen bonds are properly detected."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        if not expected["has_hydrogen_atoms"]:
            pytest.skip(f"{filename} does not have hydrogen atoms")

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure, use_file_bonds=False)

        # Count hydrogen bonds
        h_bond_count = 0
        for i in range(len(bonds)):
            atom1, atom2 = bonds.bonds[i]
            elem1, elem2 = structure.element[atom1], structure.element[atom2]
            if elem1 == "H" or elem2 == "H":
                h_bond_count += 1

        print(f"{filename}: {h_bond_count} hydrogen bonds out of {len(bonds)} total")

        # Should detect minimum expected hydrogen bonds
        if expected["hydrogen_bonds_min"] is not None:
            assert (
                h_bond_count >= expected["hydrogen_bonds_min"]
            ), f"{filename}: Expected at least {expected['hydrogen_bonds_min']} H bonds, got {h_bond_count}"

        # Should meet minimum hydrogen fraction
        if expected["hydrogen_fraction_min"] is not None and len(bonds) > 0:
            h_fraction = h_bond_count / len(bonds)
            assert (
                h_fraction >= expected["hydrogen_fraction_min"]
            ), f"{filename}: Expected H bonds >= {expected['hydrogen_fraction_min']:.1%}, got {h_fraction:.1%}"

    def test_distance_vs_file_bonds(self, test_data):
        """Test that file bonds add to distance-based detection."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        detector = DefaultBondDetector()

        # Test with all methods enabled
        bonds_all = detector.detect_bonds(structure, use_file_bonds=True)

        # Test with only distance-based detection
        bonds_distance = detector.detect_bonds(structure, use_file_bonds=False)

        # Basic validation
        assert len(bonds_all) > 0, f"Should detect bonds in {filename}"
        assert (
            len(bonds_distance) > 0
        ), f"Should detect distance-based bonds in {filename}"
        assert isinstance(bonds_all, BondList)
        assert isinstance(bonds_distance, BondList)

        # File bonds should increase total bond count (if CONECT records exist)
        if expected["has_conect_records"]:
            assert len(bonds_all) >= len(
                bonds_distance
            ), f"File bonds should add to total in {filename}"

        print(
            f"{filename}: {len(bonds_all)} total bonds, {len(bonds_distance)} distance-only bonds"
        )

    def test_residue_template_bonds(self, test_data):
        """Test that template-based bonds are detected for standard residues."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure, use_file_bonds=False)

        # Count template bonds
        template_bond_count = 0
        for i in range(len(bonds)):
            if (
                hasattr(bonds, "detection_method")
                and bonds.detection_method is not None
            ):
                if bonds.detection_method[i] == "template":
                    template_bond_count += 1

        print(f"{filename}: {template_bond_count} template bonds")

        # Validate template bond count if expected
        if expected["template_bonds"] is not None:
            min_template, max_template = expected["template_bonds"]
            assert (
                min_template <= template_bond_count <= max_template
            ), f"{filename}: Expected {min_template}-{max_template} template bonds, got {template_bond_count}"

        # Should have some template bonds for standard protein residues
        # (Note: this may be 0 if residue templates aren't fully implemented yet)
        if template_bond_count > 0:
            # For small structures, all bonds might be from templates
            if len(bonds) > 500:  # Only check for larger structures
                assert template_bond_count < len(
                    bonds
                ), "Not all bonds should be from templates"

    def test_bond_validation(self, test_data):
        """Test that detected bonds are chemically reasonable."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure, use_file_bonds=True)

        # Validate bond indices
        for i in range(len(bonds)):
            atom1, atom2 = bonds.bonds[i]
            assert 0 <= atom1 < structure.n_atoms, f"Invalid atom1 index in {filename}"
            assert 0 <= atom2 < structure.n_atoms, f"Invalid atom2 index in {filename}"
            assert atom1 != atom2, f"Self-bond detected in {filename}"
            assert atom1 < atom2, f"Bond atoms not properly ordered in {filename}"

        # Validate bond distances are reasonable
        unreasonable_count = 0
        for i in range(len(bonds)):
            atom1, atom2 = bonds.bonds[i]
            distance = np.linalg.norm(structure.coord[atom1] - structure.coord[atom2])

            # Very loose bounds for reasonable bond lengths
            if distance < 0.5 or distance > 5.0:
                unreasonable_count += 1

        # Allow some unreasonable bonds (might be metal coordination, etc.)
        unreasonable_fraction = unreasonable_count / len(bonds) if len(bonds) > 0 else 0
        assert (
            unreasonable_fraction < expected["max_unreasonable_fraction"]
        ), f"Too many unreasonable bonds in {filename}: {unreasonable_fraction:.2%}"

        print(
            f"{filename}: {unreasonable_count}/{len(bonds)} potentially unreasonable bonds"
        )

    def test_structure_integration(self, test_data):
        """Test integration with Structure.detect_bonds() method."""
        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Test default behavior
        bonds1 = structure.detect_bonds()
        assert isinstance(bonds1, BondList)
        assert len(bonds1) > 0

        # Test parameter variations
        bonds2 = structure.detect_bonds(vdw_factor=0.5, use_file_bonds=False)
        bonds3 = structure.detect_bonds(vdw_factor=0.9, use_file_bonds=True)

        assert isinstance(bonds2, BondList)
        assert isinstance(bonds3, BondList)

        # Different parameters should potentially give different results
        print(
            f"{filename} Structure.detect_bonds(): default={len(bonds1)}, "
            f"vdw=0.5/no-file={len(bonds2)}, vdw=0.9/with-file={len(bonds3)}"
        )

        # Test that bonds are stored in structure
        assert (
            structure.bonds is not None
        ), f"Bonds should be stored in {filename} structure"
        # Note: structure.bonds contains the last detection result (bonds3)
        assert len(structure.bonds) == len(
            bonds3
        ), f"Stored bonds count mismatch in {filename}"

    def test_performance_baseline(self, test_data):
        """Basic performance test to ensure reasonable execution time."""
        import time

        structure = test_data["structure"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        detector = DefaultBondDetector()

        start_time = time.time()
        bonds = detector.detect_bonds(structure, use_file_bonds=True)
        end_time = time.time()

        execution_time = end_time - start_time
        atoms_per_second = (
            structure.n_atoms / execution_time if execution_time > 0 else float("inf")
        )

        print(
            f"{filename}: {structure.n_atoms} atoms, {len(bonds)} bonds, "
            f"{execution_time:.3f}s ({atoms_per_second:.0f} atoms/s)"
        )

        # Validate performance expectations
        assert (
            execution_time < expected["max_execution_time"]
        ), f"Bond detection too slow for {filename}: {execution_time:.3f}s"
        assert (
            atoms_per_second > expected["min_atoms_per_second"]
        ), f"Bond detection throughput too low for {filename}: {atoms_per_second:.0f} atoms/s"


class TestBondDetectionEdgeCases:
    """Test edge cases and error handling."""

    def test_structure_without_file_bonds(self):
        """Test bond detection on structure without CONECT records."""
        # Create minimal structure without file bonds
        structure = Structure(3)
        structure.coord[0] = [0.0, 0.0, 0.0]
        structure.coord[1] = [1.5, 0.0, 0.0]
        structure.coord[2] = [0.0, 1.5, 0.0]

        structure.atom_name[:] = ["C1", "C2", "O1"]
        structure.element[:] = ["C", "C", "O"]
        structure.res_name[:] = "UNK"
        structure.res_id[:] = 1
        structure.chain_id[:] = "A"

        # Ensure no file bonds
        assert structure.file_bonds is None

        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(structure, use_file_bonds=True)

        # Should still work without file bonds
        assert isinstance(bonds, BondList)
        # May or may not have bonds depending on distance criteria

    def test_different_vdw_factors(self):
        """Test bond detection with different VdW factors on 6rsa.pdb."""
        pdb_file = TEST_DATA_DIR / "6rsa.pdb"
        parser = PDBParser()
        structure = parser.parse_file(str(pdb_file))

        detector_tight = DefaultBondDetector(vdw_factor=0.5)
        detector_loose = DefaultBondDetector(vdw_factor=0.9)

        bonds_tight = detector_tight.detect_bonds(structure, use_file_bonds=False)
        bonds_loose = detector_loose.detect_bonds(structure, use_file_bonds=False)

        # Loose criteria should generally find more or equal bonds
        assert len(bonds_loose) >= len(
            bonds_tight
        ), "Loose VdW factor should find >= bonds"

        print(
            f"VdW factor comparison: tight(0.5)={len(bonds_tight)}, loose(0.9)={len(bonds_loose)}"
        )


if __name__ == "__main__":
    # Allow running the test directly for development
    pytest.main([__file__, "-v", "-s"])
