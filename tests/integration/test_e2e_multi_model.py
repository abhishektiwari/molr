"""
End-to-end integration tests for multi-model PDB structure parsing.

This module tests the complete parsing pipeline for multi-model PDB files
to validate that StructureEnsemble objects are created correctly and that
individual frames can be accessed and processed.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pytest

from molr.bond_detection.default_detector import DefaultBondDetector
from molr.core.bond_list import BondList
from molr.core.structure import Structure
from molr.core.structure_ensemble import StructureEnsemble
from molr.io.pdb import PDBParser

# Get path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "example_pdb_files"


# Test data for multi-model structures
MULTI_MODEL_TEST_DATA = [
    {
        "filename": "generated/multi_model.pdb",
        "n_frames": 2,
        "n_atoms": 6,
        "elements": {"C", "N", "O"},
        "residue_types": {"GLY", "HOH"},
        "has_hydrogen": False,
        "frame_variations": True,  # Coordinates vary between frames
        "description": "Simple 2-frame structure with GLY and HOH",
        # Bond detection expectations for frame 0
        "frame0_total_bonds": (3, 10),
        "frame0_template_bonds": (3, 8),
        "frame0_distance_bonds": (0, 3),
        # Coordinate differences between frames
        "max_coord_diff": 1.5,  # Maximum expected coordinate difference (1.4 observed)
        "min_coord_diff": 0.0,  # Some atoms might not move
    },
    {
        "filename": "generated/test_multi_model_ags.pdb",
        "n_frames": 3,
        "n_atoms": 16,
        "elements": {"C", "N", "O"},
        "residue_types": {"ALA", "GLY", "SER"},
        "has_hydrogen": False,
        "frame_variations": True,  # Coordinates vary between frames
        "description": "3-frame ALA-GLY-SER tripeptide with conformational changes",
        # Bond detection expectations for frame 0
        "frame0_total_bonds": (12, 15),
        "frame0_template_bonds": (12, 15),
        "frame0_distance_bonds": (0, 2),
        # Coordinate differences between frames
        "max_coord_diff": 2.5,  # Maximum expected coordinate difference (2.465 observed)
        "min_coord_diff": 0.0,  # Some atoms might not move
    },
]


class TestMultiModelParsing:
    """Test parsing and handling of multi-model PDB files."""

    @pytest.fixture(params=MULTI_MODEL_TEST_DATA)
    def test_data(self, request):
        """Load multi-model structure and expected values."""
        data = request.param
        filename = data["filename"]
        file_path = TEST_DATA_DIR / filename
        if not file_path.exists():
            pytest.skip(f"Test file {filename} not found")

        parser = PDBParser()
        result = parser.parse_file(str(file_path))

        # Verify we got a StructureEnsemble
        if not isinstance(result, StructureEnsemble):
            pytest.fail(f"{filename} did not parse as StructureEnsemble")

        return {"ensemble": result, "expected": data, "filename": filename}

    def test_ensemble_properties(self, test_data):
        """Test basic properties of StructureEnsemble."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Check frame count
        assert (
            ensemble.n_frames == expected["n_frames"]
        ), f"{filename}: Expected {expected['n_frames']} frames, got {ensemble.n_frames}"

        # Check atom count per frame
        assert (
            ensemble.n_atoms == expected["n_atoms"]
        ), f"{filename}: Expected {expected['n_atoms']} atoms per frame, got {ensemble.n_atoms}"

        # Check that template structure exists
        assert hasattr(ensemble, "template"), f"{filename}: Missing template structure"
        assert isinstance(
            ensemble.template, Structure
        ), f"{filename}: Template should be a Structure object"

        # Verify coordinate array shape
        assert (
            ensemble.coords.shape[0] >= ensemble.n_frames
        ), f"{filename}: Coordinate array too small"
        assert (
            ensemble.coords.shape[1] == ensemble.n_atoms
        ), f"{filename}: Wrong number of atoms in coordinate array"
        assert ensemble.coords.shape[2] == 3, f"{filename}: Coordinates should be 3D"

        print(f"{filename}: {ensemble.n_frames} frames, {ensemble.n_atoms} atoms each")

    def test_frame_access(self, test_data):
        """Test accessing individual frames from ensemble."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Test accessing each frame
        for i in range(ensemble.n_frames):
            frame = ensemble[i]
            assert isinstance(
                frame, Structure
            ), f"{filename}: Frame {i} should be a Structure object"
            assert (
                frame.n_atoms == expected["n_atoms"]
            ), f"{filename}: Frame {i} has wrong atom count"

            # Check elements
            elements = set(frame.element)
            assert (
                elements == expected["elements"]
            ), f"{filename}: Frame {i} has wrong elements: {elements}"

            # Check residue types
            res_types = set(frame.res_name)
            assert (
                res_types == expected["residue_types"]
            ), f"{filename}: Frame {i} has wrong residue types: {res_types}"

            # Check hydrogen presence
            has_h = "H" in frame.element
            assert (
                has_h == expected["has_hydrogen"]
            ), f"{filename}: Frame {i} hydrogen presence mismatch"

        # Test negative indexing
        last_frame = ensemble[-1]
        assert isinstance(
            last_frame, Structure
        ), f"{filename}: Negative indexing should work"
        assert (
            last_frame.n_atoms == expected["n_atoms"]
        ), f"{filename}: Last frame has wrong atom count"

        # Test out of bounds access
        with pytest.raises(IndexError):
            _ = ensemble[ensemble.n_frames]

    def test_frame_variations(self, test_data):
        """Test that coordinates vary between frames."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        if not expected.get("frame_variations", True):
            pytest.skip(f"{filename} does not have frame variations")

        if ensemble.n_frames < 2:
            pytest.skip(f"{filename} has only one frame")

        # Compare coordinates between first and last frame
        frame0 = ensemble[0]
        frame_last = ensemble[-1]

        # Calculate coordinate differences
        coord_diff = np.abs(frame_last.coord - frame0.coord)
        max_diff = np.max(coord_diff)
        min_diff = np.min(coord_diff)
        mean_diff = np.mean(coord_diff)

        print(
            f"{filename} coordinate variations: "
            f"max={max_diff:.3f}, min={min_diff:.3f}, mean={mean_diff:.3f}"
        )

        # Some coordinates should change
        assert max_diff > 0.01, f"{filename}: No coordinate changes between frames"

        # But not too much for a reasonable structure
        if "max_coord_diff" in expected:
            assert (
                max_diff <= expected["max_coord_diff"]
            ), f"{filename}: Coordinate changes too large: {max_diff}"

    def test_bond_detection_on_frames(self, test_data):
        """Test bond detection on individual frames."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Test bond detection on first frame
        frame0 = ensemble[0]
        detector = DefaultBondDetector()
        bonds = detector.detect_bonds(frame0)

        assert isinstance(
            bonds, BondList
        ), f"{filename}: Bond detection should return BondList"

        # Check total bond count
        if "frame0_total_bonds" in expected:
            min_bonds, max_bonds = expected["frame0_total_bonds"]
            assert (
                min_bonds <= len(bonds) <= max_bonds
            ), f"{filename}: Frame 0 expected {min_bonds}-{max_bonds} bonds, got {len(bonds)}"

        # Check bond method distribution
        if hasattr(bonds, "detection_method") and bonds.detection_method is not None:
            method_counts = {}
            for i in range(len(bonds)):
                method = bonds.detection_method[i]
                method_counts[method] = method_counts.get(method, 0) + 1

            print(f"{filename} frame 0 bond methods: {method_counts}")

            # Validate template bonds if expected
            if "frame0_template_bonds" in expected:
                template_count = method_counts.get("template", 0)
                min_template, max_template = expected["frame0_template_bonds"]
                assert (
                    min_template <= template_count <= max_template
                ), f"{filename}: Frame 0 expected {min_template}-{max_template} template bonds, got {template_count}"

    def test_ensemble_slicing(self, test_data):
        """Test slicing operations on StructureEnsemble."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        if ensemble.n_frames < 2:
            pytest.skip(f"{filename} needs at least 2 frames for slicing")

        # Test simple slice
        sub_ensemble = ensemble[0:1]
        assert isinstance(
            sub_ensemble, StructureEnsemble
        ), f"{filename}: Slicing should return StructureEnsemble"
        assert (
            sub_ensemble.n_frames == 1
        ), f"{filename}: Slice [0:1] should have 1 frame"
        assert (
            sub_ensemble.n_atoms == ensemble.n_atoms
        ), f"{filename}: Sliced ensemble should have same atom count"

        # Test step slice
        if ensemble.n_frames >= 2:
            step_ensemble = ensemble[::2]
            expected_frames = len(range(0, ensemble.n_frames, 2))
            assert (
                step_ensemble.n_frames == expected_frames
            ), f"{filename}: Step slice should have {expected_frames} frames"

        # Test negative slice
        last_two = ensemble[-2:]
        expected_last = min(2, ensemble.n_frames)
        assert (
            last_two.n_frames == expected_last
        ), f"{filename}: Last two frames slice should have {expected_last} frames"

    def test_structure_detect_bonds_integration(self, test_data):
        """Test Structure.detect_bonds() method on frames."""
        ensemble = test_data["ensemble"]
        expected = test_data["expected"]
        filename = test_data["filename"]

        # Test detect_bonds on each frame
        bond_counts = []
        for i in range(min(3, ensemble.n_frames)):  # Test up to 3 frames
            frame = ensemble[i]
            bonds = frame.detect_bonds()

            assert isinstance(
                bonds, BondList
            ), f"{filename}: Frame {i} detect_bonds should return BondList"
            assert (
                frame.bonds is not None
            ), f"{filename}: Frame {i} should store bonds after detection"
            assert len(frame.bonds) == len(
                bonds
            ), f"{filename}: Frame {i} stored bonds mismatch"

            bond_counts.append(len(bonds))

        print(f"{filename} bond counts across frames: {bond_counts}")

        # For structures with variations, bond counts might differ slightly
        if expected.get("frame_variations", True) and len(bond_counts) > 1:
            # But shouldn't vary too much
            max_variation = max(bond_counts) - min(bond_counts)
            assert (
                max_variation <= 5
            ), f"{filename}: Bond count varies too much across frames: {bond_counts}"


class TestMultiModelEdgeCases:
    """Test edge cases for multi-model structure handling."""

    def test_single_frame_ensemble(self):
        """Test creating ensemble with single frame."""
        # Create a simple structure
        structure = Structure(3)
        structure.coord[0] = [0.0, 0.0, 0.0]
        structure.coord[1] = [1.5, 0.0, 0.0]
        structure.coord[2] = [0.0, 1.5, 0.0]
        structure.atom_name[:] = ["C1", "C2", "O1"]
        structure.element[:] = ["C", "C", "O"]
        structure.res_name[:] = "UNK"
        structure.res_id[:] = 1
        structure.chain_id[:] = "A"

        # Create ensemble with single frame
        ensemble = StructureEnsemble(structure, n_frames=1)
        ensemble.coords[0] = structure.coord
        ensemble.n_frames = 1

        # Test access
        frame0 = ensemble[0]
        assert isinstance(frame0, Structure)
        assert frame0.n_atoms == 3
        assert np.allclose(frame0.coord, structure.coord)

        # Test slicing
        sub = ensemble[0:1]
        assert isinstance(sub, StructureEnsemble)
        assert sub.n_frames == 1

    def test_coordinate_consistency(self):
        """Test that frame coordinates are independent."""
        # Create ensemble
        structure = Structure(2)
        structure.atom_name[:] = ["C1", "C2"]
        structure.element[:] = ["C", "C"]
        structure.res_name[:] = "UNK"
        structure.res_id[:] = 1
        structure.chain_id[:] = "A"

        ensemble = StructureEnsemble(structure, n_frames=2)
        ensemble.coords[0] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        ensemble.coords[1] = [[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        ensemble.n_frames = 2

        # Get frames
        frame0 = ensemble[0]
        frame1 = ensemble[1]

        # Modify frame0 coordinates
        frame0.coord[0] = [5.0, 5.0, 5.0]

        # Frame1 and ensemble should be unchanged
        assert frame1.coord[0, 0] == 0.0
        assert ensemble.coords[0, 0, 0] == 0.0
        assert ensemble.coords[1, 0, 0] == 0.0


class TestProblematicMultiModelFiles:
    """Test handling of problematic multi-model files."""

    def test_inconsistent_atom_count_file(self):
        """Test that parser properly handles files with inconsistent atom counts."""
        # This tests the known issue with 1bq0.pdb
        parser = PDBParser()
        pdb_file = TEST_DATA_DIR / "1bq0.pdb"

        if pdb_file.exists():
            with pytest.raises(IOError) as exc_info:
                parser.parse_file(str(pdb_file))

            assert "Model" in str(exc_info.value)
            assert "atoms" in str(exc_info.value)
            print(f"1bq0.pdb parsing error (expected): {exc_info.value}")
        else:
            pytest.skip("1bq0.pdb not found")


if __name__ == "__main__":
    # Allow running the test directly for development
    pytest.main([__file__, "-v", "-s"])
