"""
StructureEnsemble class for trajectory representation.

This module provides the StructureEnsemble class for handling molecular trajectories
with multiple frames, using memory-efficient storage where annotations are shared
across frames and only coordinates vary.
"""

from typing import Any, List, Optional, Union

import numpy as np

from .structure import Structure


class StructureEnsemble:
    """
    Ensemble of molecular structures representing trajectory data.

    This class stores multiple frames of coordinate data while sharing
    annotations (atom names, elements, etc.) across all frames for
    memory efficiency. Designed for trajectory analysis and multi-model
    PDB files.

    Memory layout:
        - coords: (n_frames, n_atoms, 3) array of coordinates
        - Annotations shared from template Structure
        - Optional time and box information per frame

    Example:
        >>> ensemble = StructureEnsemble.from_structures([struct1, struct2])
        >>> print(f"Trajectory with {ensemble.n_frames} frames")
        >>> frame0 = ensemble[0]  # Returns Structure for frame 0
    """

    def __init__(self, template: Structure, n_frames: int = 0):
        """
        Initialize StructureEnsemble from template Structure.

        Args:
            template: Template Structure with shared annotations
            n_frames: Number of frames (default: 0 for dynamic growth)
        """
        self.template = template.copy()  # Store copy of template
        self.n_atoms = template.n_atoms
        self.n_frames = n_frames
        self._capacity = max(n_frames, 10)  # Minimum capacity

        # Coordinate storage (n_frames, n_atoms, 3)
        self.coords = np.zeros((self._capacity, self.n_atoms, 3), dtype=np.float64)

        # Optional per-frame data
        self.times: Optional[np.ndarray] = None
        self.box_vectors: Optional[np.ndarray] = None  # (n_frames, 3, 3) for periodic

    @classmethod
    def from_pdb(cls, filename: str) -> "StructureEnsemble":
        """
        Create StructureEnsemble from multi-model PDB file.

        Args:
            filename: Path to multi-model PDB file

        Returns:
            StructureEnsemble with all models as frames

        Raises:
            ValueError: If PDB file contains only single model
        """
        from ..io.pdb import PDBParser

        parser = PDBParser()
        result = parser.parse_file(filename)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "PDB file contains only single model, use Structure.from_pdb() instead"
            )

    @classmethod
    def from_pdb_string(cls, pdb_content: str) -> "StructureEnsemble":
        """
        Create StructureEnsemble from multi-model PDB content string.

        Args:
            pdb_content: PDB content string with multiple models

        Returns:
            StructureEnsemble with all models as frames

        Raises:
            ValueError: If PDB content contains only single model
        """
        from ..io.pdb import PDBParser

        parser = PDBParser()
        result = parser.parse_string(pdb_content)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "PDB content contains only single model, use Structure.from_pdb_string() instead"
            )

    @classmethod
    def from_mmcif(cls, filename: str) -> "StructureEnsemble":
        """
        Create StructureEnsemble from multi-model mmCIF file.

        Args:
            filename: Path to multi-model mmCIF file

        Returns:
            StructureEnsemble with all models as frames

        Raises:
            ValueError: If mmCIF file contains only single model
        """
        from ..io.mmcif import mmCIFParser

        parser = mmCIFParser()
        result = parser.parse_file(filename)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "mmCIF file contains only single model, use Structure.from_mmcif() instead"
            )

    @classmethod
    def from_mmcif_string(cls, mmcif_content: str) -> "StructureEnsemble":
        """
        Create StructureEnsemble from multi-model mmCIF content string.

        Args:
            mmcif_content: mmCIF content string with multiple models

        Returns:
            StructureEnsemble with all models as frames

        Raises:
            ValueError: If mmCIF content contains only single model
        """
        from ..io.mmcif import mmCIFParser

        parser = mmCIFParser()
        result = parser.parse_string(mmcif_content)

        if isinstance(result, cls):
            return result
        else:
            raise ValueError(
                "mmCIF content contains only single model, use Structure.from_mmcif_string() instead"
            )

    @classmethod
    def from_structures(cls, structures: List[Structure]) -> "StructureEnsemble":
        """
        Create StructureEnsemble from list of Structure objects.

        Args:
            structures: List of Structure objects with same atoms

        Returns:
            StructureEnsemble with structures as frames

        Raises:
            ValueError: If structures have different atom counts
        """
        if not structures:
            raise ValueError("Cannot create ensemble from empty structure list")

        template = structures[0]
        n_frames = len(structures)

        # Validate all structures have same atom count
        for i, struct in enumerate(structures):
            if struct.n_atoms != template.n_atoms:
                raise ValueError(
                    f"Structure {i} has {struct.n_atoms} atoms, expected {template.n_atoms}"
                )

        ensemble = cls(template, n_frames)

        # Copy coordinates from each structure
        for i, struct in enumerate(structures):
            ensemble.coords[i] = struct.coord

        return ensemble

    def _ensure_capacity(self, required_frames: int) -> None:
        """
        Ensure sufficient capacity, growing arrays if necessary.

        Args:
            required_frames: Required minimum number of frames
        """
        if required_frames <= self._capacity:
            return

        new_capacity = max(int(self._capacity * 1.5), required_frames)

        # Resize coordinate array
        new_coords = np.zeros((new_capacity, self.n_atoms, 3), dtype=np.float64)
        new_coords[: self._capacity] = self.coords
        self.coords = new_coords

        # Resize optional arrays
        if self.times is not None:
            new_times = np.zeros(new_capacity, dtype=np.float64)
            new_times[: self._capacity] = self.times
            self.times = new_times

        if self.box_vectors is not None:
            new_box = np.zeros((new_capacity, 3, 3), dtype=np.float64)
            new_box[: self._capacity] = self.box_vectors
            self.box_vectors = new_box

        self._capacity = new_capacity

    def add_frame(self, structure: Structure, time: Optional[float] = None) -> None:
        """
        Add a new frame to the ensemble.

        Args:
            structure: Structure to add as new frame
            time: Optional time value for this frame

        Raises:
            ValueError: If structure atom count doesn't match
        """
        if structure.n_atoms != self.n_atoms:
            raise ValueError(
                f"Structure has {structure.n_atoms} atoms, expected {self.n_atoms}"
            )

        self._ensure_capacity(self.n_frames + 1)

        # Add coordinates
        self.coords[self.n_frames] = structure.coord

        # Add time if provided
        if time is not None:
            if self.times is None:
                self.times = np.full(self._capacity, np.nan, dtype=np.float64)
            self.times[self.n_frames] = time

        self.n_frames += 1

    def __getitem__(
        self, index: Union[int, slice]
    ) -> Union[Structure, "StructureEnsemble"]:
        """
        Get frame(s) from ensemble.

        Args:
            index: Frame index or slice

        Returns:
            Structure for single frame, StructureEnsemble for slice

        Examples:
            >>> frame0 = ensemble[0]  # Single frame as Structure
            >>> sub_traj = ensemble[10:20]  # Sub-trajectory as StructureEnsemble
        """
        if isinstance(index, int):
            # Return single frame as Structure
            if index < 0:
                index = self.n_frames + index
            if not (0 <= index < self.n_frames):
                raise IndexError(f"Frame index {index} out of range")

            frame = self.template.copy()
            frame.coord = self.coords[index].copy()
            return frame

        elif isinstance(index, slice):
            # Return sub-trajectory as StructureEnsemble
            start, stop, step = index.indices(self.n_frames)
            frame_indices = range(start, stop, step)

            sub_ensemble = StructureEnsemble(self.template, len(frame_indices))
            for i, frame_idx in enumerate(frame_indices):
                sub_ensemble.coords[i] = self.coords[frame_idx]

            # Copy time data if available
            if self.times is not None:
                sub_ensemble.times = np.zeros(sub_ensemble._capacity, dtype=np.float64)
                for i, frame_idx in enumerate(frame_indices):
                    sub_ensemble.times[i] = self.times[frame_idx]

            return sub_ensemble

        else:
            raise TypeError("Index must be int or slice")

    def __len__(self) -> int:
        """Return number of frames."""
        return self.n_frames

    def __iter__(self) -> Any:
        """Iterate over frames as Structure objects."""
        for i in range(self.n_frames):
            yield self[i]

    def get_frame_coords(self, frame_index: int) -> np.ndarray:
        """
        Get coordinates for specific frame.

        Args:
            frame_index: Index of frame

        Returns:
            Coordinate array (n_atoms, 3) for the frame
        """
        if not (0 <= frame_index < self.n_frames):
            raise IndexError(f"Frame index {frame_index} out of range")
        return self.coords[frame_index].copy()  # type: ignore[no-any-return]

    def set_frame_coords(self, frame_index: int, coords: np.ndarray) -> None:
        """
        Set coordinates for specific frame.

        Args:
            frame_index: Index of frame
            coords: Coordinate array (n_atoms, 3)
        """
        if not (0 <= frame_index < self.n_frames):
            raise IndexError(f"Frame index {frame_index} out of range")
        if coords.shape != (self.n_atoms, 3):
            raise ValueError(
                f"Coordinates shape {coords.shape} doesn't match expected {(self.n_atoms, 3)}"
            )

        self.coords[frame_index] = coords

    def center_frames(self, selection: Optional[str] = None) -> None:
        """
        Center all frames at origin.

        Args:
            selection: Optional selection for center calculation (default: all atoms)
        """
        if selection is None:
            # Center on all atoms
            for i in range(self.n_frames):
                center = np.mean(self.coords[i], axis=0)
                self.coords[i] -= center
        else:
            # Center on selected atoms
            mask = self.template.select(selection)
            for i in range(self.n_frames):
                center = np.mean(self.coords[i][mask], axis=0)
                self.coords[i] -= center

    def rmsd(
        self, reference_frame: int = 0, selection: Optional[str] = None
    ) -> np.ndarray:
        """
        Calculate RMSD of each frame relative to reference.

        Args:
            reference_frame: Index of reference frame
            selection: Optional atom selection for RMSD calculation

        Returns:
            Array of RMSD values for each frame
        """
        if not (0 <= reference_frame < self.n_frames):
            raise IndexError(f"Reference frame {reference_frame} out of range")

        ref_coords = self.coords[reference_frame]
        rmsd_values = np.zeros(self.n_frames)

        if selection is None:
            # Use all atoms
            for i in range(self.n_frames):
                diff = self.coords[i] - ref_coords
                rmsd_values[i] = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        else:
            # Use selected atoms
            mask = self.template.select(selection)
            ref_coords_sel = ref_coords[mask]
            for i in range(self.n_frames):
                coords_sel = self.coords[i][mask]
                diff = coords_sel - ref_coords_sel
                rmsd_values[i] = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        return rmsd_values  # type: ignore[no-any-return]

    def __repr__(self) -> str:
        """String representation."""
        return f"StructureEnsemble(n_atoms={self.n_atoms}, n_frames={self.n_frames})"

    def __str__(self) -> str:
        """Detailed string representation."""
        lines = [
            f"StructureEnsemble with {self.n_frames} frames and {self.n_atoms} atoms"
        ]

        if self.n_frames > 0:
            lines.append(f"Frame indices: 0 to {self.n_frames - 1}")

        if self.times is not None:
            time_range = f"{self.times[0]:.2f} to {self.times[self.n_frames-1]:.2f}"
            lines.append(f"Time range: {time_range}")

        return "\n".join(lines)
