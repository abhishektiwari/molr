"""
PDB file parser for the space module using pdbreader.

This module provides a clean PDB parser specifically designed for the NumPy-based
Structure class, converting pdbreader output directly to NumPy arrays for optimal
performance.
"""

import math
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from ..constants.atomic_data import AtomicData
from ..constants.pdb_constants import (
    DNA_RESIDUES,
    PROTEIN_RESIDUES,
    RNA_RESIDUES,
    WATER_MOLECULES,
)
from ..core.bond_list import BondList, BondOrder
from ..core.structure import Structure
from ..core.structure_ensemble import StructureEnsemble
from ..utilities import pdb_atom_to_element

try:
    import pdbreader
except ImportError:
    raise ImportError(
        "pdbreader package is required for PDB parsing. Install with: pip install pdbreader"
    )


def _safe_convert_int(value: Any, default: int = 0) -> int:
    """Safely convert value to integer, handling NaN and None."""
    if value is None:
        return default
    try:
        if isinstance(value, float) and math.isnan(value):
            return default
        return int(value)
    except (ValueError, TypeError):
        return default


def _safe_convert_float(value: Any, default: float = 0.0) -> float:
    """Safely convert value to float, handling NaN and None."""
    if value is None:
        return default
    try:
        float_val = float(value)
        if math.isnan(float_val):
            return default
        return float_val
    except (ValueError, TypeError):
        return default


def _safe_convert_str(value: Any, default: str = "") -> str:
    """Safely convert value to string, handling None."""
    if value is None:
        return default
    return str(value).strip()


class PDBParser:
    """
    PDB file parser for the space module.

    Designed specifically for the NumPy-based Structure class, this parser
    converts pdbreader output directly to NumPy arrays for optimal performance.

    Features:
    - Direct conversion to NumPy arrays
    - CONECT record parsing for explicit bonds
    - Multi-model support for trajectories
    - Efficient memory usage
    - Full PDB annotation support
    """

    def __init__(self) -> None:
        """Initialize the PDB parser."""
        self._current_structure: Optional[Structure] = None
        self._current_bonds: Optional[BondList] = None

    def parse_file(self, filename: str) -> Union[Structure, StructureEnsemble]:
        """
        Parse a PDB file and return a Structure.

        Args:
            filename: Path to the PDB file

        Returns:
            Structure object with all atoms and annotations

        Raises:
            IOError: If file cannot be read
            ValueError: If PDB format is invalid
        """
        try:
            # Use pdbreader to parse the file
            pdb_data = pdbreader.read_pdb(filename)
            return self._convert_pdb_data(pdb_data)
        except Exception as e:
            raise IOError(f"Error parsing PDB file '{filename}': {e}")

    def parse_string(self, pdb_content: str) -> Union[Structure, StructureEnsemble]:
        """
        Parse PDB content from a string.

        Args:
            pdb_content: PDB file content as string

        Returns:
            Structure object with all atoms and annotations
        """
        import os
        import tempfile

        try:
            # Write to temporary file since pdbreader expects file path
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".pdb", delete=False
            ) as f:
                f.write(pdb_content)
                temp_filename = f.name

            try:
                pdb_data = pdbreader.read_pdb(temp_filename)
                return self._convert_pdb_data(pdb_data)
            finally:
                os.unlink(temp_filename)

        except Exception as e:
            raise ValueError(f"Error parsing PDB content: {e}")

    def _convert_pdb_data(
        self, pdb_data: Dict[str, Any]
    ) -> Union[Structure, StructureEnsemble]:
        """
        Convert pdbreader output to Structure or StructureEnsemble.

        Args:
            pdb_data: Dictionary from pdbreader.read_pdb()

        Returns:
            Structure for single model, StructureEnsemble for multi-model
        """
        # Check if we have multiple models
        models = self._group_by_models(pdb_data)

        if len(models) == 1:
            # Single model - return Structure
            return self._create_single_structure(models[0], pdb_data)
        else:
            # Multiple models - return StructureEnsemble
            return self._create_structure_ensemble(models, pdb_data)

    def _group_by_models(self, pdb_data: Dict[str, Any]) -> List[List[Dict[str, Any]]]:
        """
        Group atom records by model number.

        Args:
            pdb_data: Dictionary from pdbreader.read_pdb()

        Returns:
            List of models, each containing atom records for that model
        """
        models: Dict[int, List] = {}  # model_id -> atom_records

        # Process ATOM records
        if "ATOM" in pdb_data and len(pdb_data["ATOM"]) > 0:
            for _, row in pdb_data["ATOM"].iterrows():
                model_id = row.get("model_id", 0)  # Default to model 0
                if model_id not in models:
                    models[model_id] = []
                models[model_id].append(self._convert_atom_record(row, is_hetero=False))

        # Process HETATM records
        if "HETATM" in pdb_data and len(pdb_data["HETATM"]) > 0:
            for _, row in pdb_data["HETATM"].iterrows():
                model_id = row.get("model_id", 0)  # Default to model 0
                if model_id not in models:
                    models[model_id] = []
                models[model_id].append(self._convert_atom_record(row, is_hetero=True))

        if not models:
            raise ValueError("No atom records found in PDB file")

        # Return models in order
        model_ids = sorted(models.keys())
        return [models[model_id] for model_id in model_ids]

    def _create_single_structure(
        self, atom_records: List[Dict[str, Any]], pdb_data: Dict[str, Any]
    ) -> Structure:
        """
        Create single Structure from atom records.

        Args:
            atom_records: List of atom record dictionaries
            pdb_data: Original PDB data for CONECT records

        Returns:
            Structure object
        """
        structure = self._create_structure_from_records(atom_records)

        # Parse CONECT records if available
        if "CONECT" in pdb_data and len(pdb_data["CONECT"]) > 0:
            bonds = self._parse_conect_records(pdb_data["CONECT"], structure)
            structure.file_bonds = bonds  # Store as file bonds

        return structure

    def _create_structure_ensemble(
        self, models: List[List[Dict[str, Any]]], pdb_data: Dict[str, Any]
    ) -> StructureEnsemble:
        """
        Create StructureEnsemble from multiple models.

        Args:
            models: List of models, each containing atom records
            pdb_data: Original PDB data for CONECT records

        Returns:
            StructureEnsemble object
        """
        if not models:
            raise ValueError("No models found")

        # Create template structure from first model
        template = self._create_structure_from_records(models[0])

        # Validate all models have same number of atoms
        n_atoms = len(models[0])
        for i, model in enumerate(models):
            if len(model) != n_atoms:
                raise ValueError(
                    f"Model {i} has {len(model)} atoms, expected {n_atoms}"
                )

        # Create ensemble
        ensemble = StructureEnsemble(template, len(models))

        # Add coordinates from each model
        for i, model in enumerate(models):
            coords = np.array([record["coord"] for record in model])
            ensemble.coords[i] = coords

        # Parse CONECT records if available (applied to template)
        if "CONECT" in pdb_data and len(pdb_data["CONECT"]) > 0:
            bonds = self._parse_conect_records(pdb_data["CONECT"], template)
            template.file_bonds = bonds  # Store as file bonds

        return ensemble

    def _convert_atom_record(self, row: Any, is_hetero: bool = False) -> Dict[str, Any]:
        """
        Convert a single atom record from pdbreader to dictionary.

        Args:
            row: Pandas row from pdbreader output
            is_hetero: Whether this is a HETATM record

        Returns:
            Dictionary with atom data
        """
        # Extract basic information using pdbreader column names
        serial = _safe_convert_int(row.get("id", 0))
        atom_name = _safe_convert_str(row.get("name", "")).upper()
        alt_loc = _safe_convert_str(row.get("loc_indicator", ""))
        res_name = _safe_convert_str(row.get("resname", "")).upper()
        chain_id = _safe_convert_str(row.get("chain", "A"))
        res_id = _safe_convert_int(row.get("resid", 0))
        insertion_code = _safe_convert_str(row.get("res_icode", ""))

        # Extract coordinates
        x = _safe_convert_float(row.get("x", 0.0))
        y = _safe_convert_float(row.get("y", 0.0))
        z = _safe_convert_float(row.get("z", 0.0))

        # Extract optional fields
        occupancy = _safe_convert_float(row.get("occupancy", 1.0))
        b_factor = _safe_convert_float(row.get("b_factor", 0.0))
        element = _safe_convert_str(row.get("element", ""))
        charge = _safe_convert_str(row.get("charge", ""))

        # Determine element if not provided
        if not element:
            element = pdb_atom_to_element(atom_name)

        return {
            "serial": serial,
            "atom_name": atom_name,
            "element": element,
            "res_name": res_name,
            "res_id": res_id,
            "chain_id": chain_id,
            "coord": np.array([x, y, z], dtype=np.float64),
            "alt_loc": alt_loc,
            "insertion_code": insertion_code,
            "occupancy": occupancy,
            "b_factor": b_factor,
            "charge": charge,
            "is_hetero": is_hetero,
        }

    def _create_structure_from_records(
        self, records: List[Dict[str, Any]]
    ) -> Structure:
        """
        Create Structure object from atom records.

        Args:
            records: List of atom record dictionaries

        Returns:
            Structure object
        """
        n_atoms = len(records)
        structure = Structure(n_atoms)

        # Fill core annotations
        for i, record in enumerate(records):
            structure.coord[i] = record["coord"]
            structure.atom_name[i] = record["atom_name"]
            structure.element[i] = record["element"]
            structure.res_name[i] = record["res_name"]
            structure.res_id[i] = record["res_id"]
            structure.chain_id[i] = record["chain_id"]

        # Initialize optional annotations that have data
        has_alt_loc = any(record["alt_loc"] for record in records)
        has_insertion = any(record["insertion_code"] for record in records)
        has_occupancy = any(record["occupancy"] != 1.0 for record in records)
        has_b_factor = any(record["b_factor"] != 0.0 for record in records)
        has_charge = any(record["charge"] for record in records)
        has_hetero = any(record["is_hetero"] for record in records)

        if has_alt_loc:
            structure.alt_loc = np.array(
                [record["alt_loc"] for record in records], dtype="U1"
            )

        if has_insertion:
            structure.insertion_code = np.array(
                [record["insertion_code"] for record in records], dtype="U1"
            )

        # Always initialize occupancy and b_factor since they're commonly used
        structure.occupancy = np.array(
            [record["occupancy"] for record in records], dtype=np.float32
        )
        structure.b_factor = np.array(
            [record["b_factor"] for record in records], dtype=np.float32
        )

        if has_charge:
            structure.charge = np.array(
                [record["charge"] for record in records], dtype=np.float32
            )

        # Always initialize serial numbers for CONECT parsing
        structure.serial = np.array(
            [record["serial"] for record in records], dtype=np.int32
        )

        return structure

    def _parse_conect_records(self, conect_data: Any, structure: Structure) -> BondList:
        """
        Parse CONECT records and create BondList.

        Args:
            conect_data: CONECT records from pdbreader
            structure: Structure object for serial number mapping

        Returns:
            BondList with bonds from CONECT records
        """
        bonds = BondList()

        # Create mapping from serial numbers to atom indices
        serial_to_index = {}
        if structure.serial is not None:
            for i, serial in enumerate(structure.serial):
                serial_to_index[serial] = i

        # Process CONECT records
        for _, conect_row in conect_data.iterrows():
            parent_serial = _safe_convert_int(conect_row.get("parent", 0))

            if parent_serial not in serial_to_index:
                continue

            parent_index = serial_to_index[parent_serial]

            # Get bonded atoms list
            bonded_serials = conect_row.get("bonds", [])
            if not isinstance(bonded_serials, list):
                bonded_serials = [bonded_serials]

            for bonded_serial in bonded_serials:
                bonded_serial = _safe_convert_int(bonded_serial, 0)
                if bonded_serial > 0 and bonded_serial in serial_to_index:
                    bonded_index = serial_to_index[bonded_serial]

                    # Add bond with CONECT detection method
                    bonds.add_bond(
                        parent_index,
                        bonded_index,
                        bond_order=1.0,
                        bond_type="covalent",
                        detection_method="conect",
                    )

        return bonds
