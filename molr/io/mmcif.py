"""
mmCIF file parser for the space module using the mmcif package.

This module provides a clean mmCIF parser specifically designed for the NumPy-based
Structure class, converting mmcif output directly to NumPy arrays for optimal
performance.
"""

import math
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
    from mmcif.api.PdbxContainers import DataContainer
    from mmcif.io.PdbxReader import PdbxReader
except ImportError:
    raise ImportError(
        "mmcif package is required for mmCIF parsing. Install with: pip install mmcif"
    )


def _safe_convert_int(value: Any, default: int = 0) -> int:
    """Safely convert value to integer, handling NaN and None."""
    if value is None or value == "?" or value == ".":
        return default
    try:
        if isinstance(value, float) and math.isnan(value):
            return default
        return int(value)
    except (ValueError, TypeError):
        return default


def _safe_convert_float(value: Any, default: float = 0.0) -> float:
    """Safely convert value to float, handling NaN and None."""
    if value is None or value == "?" or value == ".":
        return default
    try:
        float_val = float(value)
        if math.isnan(float_val):
            return default
        return float_val
    except (ValueError, TypeError):
        return default


def _safe_convert_str(value: Any, default: str = "") -> str:
    """Safely convert value to string, handling None and mmCIF null values."""
    if value is None or value == "?" or value == ".":
        return default
    return str(value).strip()


class mmCIFParser:
    """
    mmCIF file parser for the space module.

    Designed specifically for the NumPy-based Structure class, this parser
    converts mmcif output directly to NumPy arrays for optimal performance.

    Features:
    - Direct conversion to NumPy arrays
    - Multi-model support for trajectories
    - Efficient memory usage
    - Full mmCIF annotation support
    - Chemical bond information from mmCIF data
    """

    def __init__(self) -> None:
        """Initialize the mmCIF parser."""
        self._current_structure: Optional[Structure] = None
        self._current_bonds: Optional[BondList] = None

    def parse_file(self, filename: str) -> Union[Structure, StructureEnsemble]:
        """
        Parse an mmCIF file and return a Structure or StructureEnsemble.

        Args:
            filename: Path to the mmCIF file

        Returns:
            Structure object for single model, StructureEnsemble for multi-model

        Raises:
            IOError: If file cannot be read
            ValueError: If mmCIF format is invalid
        """
        try:
            with open(filename, "r") as file_handle:
                reader = PdbxReader(file_handle)
                data_containers: List[Any] = []
                reader.read(data_containers)

                if not data_containers:
                    raise ValueError("No data found in mmCIF file")

                return self._convert_mmcif_data(data_containers[0])

        except Exception as e:
            raise IOError(f"Error parsing mmCIF file '{filename}': {e}")

    def parse_string(self, mmcif_content: str) -> Union[Structure, StructureEnsemble]:
        """
        Parse mmCIF content from a string.

        Args:
            mmcif_content: mmCIF file content as string

        Returns:
            Structure object with all atoms and annotations
        """
        import os
        import tempfile

        try:
            # Write to temporary file since mmcif expects file handle
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".cif", delete=False
            ) as f:
                f.write(mmcif_content)
                temp_filename = f.name

            try:
                return self.parse_file(temp_filename)
            finally:
                os.unlink(temp_filename)

        except Exception as e:
            raise ValueError(f"Error parsing mmCIF content: {e}")

    def _convert_mmcif_data(
        self, data_container: DataContainer
    ) -> Union[Structure, StructureEnsemble]:
        """
        Convert mmCIF data container to Structure or StructureEnsemble.

        Args:
            data_container: DataContainer from mmcif package

        Returns:
            Structure for single model, StructureEnsemble for multi-model
        """
        # Get atom site data
        atom_site_obj = data_container.getObj("atom_site")
        if not atom_site_obj:
            raise ValueError("No atom_site data found in mmCIF file")

        # Check if we have multiple models
        models = self._group_by_models(atom_site_obj)

        if len(models) == 1:
            # Single model - return Structure
            return self._create_single_structure(models[0], data_container)
        else:
            # Multiple models - return StructureEnsemble
            return self._create_structure_ensemble(models, data_container)

    def _group_by_models(self, atom_site_obj: Any) -> List[List[Dict[str, Any]]]:
        """
        Group atom records by model number.

        Args:
            atom_site_obj: atom_site object from mmCIF data

        Returns:
            List of models, each containing atom records for that model
        """
        models: Dict[int, List] = {}  # model_id -> atom_records

        # Get the number of rows
        row_count = atom_site_obj.getRowCount()

        for row_idx in range(row_count):
            # Get model information - use safe access
            available_attrs = atom_site_obj.getAttributeList()
            if "pdbx_PDB_model_num" in available_attrs:
                model_num = _safe_convert_int(
                    atom_site_obj.getValue("pdbx_PDB_model_num", row_idx), 1
                )
            else:
                model_num = 1  # Default to single model

            if model_num not in models:
                models[model_num] = []

            # Convert atom record
            atom_record = self._convert_atom_record(atom_site_obj, row_idx)
            models[model_num].append(atom_record)

        if not models:
            raise ValueError("No atom records found in mmCIF file")

        # Return models in order
        model_ids = sorted(models.keys())
        return [models[model_id] for model_id in model_ids]

    def _create_single_structure(
        self, atom_records: List[Dict[str, Any]], data_container: DataContainer
    ) -> Structure:
        """
        Create single Structure from atom records.

        Args:
            atom_records: List of atom record dictionaries
            data_container: Original mmCIF data container for bond information

        Returns:
            Structure object
        """
        structure = self._create_structure_from_records(atom_records)

        # Parse bond information if available
        bonds = self._parse_bonds(data_container, structure)
        if bonds is not None:
            structure.file_bonds = bonds  # Store as file bonds

        return structure

    def _create_structure_ensemble(
        self, models: List[List[Dict[str, Any]]], data_container: DataContainer
    ) -> StructureEnsemble:
        """
        Create StructureEnsemble from multiple models.

        Args:
            models: List of models, each containing atom records
            data_container: Original mmCIF data container for bond information

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

        # Parse bond information if available (applied to template)
        bonds = self._parse_bonds(data_container, template)
        if bonds is not None:
            template.file_bonds = bonds  # Store as file bonds

        return ensemble

    def _convert_atom_record(self, atom_site_obj: Any, row_idx: int) -> Dict[str, Any]:
        """
        Convert a single atom record from mmCIF atom_site data.

        Args:
            atom_site_obj: atom_site object from mmCIF
            row_idx: Row index in the atom_site table

        Returns:
            Dictionary with atom data
        """
        # Get available attributes to handle missing fields gracefully
        available_attrs = atom_site_obj.getAttributeList()

        def safe_get_value(attr_name: str, default: Any = None) -> Any:
            """Safely get value from atom_site object."""
            if attr_name in available_attrs:
                return atom_site_obj.getValue(attr_name, row_idx)
            return default

        # Extract basic information using mmCIF column names
        atom_id = _safe_convert_int(safe_get_value("id"), 0)
        atom_name = _safe_convert_str(safe_get_value("label_atom_id")).upper()
        alt_loc = _safe_convert_str(safe_get_value("label_alt_id"))
        res_name = _safe_convert_str(safe_get_value("label_comp_id")).upper()
        chain_id = _safe_convert_str(safe_get_value("label_asym_id"))
        res_id = _safe_convert_int(safe_get_value("label_seq_id"), 0)
        insertion_code = _safe_convert_str(safe_get_value("pdbx_PDB_ins_code"))

        # Extract coordinates
        x = _safe_convert_float(safe_get_value("Cartn_x"), 0.0)
        y = _safe_convert_float(safe_get_value("Cartn_y"), 0.0)
        z = _safe_convert_float(safe_get_value("Cartn_z"), 0.0)

        # Extract optional fields
        occupancy = _safe_convert_float(safe_get_value("occupancy"), 1.0)
        b_factor = _safe_convert_float(safe_get_value("B_iso_or_equiv"), 0.0)
        element = _safe_convert_str(safe_get_value("type_symbol"))
        formal_charge = _safe_convert_str(safe_get_value("pdbx_formal_charge"))

        # Get group type (ATOM or HETATM)
        group_pdb = _safe_convert_str(safe_get_value("group_PDB"), "ATOM")
        is_hetero = group_pdb == "HETATM"

        # Determine element if not provided
        if not element:
            element = pdb_atom_to_element(atom_name)

        # Handle charge - convert to float if possible
        try:
            charge_value = float(formal_charge) if formal_charge else 0.0
        except (ValueError, TypeError):
            charge_value = 0.0

        return {
            "serial": atom_id,
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
            "charge": charge_value,
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
        has_charge = any(record["charge"] != 0.0 for record in records)
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

        # Always initialize serial numbers for bond parsing
        structure.serial = np.array(
            [record["serial"] for record in records], dtype=np.int32
        )

        return structure

    def _parse_bonds(
        self, data_container: DataContainer, structure: Structure
    ) -> Optional[BondList]:
        """
        Parse bond information from mmCIF data.

        Args:
            data_container: mmCIF data container
            structure: Structure object for atom index mapping

        Returns:
            BondList with bonds from mmCIF data, or None if no bond data
        """
        # Try different bond tables that might be present
        bond_tables = ["chem_comp_bond", "struct_conn", "geom_bond"]

        for table_name in bond_tables:
            bond_obj = data_container.getObj(table_name)
            if bond_obj and bond_obj.getRowCount() > 0:
                return self._parse_bond_table(bond_obj, table_name, structure)

        return None

    def _parse_bond_table(
        self, bond_obj: Any, table_name: str, structure: Structure
    ) -> BondList:
        """
        Parse specific bond table from mmCIF data.

        Args:
            bond_obj: Bond table object from mmCIF
            table_name: Name of the bond table
            structure: Structure object for atom mapping

        Returns:
            BondList with parsed bonds
        """
        bonds = BondList()

        # Create mapping from atom identifiers to indices
        atom_to_index = self._create_atom_mapping(structure)

        row_count = bond_obj.getRowCount()

        for row_idx in range(row_count):
            if table_name == "chem_comp_bond":
                # Chemical component bonds
                atom1_id = _safe_convert_str(bond_obj.getValue("atom_id_1", row_idx))
                atom2_id = _safe_convert_str(bond_obj.getValue("atom_id_2", row_idx))
                bond_order_str = _safe_convert_str(
                    bond_obj.getValue("value_order", row_idx), "SING"
                )

            elif table_name == "struct_conn":
                # Structural connections
                atom1_id = _safe_convert_str(
                    bond_obj.getValue("ptnr1_label_atom_id", row_idx)
                )
                atom2_id = _safe_convert_str(
                    bond_obj.getValue("ptnr2_label_atom_id", row_idx)
                )
                bond_order_str = "SING"  # struct_conn doesn't usually specify order

            elif table_name == "geom_bond":
                # Geometric bonds
                atom1_id = _safe_convert_str(bond_obj.getValue("atom_id_1", row_idx))
                atom2_id = _safe_convert_str(bond_obj.getValue("atom_id_2", row_idx))
                bond_order_str = "SING"  # geom_bond doesn't specify order

            else:
                continue

            # Map bond order
            bond_order = self._map_bond_order(bond_order_str)

            # Find atom indices
            atom1_idx = atom_to_index.get(atom1_id)
            atom2_idx = atom_to_index.get(atom2_id)

            if atom1_idx is not None and atom2_idx is not None:
                bonds.add_bond(
                    atom1_idx,
                    atom2_idx,
                    bond_order=bond_order,
                    bond_type="covalent",
                    detection_method="mmcif",
                )

        return bonds

    def _create_atom_mapping(self, structure: Structure) -> Dict[str, int]:
        """
        Create mapping from atom identifiers to structure indices.

        Args:
            structure: Structure object

        Returns:
            Dictionary mapping atom names to indices
        """
        mapping = {}
        for i in range(structure.n_atoms):
            # Use atom name as identifier (could be extended for more complex mapping)
            atom_id = structure.atom_name[i]
            mapping[atom_id] = i
        return mapping

    def _map_bond_order(self, bond_order_str: str) -> float:
        """
        Map mmCIF bond order string to numeric value.

        Args:
            bond_order_str: Bond order string from mmCIF

        Returns:
            Numeric bond order
        """
        order_mapping = {
            "SING": 1.0,
            "DOUB": 2.0,
            "TRIP": 3.0,
            "QUAD": 4.0,
            "AROM": 1.5,
            "POLY": 1.0,
            "DELO": 0.5,
            "PI": 1.0,
        }

        return order_mapping.get(bond_order_str.upper(), 1.0)
