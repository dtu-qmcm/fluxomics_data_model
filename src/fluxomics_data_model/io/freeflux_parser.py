"""
Freeflux tabular format parser for tsv/csv/xlsx files.

Freeflux is a Python package for 13C metabolic flux analysis that uses
tabular input files. This parser supports the following file types:
    - reactions: Network definition with atom mappings
    - fluxes: Flux values (simulation/reference)
    - concentrations: Metabolite pool sizes
    - measured_MDVs: Steady-state mass distribution vector measurements
    - measured_fluxes: Measured flux values with uncertainties
    - measured_inst_MDVs: Time-course MDV measurements

Supports .tsv, .csv, and .xlsx file formats.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from collections import defaultdict
import pandas as pd

from ..core.core import FluxomicsDataModel, Metadata, Model, Experiments
from ..core.common import DictList, TextualOrMath, ErrorModel
from ..model.metabolite import Metabolite
from ..model.reaction import Reaction
from ..model.atom_mapping import AtomMapping
from ..experiment.tracer import Tracers
from ..experiment.measurement import (
    Measurement,
    MeasurementModel,
    MeasurementData,
    LabelingMeasurement,
    FluxMeasurement,
    Group,
    NetFlux,
    Datum,
)
from ..output.simulation import Simulation, Variables, FluxValue, MetaboliteSizeValue


class FreefluxParser:
    """Parser for Freeflux tabular format files (tsv/csv/xlsx).

    Freeflux uses a set of tabular files to define metabolic networks and
    experimental data for 13C metabolic flux analysis.

    Example usage:
        parser = FreefluxParser()
        model = parser.parse("data/freeflux_test/toy")
        # or
        model = parser.parse("data/freeflux_test/ecoli/synthetic_data")
    """

    # File patterns for each data type
    FILE_PATTERNS = {
        "reactions": "reactions",
        "fluxes": "fluxes",
        "concentrations": "concentrations",
        "measured_MDVs": "measured_MDVs",
        "measured_fluxes": "measured_fluxes",
        "measured_inst_MDVs": "measured_inst_MDVs",
    }

    # Supported file extensions
    EXTENSIONS = [".tsv", ".csv", ".xlsx"]

    def __init__(self):
        self._metabolites: Dict[str, Metabolite] = {}
        self._reactions: Dict[str, Reaction] = {}
        self._atom_mappings: Dict[str, AtomMapping] = {}

    def parse(self, base_path: str) -> FluxomicsDataModel:
        """
        Parse Freeflux files from a directory.

        Args:
            base_path: Path to directory containing Freeflux files.

        Returns:
            FluxomicsDataModel containing the parsed data
        """
        base_path = Path(base_path)

        # If given a file, use its parent directory
        if base_path.is_file():
            base_path = base_path.parent

        # Find and parse reactions file (required)
        reactions_path = self._find_file(base_path, "reactions")
        if not reactions_path:
            raise FileNotFoundError(
                f"Reactions file not found in: {base_path}"
            )

        self._parse_reactions(reactions_path)

        # Build Model from parsed network
        model = Model(
            metabolites=DictList(list(self._metabolites.values())),
            reactions=DictList(list(self._reactions.values())),
            atom_mappings=self._atom_mappings,
        )

        # Parse optional files
        fluxes = self._parse_fluxes(base_path)
        concentrations = self._parse_concentrations(base_path)
        measurement = self._parse_measurements(base_path)

        # Create experiment if we have any experimental data
        experiments = []
        simulation = None

        # Build simulation from fluxes and concentrations
        if fluxes or concentrations:
            flux_values = fluxes or []
            metabolitesize_values = concentrations or []
            simulation = Simulation(
                variables=Variables(
                    flux_values=flux_values,
                    metabolitesize_values=metabolitesize_values,
                ),
            )

        if measurement or simulation:
            # Determine if we have inst_MDVs (non-stationary)
            has_inst_mdvs = self._find_file(base_path, "measured_inst_MDVs") is not None

            experiment = Experiments(
                name=base_path.name,
                stationary=not has_inst_mdvs,
                tracers=[],
                measurement=measurement,
                simulation=simulation,
            )
            experiments.append(experiment)

        # Create metadata
        metadata = Metadata(
            name=base_path.name,
            comment=f"Imported from Freeflux format: {base_path.name}",
        )

        return FluxomicsDataModel(
            model=model,
            info=metadata,
            experiments=experiments,
        )

    def _find_file(self, base_path: Path, file_type: str) -> Optional[Path]:
        """Find a file matching the pattern with any supported extension."""
        pattern = self.FILE_PATTERNS[file_type]
        for ext in self.EXTENSIONS:
            file_path = base_path / f"{pattern}{ext}"
            if file_path.exists():
                return file_path
        return None

    def _read_tabular(self, file_path: Path) -> pd.DataFrame:
        """Read tabular file (tsv/csv/xlsx) into DataFrame."""
        suffix = file_path.suffix.lower()

        if suffix == ".xlsx":
            df = pd.read_excel(file_path, comment="#")
        elif suffix == ".csv":
            df = pd.read_csv(file_path, comment="#")
        else:  # .tsv
            df = pd.read_csv(file_path, sep="\t", comment="#")

        # Clean column names - remove leading # and strip whitespace
        df.columns = [
            col.lstrip("#").strip() if isinstance(col, str) else col
            for col in df.columns
        ]

        return df

    def _parse_reactions(self, filepath: Path) -> None:
        """
        Parse reactions file containing reaction network with atom mappings.

        Format:
        - #reaction_ID: Reaction identifier (section headers start with #)
        - substrate_IDs(atom) or reactant_IDs(atom): Substrates with atom mapping
        - product_IDs(atom): Products with atom mapping
        - reversibility: 0 (irreversible) or 1 (reversible)
        """
        self._metabolites = {}
        self._reactions = {}
        self._atom_mappings = {}

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "reaction" in col_lower and "id" in col_lower:
                col_map["reaction_id"] = col
            elif "substrate" in col_lower or "reactant" in col_lower:
                col_map["substrates"] = col
            elif "product" in col_lower:
                col_map["products"] = col
            elif "reversib" in col_lower:
                col_map["reversibility"] = col

        for _, row in df.iterrows():
            rxn_id = str(row.get(col_map.get("reaction_id", ""), "")).strip()

            # Skip empty rows or section headers (start with #)
            if not rxn_id or rxn_id.startswith("#") or pd.isna(row.get(col_map.get("reaction_id", ""))):
                continue

            substrates_str = str(row.get(col_map.get("substrates", ""), ""))
            products_str = str(row.get(col_map.get("products", ""), ""))
            rev_value = row.get(col_map.get("reversibility", ""), 0)

            # Skip if no substrates or products
            if pd.isna(substrates_str) or substrates_str == "nan" or not substrates_str.strip():
                continue
            if pd.isna(products_str) or products_str == "nan" or not products_str.strip():
                continue

            # Parse reversibility
            reversible = False
            if not pd.isna(rev_value):
                if isinstance(rev_value, (int, float)):
                    reversible = int(rev_value) == 1
                elif isinstance(rev_value, str):
                    reversible = rev_value.strip().lower() in ("1", "true", "yes")

            # Parse substrates and products with atom mappings
            reactants, reactant_atoms = self._parse_compounds(substrates_str)
            products, product_atoms = self._parse_compounds(products_str)

            # Create metabolites if not already defined
            for cpd_id in reactants + products:
                if cpd_id not in self._metabolites:
                    self._metabolites[cpd_id] = Metabolite(id=cpd_id)

            # Create reaction
            reaction = Reaction(
                id=rxn_id,
                reactants=reactants,
                products=products,
                reversibility=reversible,
            )
            self._reactions[rxn_id] = reaction

            # Create atom mapping if atoms are specified
            if reactant_atoms and product_atoms:
                atom_mapping = self._create_atom_mapping(
                    rxn_id, reactants, products, reactant_atoms, product_atoms
                )
                if atom_mapping:
                    self._atom_mappings[rxn_id] = atom_mapping

    def _parse_compounds(
        self, compounds_str: str
    ) -> Tuple[List[str], List[Tuple[str, str]]]:
        """
        Parse compounds string with optional atom mappings.

        Args:
            compounds_str: String like "G6P(abcdef)+AcCoA(gh)" or "G6P(abcdef) + AcCoA(gh)"

        Returns:
            Tuple of (compound_ids, [(compound_id, atoms), ...])
        """
        compound_ids = []
        compound_atoms = []

        # Split by + and parse each compound
        parts = compounds_str.split("+")

        for part in parts:
            part = part.strip()
            if not part:
                continue

            # Match compound with optional atom mapping: "Compound(atoms)" or "123Compound(atoms)"
            # Also handles stoichiometry like "2NADPH" or "0.5O2"
            match = re.match(r"^([\d.]*)?(\S+?)(?:\(([^)]+)\))?$", part)
            if match:
                stoich = match.group(1) or ""
                cpd_id = match.group(2)
                atoms = match.group(3) or ""

                # Handle symmetric compounds with comma-separated atoms like "abcd,dcba"
                # Just use the first variant for the compound id
                atoms_clean = atoms.split(",")[0] if atoms else ""

                compound_ids.append(cpd_id)
                if atoms_clean:
                    compound_atoms.append((cpd_id, atoms))

        return compound_ids, compound_atoms

    def _create_atom_mapping(
        self,
        rxn_id: str,
        reactants: List[str],
        products: List[str],
        reactant_atoms: List[Tuple[str, str]],
        product_atoms: List[Tuple[str, str]],
    ) -> Optional[AtomMapping]:
        """
        Create AtomMapping from letter notation.

        Handles symmetric compounds with comma-separated atom variants.
        """
        # Check for symmetric compounds (comma-separated atom strings)
        has_variants = any("," in atoms for _, atoms in reactant_atoms + product_atoms)

        if has_variants:
            # Generate all variant combinations
            return self._create_variant_atom_mapping(
                rxn_id, reactants, products, reactant_atoms, product_atoms
            )
        else:
            # Simple case: single atom mapping
            atom_map = AtomMapping.parse_letter_notation(
                reactant_atoms, product_atoms
            )
            return AtomMapping(
                reaction_id=rxn_id,
                reactants=reactants,
                products=products,
                maps={rxn_id: atom_map},
                weights={rxn_id: 1.0},
            )

    def _create_variant_atom_mapping(
        self,
        rxn_id: str,
        reactants: List[str],
        products: List[str],
        reactant_atoms: List[Tuple[str, str]],
        product_atoms: List[Tuple[str, str]],
    ) -> AtomMapping:
        """Create atom mapping with variants for symmetric compounds."""
        from itertools import product as cartesian_product

        # Extract variants for each compound
        reactant_variants = []
        for cpd, atoms in reactant_atoms:
            variants = atoms.split(",")
            reactant_variants.append([(cpd, v.strip()) for v in variants])

        product_variants = []
        for cpd, atoms in product_atoms:
            variants = atoms.split(",")
            product_variants.append([(cpd, v.strip()) for v in variants])

        # Generate all combinations
        all_reactant_combos = list(cartesian_product(*reactant_variants)) if reactant_variants else [()]
        all_product_combos = list(cartesian_product(*product_variants)) if product_variants else [()]

        maps_dict = {}
        variant_idx = 1

        for r_combo in all_reactant_combos:
            for p_combo in all_product_combos:
                r_items = list(r_combo) if r_combo else reactant_atoms
                p_items = list(p_combo) if p_combo else product_atoms

                try:
                    atom_map = AtomMapping.parse_letter_notation(r_items, p_items)
                    map_id = f"{rxn_id}___{variant_idx}" if len(all_reactant_combos) * len(all_product_combos) > 1 else rxn_id
                    maps_dict[map_id] = atom_map
                    variant_idx += 1
                except Exception:
                    # Skip invalid combinations
                    continue

        if not maps_dict:
            return None

        # Calculate uniform weights
        uniform_weight = 1.0 / len(maps_dict)
        weights_dict = {k: uniform_weight for k in maps_dict.keys()}

        # Store variant IDs in reaction if multiple variants
        if len(maps_dict) > 1:
            if rxn_id in self._reactions:
                # Update reaction with variant IDs
                old_rxn = self._reactions[rxn_id]
                self._reactions[rxn_id] = Reaction(
                    id=rxn_id,
                    reactants=old_rxn.reactants,
                    products=old_rxn.products,
                    reversibility=old_rxn.reversibility,
                    atom_mapping_ids=list(maps_dict.keys()),
                )

        return AtomMapping(
            reaction_id=rxn_id,
            reactants=reactants,
            products=products,
            maps=maps_dict,
            weights=weights_dict,
        )

    def _parse_fluxes(self, base_path: Path) -> Optional[List[FluxValue]]:
        """
        Parse fluxes file containing flux values.

        Format:
        - #flux_ID: Flux identifier (e.g., "v1" or "v1_f", "v1_b" for reversible)
        - value: Flux value
        """
        filepath = self._find_file(base_path, "fluxes")
        if not filepath:
            return None

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "flux" in col_lower and "id" in col_lower:
                col_map["flux_id"] = col
            elif col_lower == "value":
                col_map["value"] = col

        flux_values = []
        for _, row in df.iterrows():
            flux_id = str(row.get(col_map.get("flux_id", ""), "")).strip()
            value = row.get(col_map.get("value", ""), None)

            if not flux_id or pd.isna(value):
                continue

            # Determine flux type from ID suffix (_f for forward, _b for backward)
            flux_type = "net"
            base_id = flux_id
            if flux_id.endswith("_f"):
                flux_type = "net"  # forward flux
                base_id = flux_id[:-2]
            elif flux_id.endswith("_b"):
                flux_type = "xch"  # backward as exchange
                base_id = flux_id[:-2]

            flux_value = FluxValue(
                flux=flux_id,
                type=flux_type,
                lo=float(value),
                hi=float(value),
                value=float(value),
            )
            flux_values.append(flux_value)

        return flux_values if flux_values else None

    def _parse_concentrations(self, base_path: Path) -> Optional[List[MetaboliteSizeValue]]:
        """
        Parse concentrations file containing metabolite pool sizes.

        Format:
        - #metabolite_ID or #metab_ID: Metabolite identifier
        - value or value (umol/gCDW): Concentration value
        """
        filepath = self._find_file(base_path, "concentrations")
        if not filepath:
            return None

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "metabolite" in col_lower or "metab" in col_lower:
                col_map["metab_id"] = col
            elif "value" in col_lower:
                col_map["value"] = col

        metabolitesize_values = []
        for _, row in df.iterrows():
            metab_id = str(row.get(col_map.get("metab_id", ""), "")).strip()
            value = row.get(col_map.get("value", ""), None)

            if not metab_id or pd.isna(value):
                continue

            metabolitesize_value = MetaboliteSizeValue(
                metabolite=metab_id,
                lo=float(value),
                hi=float(value),
                value=float(value),
            )
            metabolitesize_values.append(metabolitesize_value)

        return metabolitesize_values if metabolitesize_values else None

    def _parse_measurements(self, base_path: Path) -> Optional[Measurement]:
        """Parse measurement files (measured_MDVs, measured_fluxes, measured_inst_MDVs)."""
        # Parse steady-state MDV measurements
        labeling_measurement, labeling_data = self._parse_measured_mdvs(base_path)

        # Parse time-course MDV measurements
        inst_labeling, inst_data = self._parse_measured_inst_mdvs(base_path)

        # Parse flux measurements
        flux_measurement, flux_data = self._parse_measured_fluxes(base_path)

        # Merge labeling measurements
        if inst_labeling and labeling_measurement:
            # Combine groups
            all_groups = list(labeling_measurement.groups) + list(inst_labeling.groups)
            labeling_measurement = LabelingMeasurement(groups=all_groups)
        elif inst_labeling:
            labeling_measurement = inst_labeling

        # Combine data
        all_data = []
        if labeling_data:
            all_data.extend(labeling_data)
        if inst_data:
            all_data.extend(inst_data)
        if flux_data:
            all_data.extend(flux_data)

        if not labeling_measurement and not flux_measurement:
            return None

        return Measurement(
            model=MeasurementModel(
                labeling_measurement=labeling_measurement,
                flux_measurement=flux_measurement,
            ),
            data=MeasurementData(data=all_data),
        )

    def _parse_measured_mdvs(
        self, base_path: Path
    ) -> Tuple[Optional[LabelingMeasurement], Optional[List[Datum]]]:
        """
        Parse measured_MDVs file containing steady-state MDV measurements.

        Format:
        - #fragment_ID: Fragment identifier like "Glu_12345" (metabolite_positions)
        - mean: Comma-separated MDV values
        - sd: Comma-separated standard deviations
        """
        filepath = self._find_file(base_path, "measured_MDVs")
        if not filepath:
            return None, None

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "fragment" in col_lower and "id" in col_lower:
                col_map["fragment_id"] = col
            elif col_lower == "mean":
                col_map["mean"] = col
            elif col_lower == "sd":
                col_map["sd"] = col

        groups: Dict[str, Group] = {}
        data: List[Datum] = []

        for _, row in df.iterrows():
            fragment_id = str(row.get(col_map.get("fragment_id", ""), "")).strip()
            mean_str = str(row.get(col_map.get("mean", ""), ""))
            sd_str = str(row.get(col_map.get("sd", ""), ""))

            if not fragment_id or pd.isna(mean_str) or mean_str == "nan":
                continue

            # Parse fragment_id to extract metabolite and positions
            # Format: "Metabolite_positions" e.g., "Glu_12345"
            parts = fragment_id.rsplit("_", 1)
            if len(parts) == 2:
                metabolite = parts[0]
                positions = parts[1]
            else:
                metabolite = fragment_id
                positions = ""

            # Create group if not exists
            if fragment_id not in groups:
                # Build expression for MS measurement
                if positions:
                    pos_list = "+".join(positions)
                    expression = f"{metabolite}[{pos_list}]"
                else:
                    expression = metabolite

                groups[fragment_id] = Group(
                    id=fragment_id,
                    scale="auto",
                    expression=TextualOrMath(textual=expression),
                )

            # Parse mean and sd values
            mean_values = self._parse_comma_values(mean_str)
            sd_values = self._parse_comma_values(sd_str)

            # Create datum for each MDV component (M0, M1, M2, ...)
            for i, (mean, sd) in enumerate(zip(mean_values, sd_values)):
                datum_id = f"{fragment_id}:M{i}"
                datum = Datum(
                    id=datum_id,
                    value=mean,
                    stddev=sd if sd else 0.01,
                    pos=i,
                )
                data.append(datum)

        if not groups:
            return None, None

        return (
            LabelingMeasurement(groups=list(groups.values())),
            data
        )

    def _parse_measured_inst_mdvs(
        self, base_path: Path
    ) -> Tuple[Optional[LabelingMeasurement], Optional[List[Datum]]]:
        """
        Parse measured_inst_MDVs file containing time-course MDV measurements.

        Format:
        - #fragment_ID: Fragment identifier
        - time or time (s): Time point
        - mean: Comma-separated MDV values
        - sd: Comma-separated standard deviations
        """
        filepath = self._find_file(base_path, "measured_inst_MDVs")
        if not filepath:
            return None, None

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "fragment" in col_lower and "id" in col_lower:
                col_map["fragment_id"] = col
            elif "time" in col_lower:
                col_map["time"] = col
            elif col_lower == "mean":
                col_map["mean"] = col
            elif col_lower == "sd":
                col_map["sd"] = col

        groups: Dict[str, Group] = {}
        data: List[Datum] = []
        time_points: Dict[str, List[float]] = defaultdict(list)

        for _, row in df.iterrows():
            fragment_id = str(row.get(col_map.get("fragment_id", ""), "")).strip()
            time_val = row.get(col_map.get("time", ""), None)
            mean_str = str(row.get(col_map.get("mean", ""), ""))
            sd_str = str(row.get(col_map.get("sd", ""), ""))

            if not fragment_id or pd.isna(mean_str) or mean_str == "nan":
                continue

            # Parse time
            time_float = float(time_val) if not pd.isna(time_val) else None

            # Track time points for this fragment
            if time_float is not None:
                time_points[fragment_id].append(time_float)

            # Parse fragment_id
            parts = fragment_id.rsplit("_", 1)
            if len(parts) == 2:
                metabolite = parts[0]
                positions = parts[1]
            else:
                metabolite = fragment_id
                positions = ""

            # Create group if not exists (will update times later)
            if fragment_id not in groups:
                if positions:
                    pos_list = "+".join(positions)
                    expression = f"{metabolite}[{pos_list}]"
                else:
                    expression = metabolite

                groups[fragment_id] = Group(
                    id=fragment_id,
                    scale="auto",
                    expression=TextualOrMath(textual=expression),
                )

            # Parse mean and sd values
            mean_values = self._parse_comma_values(mean_str)
            sd_values = self._parse_comma_values(sd_str)

            # Create datum for each MDV component
            for i, (mean, sd) in enumerate(zip(mean_values, sd_values)):
                datum_id = f"{fragment_id}:M{i}"
                datum = Datum(
                    id=datum_id,
                    value=mean,
                    stddev=sd if sd else 0.01,
                    time=time_float,
                    pos=i,
                )
                data.append(datum)

        if not groups:
            return None, None

        # Update groups with time information
        updated_groups = []
        for group_id, group in groups.items():
            times = sorted(set(time_points.get(group_id, [])))
            times_str = ",".join(str(t) for t in times) if times else None
            updated_groups.append(Group(
                id=group.id,
                times=times_str,
                scale=group.scale,
                expression=group.expression,
            ))

        return (
            LabelingMeasurement(groups=updated_groups),
            data
        )

    def _parse_measured_fluxes(
        self, base_path: Path
    ) -> Tuple[Optional[FluxMeasurement], Optional[List[Datum]]]:
        """
        Parse measured_fluxes file containing flux measurements.

        Format:
        - #reaction_ID: Reaction identifier
        - mean: Mean flux value
        - sd: Standard deviation
        """
        filepath = self._find_file(base_path, "measured_fluxes")
        if not filepath:
            return None, None

        df = self._read_tabular(filepath)

        # Standardize column names
        col_map = {}
        for col in df.columns:
            col_lower = col.lower()
            if "reaction" in col_lower and "id" in col_lower:
                col_map["reaction_id"] = col
            elif col_lower == "mean":
                col_map["mean"] = col
            elif col_lower == "sd":
                col_map["sd"] = col

        net_fluxes: List[NetFlux] = []
        data: List[Datum] = []

        for _, row in df.iterrows():
            rxn_id = str(row.get(col_map.get("reaction_id", ""), "")).strip()
            mean = row.get(col_map.get("mean", ""), None)
            sd = row.get(col_map.get("sd", ""), 0.01)

            if not rxn_id or pd.isna(mean):
                continue

            # Create net flux measurement
            sd_value = float(sd) if not pd.isna(sd) else 0.01
            net_flux = NetFlux(
                id=rxn_id,
                expression=TextualOrMath(textual=rxn_id),
                errormodel=ErrorModel(
                    expression=TextualOrMath(textual=str(sd_value))
                ),
            )
            net_fluxes.append(net_flux)

            # Create datum
            datum = Datum(
                id=rxn_id,
                value=float(mean),
                stddev=sd_value,
            )
            data.append(datum)

        if not net_fluxes:
            return None, None

        return (
            FluxMeasurement(net_fluxes=net_fluxes),
            data
        )

    def _parse_comma_values(self, value_str: str) -> List[float]:
        """Parse comma-separated values into list of floats."""
        if not value_str or value_str == "nan":
            return []

        # Remove quotes if present
        value_str = value_str.strip('"\'')

        values = []
        for v in value_str.split(","):
            v = v.strip()
            try:
                values.append(float(v))
            except ValueError:
                values.append(0.0)
        return values


def parse_freeflux(base_path: str) -> FluxomicsDataModel:
    """
    Convenience function to parse Freeflux files.

    Args:
        base_path: Path to directory containing Freeflux files.

    Returns:
        FluxomicsDataModel containing the parsed data
    """
    parser = FreefluxParser()
    return parser.parse(base_path)