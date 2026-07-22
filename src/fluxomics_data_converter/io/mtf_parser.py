"""MTF (Metabolic Text Format) parser for influx_si format.

MTF is a multi-file format used by influx_si for 13C metabolic flux analysis.
Each model consists of multiple files with a shared basename:
    - .netw: Network definition (reactions with atom transitions)
    - .linp: Label input (tracer specifications)
    - .miso: MS isotopomer measurements
    - .mflux: Flux measurements
    - .cnstr: Constraints (NET and XCH)
    - .tvar: Variable types and starting values
    - .opt: Options/command arguments
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import pandas as pd

from ..core.core import (
    FluxomicsData,
    Metadata,
    MetabolicNetworkModel,
    LabelingExperiments,
)
from ..core.common import DictList, TextualOrMath, ErrorModel
from ..model.metabolite import Metabolite
from ..model.reaction import Reaction
from ..model.atom_mapping import AtomTransition
from ..model.constraint import (
    Constraints,
    NetConstraints,
    ExchangeConstraints,
    ConstraintFormula,
)
from ..experiment.tracer import Tracers, LabelComposition
from ..experiment.measurement import (
    Measurement,
    MeasurementModel,
    MeasurementData,
    LabelingMeasurement,
    FluxMeasurement,
    MetaboliteSizeMeasurement,
    MetaboliteSize,
    Group,
    NetFlux,
    Datum,
)
from ..output.simulation import (
    Simulation,
    Variables,
    FluxValue,
    MetaboliteSizeValue,
)


class MTFParser:
    """Parser for MTF (Metabolic Text Format) files used by influx_si.

    MTF is a multi-file format. Each model is represented as a set of plain
    text files sharing the same base name but with different extensions::

        model.netw   — reaction network and atom transitions  (required)
        model.linp   — label input (tracer fractions)
        model.miso   — MS isotopomer measurements
        model.mflux  — flux measurements
        model.mmet   — metabolite pool-size measurements
        model.cnstr  — stoichiometric constraints
        model.tvar   — flux/metabolite variable types and initial values
        model.opt    — solver options (parsed but not used)

    Example::

        parser = MTFParser()
        model = parser.parse("data/ecoli")          # base path
        model = parser.parse("data/ecoli.netw")     # or the .netw file itself

    Reaction arrow conventions
    --------------------------
    - ``->``  : irreversible (net flux unconstrained)
    - ``->>`` : irreversible, net flux ≥ 0 (auto-generates a constraint)
    - ``<->`` : reversible (bidirectional)
    - ``<->>`` : reversible, net flux ≥ 0 (auto-generates a constraint)
    """

    # File extensions in MTF format
    EXTENSIONS = {
        "netw": ".netw",
        "linp": ".linp",
        "miso": ".miso",
        "mflux": ".mflux",
        "mmet": ".mmet",
        "cnstr": ".cnstr",
        "tvar": ".tvar",
        "opt": ".opt",
    }

    def __init__(self):
        """Initialise parser state; call parse() to read MTF files."""
        self._metabolites: Dict[str, Metabolite] = {}
        self._reactions: Dict[str, Reaction] = {}
        self._atom_mappings: Dict[str, AtomTransition] = {}
        # Track reactions with <->> arrow (reversible but net flux >= 0)
        self._net_positive_reactions: List[str] = []

    def parse(self, base_path: str) -> FluxomicsData:
        """Parse MTF files and return a FluxomicsData.

        Args:
            base_path: Path to base filename (without extension) or
                any MTF file.
                       For example: "model/e_coli" or "model/e_coli.netw"

        Returns:
            FluxomicsData containing the parsed data
        """
        base_path = Path(base_path)

        # If given a file with extension, strip it to get base path
        if base_path.suffix in self.EXTENSIONS.values():
            base_path = base_path.with_suffix("")

        # Check that at least .netw file exists
        netw_path = base_path.with_suffix(".netw")
        if not netw_path.exists():
            raise FileNotFoundError(f"Network file not found: {netw_path}")

        # Parse network file (required)
        self._parse_network(netw_path)

        # Build Model from parsed network
        model = MetabolicNetworkModel(
            metabolites=DictList(list(self._metabolites.values())),
            reactions=DictList(list(self._reactions.values())),
            atom_mappings=self._atom_mappings,
        )

        # Parse optional files
        constraints = self._parse_constraints(base_path)
        tracers = self._parse_label_input(base_path)
        measurement = self._parse_measurements(base_path)
        simulation = self._parse_variables(base_path)

        # Create experiment if we have any experimental data
        experiments = []
        if tracers or measurement or simulation:
            experiment = LabelingExperiments(
                name=base_path.stem,
                stationary=True,
                tracers=tracers or [],
                constraints=constraints,
                measurement=measurement,
                simulation=simulation,
            )
            experiments.append(experiment)

        # Create metadata
        metadata = Metadata(
            name=base_path.stem,
            comment=f"Imported from MTF format: {base_path.name}",
        )

        return FluxomicsData(
            model=model,
            metadata=metadata,
            constraints=constraints,
            experiments=experiments,
        )

    def _parse_network(self, filepath: Path) -> None:
        r"""Parse .netw file containing reaction network with atom transitions.

        Format: reaction_id:\tsubstrate (ATOMS) + substrate (atoms) ->
        product (ATOMS)
        - -> for irreversible reactions
        - <-> for reversible reactions
        - Atom transitions in parentheses using letter notation
        """
        self._metabolites = {}
        self._reactions = {}
        self._atom_mappings = {}
        self._net_positive_reactions = []

        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith("#"):
                    continue

                # Parse reaction line
                reaction = self._parse_reaction_line(line)
                if reaction:
                    rxn_obj, atom_mapping = reaction
                    self._reactions[rxn_obj.id] = rxn_obj
                    if atom_mapping:
                        self._atom_mappings[rxn_obj.id] = atom_mapping

    def _parse_reaction_line(
        self, line: str
    ) -> Optional[Tuple[Reaction, Optional[AtomTransition]]]:
        """Parse a single reaction line from .netw file.

        Returns:
            Tuple of (Reaction, AtomTransition) or None if line is invalid

        Arrow types (from influx_si documentation):
            ->   : non-reversible (exchange flux = 0, net flux can be +/-)
            ->>  : non-reversible with non-negative net flux (net flux >= 0)
            <->  : reversible
            <->> : reversible but net flux must be non-negative
        """
        # Match reaction format: id:\tsubstrates -> products
        # Order matters: match longer patterns first (<->>, ->>)
        # before shorter (<->, ->)
        # The product side uses ``.*`` (not ``.+``) so drain/output reactions
        # written with an empty product side (e.g. ``R64: PDO <->``) still
        # parse; otherwise they are dropped and any flux variable that
        # references them fails validation on round-trip.
        match = re.match(r"^(\S+):\s*(.+?)\s*(<->>|<->|->>|->)\s*(.*)$", line)
        if not match:
            return None

        rxn_id = match.group(1)
        reactants_str = match.group(2)
        arrow = match.group(3)
        products_str = match.group(4)

        # <-> and <->> are reversible; -> and ->> are not
        reversible = arrow in ("<->", "<->>")

        # Track reactions with <->> arrow (reversible but net flux >= 0)
        # These need an automatic constraint: rxn_id >= 0
        if arrow == "<->>":
            self._net_positive_reactions.append(rxn_id)

        # Parse reactants and products with atom transitions
        reactants, reactant_atoms = self._parse_compounds(reactants_str)
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

        # Create atom transition if atoms are specified
        atom_mapping = None
        if reactant_atoms and product_atoms:
            atom_mapping = self._create_atom_mapping(
                rxn_id, reactants, products, reactant_atoms, product_atoms
            )

        return reaction, atom_mapping

    def _parse_compounds(
        self, compounds_str: str
    ) -> Tuple[List[str], List[Tuple[str, str]]]:
        """Parse compounds string with optional atom transitions.

        Args:
            compounds_str: String like "A (abc) + B (def)"

        Returns:
            Tuple of (compound_ids, [(compound_id, atoms), ...])
        """
        compound_ids = []
        compound_atoms = []

        # Split by + and parse each compound
        parts = compounds_str.split("+")

        for part in parts:
            part = part.strip()

            # Match compound with optional atom transition:
            # "Compound(atoms)" or "Compound (atoms)"
            match = re.match(r"^(\S+?)(?:\s*\(([^)]+)\))?$", part)
            if match:
                cpd_id = match.group(1)
                atoms = match.group(2) or ""

                compound_ids.append(cpd_id)
                if atoms:
                    compound_atoms.append((cpd_id, atoms))

        return compound_ids, compound_atoms

    def _create_atom_mapping(
        self,
        rxn_id: str,
        reactants: List[str],
        products: List[str],
        reactant_atoms: List[Tuple[str, str]],
        product_atoms: List[Tuple[str, str]],
    ) -> AtomTransition:
        """Create AtomTransition from letter notation.

        Args:
            rxn_id: Reaction ID
            reactants: Ordered list of reactant IDs
            products: Ordered list of product IDs
            reactant_atoms: List of (compound_id, atom_letters) tuples
            product_atoms: List of (compound_id, atom_letters) tuples

        Returns:
            AtomTransition object
        """
        # Use the static method from AtomTransition to parse letter notation
        atom_map = AtomTransition.parse_letter_notation(
            reactant_atoms, product_atoms
        )

        # Create single-map AtomTransition
        return AtomTransition(
            reaction_id=rxn_id,
            reactants=reactants,
            products=products,
            maps={rxn_id: atom_map},
            weights={rxn_id: 1.0},
        )

    def _parse_constraints(self, base_path: Path) -> Optional[Constraints]:
        r"""Parse .cnstr file containing constraints.

        Format (TSV): Id\tComment\tKind\tFormula\tOperator\tValue
        - Kind: NET or XCH
        - Operator: ==, >=, <=

        Also adds automatic constraints for <->> reactions (net flux >= 0).
        """
        net_formulas = []
        xch_formulas = []

        cnstr_path = base_path.with_suffix(".cnstr")
        if cnstr_path.exists():
            try:
                df = pd.read_csv(
                    cnstr_path,
                    sep="\t",
                    comment="#",
                    header=0,
                    names=[
                        "Id",
                        "Comment",
                        "Kind",
                        "Formula",
                        "Operator",
                        "Value",
                    ],
                )

                for _, row in df.iterrows():
                    kind = str(row.get("Kind", "")).strip().upper()
                    formula = str(row.get("Formula", "")).strip()
                    operator = str(row.get("Operator", "")).strip()
                    value = row.get("Value", "")

                    if not formula or pd.isna(formula):
                        continue

                    # Build expression string
                    op_map = {"==": "=", ">=": ">=", "<=": "<="}
                    op_str = op_map.get(operator, "=")

                    # Handle value
                    if pd.isna(value) or value == "":
                        expression = formula
                    else:
                        expression = f"{formula} {op_str} {value}"

                    constraint = ConstraintFormula(expression=expression)

                    if kind == "NET":
                        net_formulas.append(constraint)
                    elif kind == "XCH":
                        xch_formulas.append(constraint)

            except Exception as e:
                raise ValueError(
                    f"Error parsing constraints file {cnstr_path}: {e}"
                )

        # Add automatic constraints for <->> reactions (net flux >= 0)
        for rxn_id in self._net_positive_reactions:
            constraint = ConstraintFormula(
                name=f"net_positive_{rxn_id}",
                expression=f"{rxn_id} >= 0",
            )
            net_formulas.append(constraint)

        if not net_formulas and not xch_formulas:
            return None

        return Constraints(
            net=NetConstraints(formulas=net_formulas) if net_formulas else None,
            xch=ExchangeConstraints(formulas=xch_formulas)
            if xch_formulas
            else None,
        )

    def _parse_label_input(self, base_path: Path) -> Optional[List[Tracers]]:
        r"""Parse .linp file containing tracer specifications.

        Format (TSV): Id\tComment\tMetabolite\tIsotopomer\tValue
        - Isotopomer: Binary string (e.g., "111111", "100000")
        - Value: Fraction (0-1)
        """
        linp_path = base_path.with_suffix(".linp")
        if not linp_path.exists():
            return None

        tracers_dict: Dict[str, List[LabelComposition]] = defaultdict(list)

        try:
            # Force Isotopomer to be read as string to preserve
            # leading zeros (e.g., "000000")
            df = pd.read_csv(
                linp_path,
                sep="\t",
                comment="#",
                header=0,
                names=["Id", "Comment", "Metabolite", "Isotopomer", "Value"],
                dtype={"Isotopomer": str},
            )

            for _, row in df.iterrows():
                metabolite = str(row.get("Metabolite", "")).strip()
                isotopomer = (
                    str(row.get("Isotopomer", "")).strip()
                    if pd.notna(row.get("Isotopomer"))
                    else ""
                )
                value = row.get("Value", 1.0)

                if not metabolite or pd.isna(metabolite):
                    continue

                # Convert isotopomer to labeled_pattern (0/1 string)
                if isotopomer and not pd.isna(isotopomer):
                    labeled_pattern = isotopomer
                else:
                    continue

                fraction = float(value) if not pd.isna(value) else 1.0

                label = LabelComposition(
                    labeled_pattern=labeled_pattern,
                    fraction=fraction,
                )
                tracers_dict[metabolite].append(label)

        except Exception as e:
            raise ValueError(f"Error parsing label input file {linp_path}: {e}")

        if not tracers_dict:
            return None

        tracers = []
        for metabolite, labels in tracers_dict.items():
            tracer = Tracers(
                metabolite=metabolite,
                labels=labels,
            )
            tracers.append(tracer)

        return tracers

    def _parse_measurements(self, base_path: Path) -> Optional[Measurement]:
        """Parse measurement files (.miso, .mflux, .mmet)."""
        # Parse MS isotopomer measurements
        labeling_measurement, labeling_data = self._parse_miso(base_path)

        # Parse flux measurements
        flux_measurement, flux_data = self._parse_mflux(base_path)

        # Parse metabolite size measurements
        metabolitesize_measurement, mmet_data = self._parse_mmet(base_path)

        if (
            not labeling_measurement
            and not flux_measurement
            and not metabolitesize_measurement
        ):
            return None

        # Combine data
        all_data = (labeling_data or []) + (flux_data or []) + (mmet_data or [])

        return Measurement(
            model=MeasurementModel(
                labeling_measurement=labeling_measurement,
                flux_measurement=flux_measurement,
                metabolitesize_measurement=metabolitesize_measurement,
            ),
            data=MeasurementData(data=all_data),
        )

    def _parse_miso(
        self, base_path: Path
    ) -> Tuple[Optional[LabelingMeasurement], Optional[List[Datum]]]:
        r"""Parse .miso file containing MS isotopomer measurements.

        Format (TSV):
        Id\tComment\tSpecie\tFragment\tDataset\tIsospecies\tValue\tSD\tTime
        - Fragment: Comma-separated atom positions (e.g., "1,2,3,4")
        - Isospecies: Mass isotopomer label (M0, M1, M2, ...)

        Output FluxML format: <group id="MS-1"><textual>
        Suc[1-4]#M0,1,2,3,4</textual></group>
        """
        miso_path = base_path.with_suffix(".miso")
        if not miso_path.exists():
            return None, None

        # First pass: collect all data by group (dataset)
        group_info: Dict[
            str, Dict
        ] = {}  # dataset -> {specie, fragment, isospecies_set}
        data: List[Datum] = []

        try:
            df = pd.read_csv(
                miso_path,
                sep="\t",
                comment="#",
                header=0,
            )

            # Standardize column names (case-insensitive, partial match)
            col_map = {}
            for col in df.columns:
                col_lower = col.strip().lower()
                if col_lower == "specie" or col_lower == "species":
                    col_map[col] = "Specie"
                elif col_lower == "fragment":
                    col_map[col] = "Fragment"
                elif col_lower == "dataset":
                    col_map[col] = "Dataset"
                elif col_lower == "isospecies":
                    col_map[col] = "Isospecies"
                elif col_lower == "value":
                    col_map[col] = "Value"
                elif col_lower == "sd" or col_lower == "std":
                    col_map[col] = "SD"
                elif col_lower == "time":
                    col_map[col] = "Time"
            df = df.rename(columns=col_map)

            for _, row in df.iterrows():
                specie = str(row.get("Specie", "")).strip()
                _frag_raw = row.get("Fragment")
                if pd.isna(_frag_raw):
                    fragment = ""
                else:
                    fragment = str(_frag_raw).strip()
                    try:
                        fragment = str(int(float(fragment)))
                    except (ValueError, OverflowError):
                        pass
                dataset = str(row.get("Dataset", "")).strip()
                isospecies = str(row.get("Isospecies", "")).strip()
                value = row.get("Value")
                sd = row.get("SD", 0.01)
                time = row.get("Time")

                if not specie or pd.isna(specie):
                    continue

                # Use dataset as group ID (e.g., "MS-1")
                group_id = (
                    dataset if dataset and not pd.isna(dataset) else specie
                )

                # Collect group info for building expression later
                if group_id not in group_info:
                    group_info[group_id] = {
                        "specie": specie,
                        "fragment": fragment,
                        "isospecies": set(),
                    }

                # Collect isospecies for this group (e.g., M0, M1, M2, ...)
                if isospecies:
                    group_info[group_id]["isospecies"].add(isospecies)

                # Create datum
                if not pd.isna(value):
                    datum_id = (
                        f"{group_id}:{isospecies}" if isospecies else group_id
                    )
                    datum = Datum(
                        id=datum_id,
                        value=float(value),
                        stddev=float(sd) if not pd.isna(sd) else 0.01,
                        time=float(time)
                        if time and not pd.isna(time)
                        else None,
                    )
                    data.append(datum)

            # Second pass: build groups with proper textual expressions
            groups: Dict[str, Group] = {}
            for group_id, info in group_info.items():
                specie = info["specie"]
                fragment = info["fragment"]
                isospecies_set = info["isospecies"]

                # Build atom fragment part
                if fragment:
                    is_cumomer = bool(
                        re.match(r"^[01xX]+$", fragment.replace(",", ""))
                    )
                    if is_cumomer:
                        atom_str = f":{fragment}"
                    else:
                        try:
                            positions = [int(p) for p in fragment.split(",")]
                            positions.sort()
                            if positions == list(
                                range(positions[0], positions[-1] + 1)
                            ):
                                if len(positions) > 2:
                                    atom_str = (
                                        f"[{positions[0]}-{positions[-1]}]"
                                    )
                                else:
                                    atom_str = f"[{
                                        (','.join(str(p) for p in positions))
                                    }]"
                            else:
                                atom_str = (
                                    f"[{','.join(str(p) for p in positions)}]"
                                )
                        except ValueError:
                            atom_str = f"[{fragment}]"
                else:
                    atom_str = ""

                # Build mass isotopomer part
                if isospecies_set:
                    all_mass = all(
                        iso.startswith("M") and iso[1:].isdigit()
                        for iso in isospecies_set
                    )
                    if all_mass:
                        sorted_iso = sorted(
                            isospecies_set,
                            key=lambda x: int(x[1:]),
                        )
                        mass_numbers = [iso[1:] for iso in sorted_iso]
                        mass_str = f"#M{','.join(mass_numbers)}"
                    else:
                        sorted_iso = sorted(isospecies_set)
                        mass_str = f":{','.join(sorted_iso)}"
                else:
                    mass_str = ""

                expression = f"{specie}{atom_str}{mass_str}"

                groups[group_id] = Group(
                    id=group_id,
                    scale="auto",
                    expression=TextualOrMath(textual=expression),
                )

        except Exception as e:
            raise ValueError(f"Error parsing miso file {miso_path}: {e}")

        if not groups:
            return None, None

        return (LabelingMeasurement(groups=list(groups.values())), data)

    def _parse_mflux(
        self, base_path: Path
    ) -> Tuple[Optional[FluxMeasurement], Optional[List[Datum]]]:
        r"""Parse .mflux file containing flux measurements.

        Format (TSV): Id\tComment\tFlux\tValue\tSD
        """
        mflux_path = base_path.with_suffix(".mflux")
        if not mflux_path.exists():
            return None, None

        net_fluxes: List[NetFlux] = []
        data: List[Datum] = []

        try:
            df = pd.read_csv(
                mflux_path,
                sep="\t",
                comment="#",
                header=0,
                names=["Id", "Comment", "Flux", "Value", "SD"],
            )

            for _, row in df.iterrows():
                flux_id = str(row.get("Flux", "")).strip()
                value = row.get("Value")
                sd = row.get("SD", 0.01)

                if not flux_id or pd.isna(flux_id):
                    continue

                # Create net flux measurement
                # ErrorModel uses an expression to define the error
                # For absolute error with SD, use the SD value as the
                # error expression
                sd_value = float(sd) if not pd.isna(sd) else 0.01
                net_flux = NetFlux(
                    id=flux_id,
                    expression=TextualOrMath(textual=flux_id),
                    errormodel=ErrorModel(
                        expression=TextualOrMath(textual=str(sd_value))
                    ),
                )
                net_fluxes.append(net_flux)

                # Create datum
                if not pd.isna(value):
                    datum = Datum(
                        id=flux_id,
                        value=float(value),
                        stddev=float(sd) if not pd.isna(sd) else 0.01,
                    )
                    data.append(datum)

        except Exception as e:
            raise ValueError(f"Error parsing mflux file {mflux_path}: {e}")

        if not net_fluxes:
            return None, None

        return (FluxMeasurement(net_fluxes=net_fluxes), data)

    def _parse_mmet(
        self, base_path: Path
    ) -> Tuple[Optional[MetaboliteSizeMeasurement], Optional[List[Datum]]]:
        r"""Parse .mmet file containing metabolite concentration measurements.

        Format (TSV): Id\tComment\tSpecie\tValue\tSD
        """
        mmet_path = base_path.with_suffix(".mmet")
        if not mmet_path.exists():
            return None, None

        metabolite_sizes: List[MetaboliteSize] = []
        data: List[Datum] = []

        try:
            df = pd.read_csv(
                mmet_path,
                sep="\t",
                comment="#",
                header=0,
            )

            for _, row in df.iterrows():
                specie = str(row.get("Specie", "")).strip()
                value = row.get("Value")
                sd = row.get("SD", 0.01)

                if not specie or pd.isna(specie):
                    continue

                met_size_id = f"ps_{specie}"
                sd_value = float(sd) if not pd.isna(sd) else 0.01

                met_size = MetaboliteSize(
                    id=met_size_id,
                    expression=TextualOrMath(textual=specie),
                    errormodel=ErrorModel(
                        expression=TextualOrMath(textual=str(sd_value))
                    ),
                )
                metabolite_sizes.append(met_size)

                if not pd.isna(value):
                    datum = Datum(
                        id=met_size_id,
                        value=float(value),
                        stddev=sd_value,
                    )
                    data.append(datum)

        except Exception as e:
            raise ValueError(f"Error parsing mmet file {mmet_path}: {e}")

        if not metabolite_sizes:
            return None, None

        return (
            MetaboliteSizeMeasurement(metabolite_sizes=metabolite_sizes),
            data,
        )

    def _parse_variables(self, base_path: Path) -> Optional[Simulation]:
        r"""Parse .tvar file containing variable types and starting values.

        Format (TSV): Id\tComment\tName\tKind\tType\tValue

        - Kind: NET or XCH
        - Type: F (Free), D (Dependent), C (Constrained)
        """
        tvar_path = base_path.with_suffix(".tvar")
        if not tvar_path.exists():
            return None

        flux_values: List[FluxValue] = []
        metab_values: List[MetaboliteSizeValue] = []

        try:
            df = pd.read_csv(
                tvar_path,
                sep="\t",
                comment="#",
                header=0,
                names=["Id", "Comment", "Name", "Kind", "Type", "Value"],
            )

            for _, row in df.iterrows():
                name = str(row.get("Name", "")).strip()
                kind = str(row.get("Kind", "")).strip().upper()
                var_type = str(row.get("Type", "")).strip().upper()
                value = row.get("Value")

                if not name or pd.isna(name):
                    continue

                if kind == "METAB":
                    if var_type == "C" and not pd.isna(value):
                        metab_values.append(
                            MetaboliteSizeValue(
                                metabolite=name,
                                lo=float(value),
                                hi=float(value),
                                value=float(value),
                                type="C",
                            )
                        )
                    elif var_type == "F" and not pd.isna(value):
                        metab_values.append(
                            MetaboliteSizeValue(
                                metabolite=name,
                                lo=0.0,
                                hi=1e6,
                                value=float(value),
                                type="F",
                            )
                        )
                    elif var_type == "D":
                        metab_values.append(
                            MetaboliteSizeValue(
                                metabolite=name,
                                type="D",
                            )
                        )
                    continue

                if kind != "NET":
                    continue

                # Determine bounds based on type
                # F = Free, D = Dependent, C = Constrained
                if var_type == "C" and not pd.isna(value):
                    # Constrained: fixed value
                    flux_value = FluxValue(
                        flux=name,
                        lo=float(value),
                        hi=float(value),
                        value=float(value),
                        type="C",
                    )
                elif var_type == "F" and not pd.isna(value):
                    # Free: starting value, wide bounds
                    flux_value = FluxValue(
                        flux=name,
                        lo=-1000.0,
                        hi=1000.0,
                        value=float(value),
                        type="F",
                    )
                elif var_type == "D":
                    # Dependent: calculated from other fluxes
                    flux_value = FluxValue(
                        flux=name,
                        type="D",
                    )
                else:
                    continue

                flux_values.append(flux_value)

        except Exception as e:
            raise ValueError(f"Error parsing tvar file {tvar_path}: {e}")

        if not flux_values and not metab_values:
            return None

        return Simulation(
            variables=Variables(
                flux_values=flux_values,
                metabolitesize_values=metab_values,
            ),
        )


def parse_mtf(base_path: str) -> FluxomicsData:
    """Convenience function to parse MTF files.

    Args:
        base_path: Path to base filename (without extension) or any MTF file.

    Returns:
        FluxomicsData containing the parsed data
    """
    parser = MTFParser()
    return parser.parse(base_path)
