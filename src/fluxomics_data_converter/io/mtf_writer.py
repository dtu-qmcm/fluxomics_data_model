"""MTF (Metabolic Text Format) writer for influx_si format.

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

from pathlib import Path
from typing import Optional, Dict, List, Tuple
from collections import defaultdict

from ..core.core import FluxomicsData, LabelingExperiments
from ..model.atom_mapping import AtomTransition


class MTFWriter:
    """Serialiser for influx_si MTF (Metabolic Text Format) files.

    Writes a :class:`~fluxomics_data_converter.FluxomicsData` to the
    multi-file MTF format.  One experiment is written per call; when the
    model contains more than one experiment ``experiment_name`` must be
    provided.

    Known limitations
    -----------------
    - Reactions are grouped under a single ``"Other"`` pathway heading
      because ``Reaction`` does not carry a ``pathway`` attribute.
    - Irreversible reactions are written with the ``"->>"`  arrow (net ≥ 0),
      not ``"->"`` (unconstrained net).
    - ``.tvar`` always writes ``kind = NET``; XCH kind is not preserved.

    Example::

        writer = MTFWriter()
        writer.write(model, "output/ecoli", experiment_name="default")
        # or use the convenience function:
        from fluxomics_data_converter.io import write_mtf
        write_mtf(model, "output/ecoli", experiment_name="default")
    """

    def __init__(self):
        """Initialise the writer with empty state; call write() to serialise."""
        self._model: Optional[FluxomicsData] = None
        self._experiment: Optional[LabelingExperiments] = None
        self._metabolite_atom_counts: Dict[str, int] = {}

    def _compute_metabolite_atom_counts(self) -> Dict[str, int]:
        """Infer atom counts for each metabolite from atom transitions."""
        counts: Dict[str, int] = {}
        for atom_mapping in self._model.model.atom_mappings.values():
            if not atom_mapping.maps:
                continue
            atom_map = next(iter(atom_mapping.maps.values()))
            for src_addr in atom_map.mapping.values():
                counts[src_addr.mol] = max(
                    counts.get(src_addr.mol, 0), src_addr.index
                )
            for prod_addr in atom_map.mapping.keys():
                counts[prod_addr.mol] = max(
                    counts.get(prod_addr.mol, 0), prod_addr.index
                )
        return counts

    def write(
        self,
        model: FluxomicsData,
        base_path: str,
        experiment_name: Optional[str] = None,
    ) -> None:
        """Write FluxomicsData to MTF files.

        Args:
            model: FluxomicsData to write
            base_path: Path to base filename (without extension).
                       For example: "output/e_coli" will create
                       "output/e_coli.netw", "output/e_coli.linp", etc.
            experiment_name: Name of experiment to write. If None and model
                            has a single experiment, uses that. If None and
                            model has multiple experiments, raises ValueError.
        """
        self._model = model
        base_path = Path(base_path)
        self._metabolite_atom_counts = self._compute_metabolite_atom_counts()

        # Select experiment
        if experiment_name:
            self._experiment = model.get_experiments(experiment_name)
            if not self._experiment:
                raise ValueError(f"Experiment '{experiment_name}' not found")
        elif len(model.experiments) == 1:
            self._experiment = model.experiments[0]
        elif len(model.experiments) > 1:
            raise ValueError(
                "Multiple experiments exist. Please specify experiment_name."
            )
        else:
            self._experiment = None

        # Ensure output directory exists
        base_path.parent.mkdir(parents=True, exist_ok=True)

        # Write each file
        self._write_netw(base_path.with_suffix(".netw"))

        if self._experiment and self._experiment.tracers:
            self._write_linp(base_path.with_suffix(".linp"))

        if self._experiment and self._experiment.measurement:
            if self._experiment.measurement.model.labeling_measurement:
                self._write_miso(base_path.with_suffix(".miso"))
            if self._experiment.measurement.model.flux_measurement:
                self._write_mflux(base_path.with_suffix(".mflux"))

        if model.constraints or (
            self._experiment and self._experiment.constraints
        ):
            self._write_cnstr(base_path.with_suffix(".cnstr"))

        if self._experiment and self._experiment.simulation:
            self._write_tvar(base_path.with_suffix(".tvar"))

        self._write_opt(base_path.with_suffix(".opt"))

    def _write_opt(self, filepath: Path) -> None:
        """Write .opt file with influx_si command-line options.

        --clownr and --zc handle near-zero fluxes at the starting point
        (e.g. when the stored simulation values are not mass-balanced for
        the chosen F/D parameterisation).  These are standard options used
        by influx_si reference examples.
        """
        if filepath.exists():
            return
        lines = [
            "Id\tComment\tName\tValue",
            "\t\tcommandArgs\t--clownr 1.e-3 --zc 1.e-4",
        ]
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _write_netw(self, filepath: Path) -> None:
        """Write .netw file containing reaction network with atom transitions.

        Format: reaction_id: substrate (ATOMS) + substrate (atoms) ->
        product (ATOMS)
        """
        lines = []
        lines.append("# Network definition")
        lines.append("# Generated by fluxomics_data_converter")
        lines.append("#")

        # Group reactions by pathway if available
        reactions_by_pathway: Dict[str, List] = defaultdict(list)
        for reaction in self._model.model.reactions:
            pathway = getattr(reaction, "pathway", None) or "Other"
            reactions_by_pathway[pathway].append(reaction)

        for pathway, reactions in reactions_by_pathway.items():
            lines.append(f"# {pathway}")
            lines.append("# " + "-" * 60)

            for reaction in reactions:
                line = self._format_reaction(reaction)
                lines.append(line)

            lines.append("#")

        # Add drain reactions for sink metabolites (produced but never consumed).
        # influx_si requires every measured metabolite to have at least one
        # consuming reaction to be treated as "internal" — without a drain,
        # the metabolite is external and its measurements are silently ignored.
        consumed: set = set()
        produced: set = set()
        for rxn in self._model.model.reactions:
            consumed.update(rxn.reactants)
            produced.update(rxn.products)
        input_pools = {
            tracer.metabolite
            for exp in self._model.experiments
            for tracer in exp.tracers
        }
        sink_metabolites = sorted(produced - consumed - input_pools)
        if sink_metabolites:
            lines.append("# Drain reactions (auto-generated for sink metabolites)")
            lines.append("# " + "-" * 60)
            for metab_id in sink_metabolites:
                metab = next(
                    (m for m in self._model.model.metabolites if m.id == metab_id),
                    None,
                )
                atom_count = (
                    (metab.atoms if metab else None)
                    or self._metabolite_atom_counts.get(metab_id)
                    or 0
                )
                if atom_count:
                    atoms = "".join(
                        chr(ord("a") + i) for i in range(min(atom_count, 26))
                    )
                    # influx_si requires balanced atom transitions: add an external
                    # product ({metab_id}_ext) so both sides carry the same atoms.
                    # {metab_id}_ext only appears as a product → treated as external.
                    lines.append(
                        f"{metab_id}_out: {metab_id} ({atoms}) ->> "
                        f"{metab_id}_ext ({atoms})"
                    )
                else:
                    lines.append(f"{metab_id}_out: {metab_id} ->>")
            lines.append("#")

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _format_reaction(self, reaction) -> str:
        """Format a single reaction line for .netw file."""
        rxn_id = reaction.id
        reversible = reaction.reversibility

        # Get atom transition if exists
        atom_mapping = self._model.model.atom_mappings.get(rxn_id)

        # Build reactant and product strings
        reactants_str = self._format_compounds(
            reaction.reactants, atom_mapping, is_reactant=True
        )
        products_str = self._format_compounds(
            reaction.products, atom_mapping, is_reactant=False
        )

        # Build arrow based on reversibility
        arrow = "<->" if reversible else "->>"

        # Handle reactions with no products (drains)
        if not reaction.products:
            return f"{rxn_id}: {reactants_str} {arrow}"

        return f"{rxn_id}: {reactants_str} {arrow} {products_str}"

    def _format_compounds(
        self,
        compounds: List[str],
        atom_mapping: Optional[AtomTransition],
        is_reactant: bool,
    ) -> str:
        """Format compounds with atom transitions for .netw file."""
        if not compounds:
            return ""

        parts = []
        compound_counts: Dict[str, int] = {}

        for cpd_id in compounds:
            # Track instance for repeated compounds
            compound_counts[cpd_id] = compound_counts.get(cpd_id, 0) + 1
            instance = compound_counts[cpd_id]

            # Get atom letters if atom transition exists
            atoms_str = ""
            if atom_mapping and atom_mapping.maps:
                # Get the first (or only) atom map
                atom_map = next(iter(atom_mapping.maps.values()))

                if is_reactant:
                    # For reactants, collect atoms that are sources
                    atoms = self._get_reactant_atoms(atom_map, cpd_id, instance)
                else:
                    # For products, collect atoms that are targets
                    atoms = self._get_product_atoms(atom_map, cpd_id, instance)

                if atoms:
                    atoms_str = f" ({atoms})"

            parts.append(f"{cpd_id}{atoms_str}")

        return " + ".join(parts)

    def _get_reactant_atoms(self, atom_map, cpd_id: str, instance: int) -> str:
        """Get atom letter notation for a reactant compound."""
        # Build a mapping from (cpd, instance, atom_idx) -> letter
        # based on the product side mappings

        # First, collect all reactant atoms used in products
        reactant_atoms: Dict[Tuple[str, int, int], str] = {}
        letters = (
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        )
        letter_idx = 0

        # Iterate through all product atoms in order
        sorted_product_atoms = sorted(
            atom_map.mapping.keys(), key=lambda a: (a.mol, a.instance, a.index)
        )

        # Assign letters based on first appearance in products
        for prod_addr in sorted_product_atoms:
            src_addr = atom_map.mapping[prod_addr]
            key = (src_addr.mol, src_addr.instance, src_addr.index)
            if key not in reactant_atoms:
                if letter_idx < len(letters):
                    reactant_atoms[key] = letters[letter_idx]
                    letter_idx += 1

        # Now collect atoms for this specific reactant
        atoms = []
        for (mol, inst, idx), letter in sorted(
            reactant_atoms.items(), key=lambda x: x[0][2]
        ):
            if mol == cpd_id and inst == instance:
                atoms.append((idx, letter))

        # Sort by atom index and return letters
        atoms.sort(key=lambda x: x[0])
        return "".join(letter for _, letter in atoms)

    def _get_product_atoms(self, atom_map, cpd_id: str, instance: int) -> str:
        """Get atom letter notation for a product compound."""
        # Same letter assignment as reactants
        reactant_atoms: Dict[Tuple[str, int, int], str] = {}
        letters = (
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        )
        letter_idx = 0

        sorted_product_atoms = sorted(
            atom_map.mapping.keys(), key=lambda a: (a.mol, a.instance, a.index)
        )

        for prod_addr in sorted_product_atoms:
            src_addr = atom_map.mapping[prod_addr]
            key = (src_addr.mol, src_addr.instance, src_addr.index)
            if key not in reactant_atoms:
                if letter_idx < len(letters):
                    reactant_atoms[key] = letters[letter_idx]
                    letter_idx += 1

        # Now collect atoms for this specific product in product order
        atoms = []
        for prod_addr in sorted_product_atoms:
            if prod_addr.mol == cpd_id and prod_addr.instance == instance:
                src_addr = atom_map.mapping[prod_addr]
                key = (src_addr.mol, src_addr.instance, src_addr.index)
                letter = reactant_atoms.get(key, "?")
                atoms.append((prod_addr.index, letter))

        # Sort by atom index and return letters
        atoms.sort(key=lambda x: x[0])
        return "".join(letter for _, letter in atoms)

    def _write_linp(self, filepath: Path) -> None:
        r"""Write .linp file containing tracer specifications.

        Format (TSV): Id\tComment\tSpecie\tIsotopomer\tValue
        """
        lines = []
        lines.append("Id\tComment\tSpecie\tIsotopomer\tValue")

        if self._experiment and self._experiment.tracers:
            for tracer in self._experiment.tracers:
                metabolite = tracer.metabolite
                for label in tracer.labels:
                    isotopomer = label.labeled_pattern
                    fraction = label.fraction
                    lines.append(f"\t\t{metabolite}\t{isotopomer}\t{fraction}")

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _write_miso(self, filepath: Path) -> None:
        r"""Write .miso file containing MS isotopomer measurements.

        Format (TSV):
        Id\tComment\tSpecie\tFragment\tDataset\tIsospecies\tValue\tSD\tTime
        """
        lines = []
        lines.append(
            "Id\tComment\tSpecie\tFragment\tDataset\tIsospecies\tValue\tSD\tTime"
        )

        if (
            self._experiment
            and self._experiment.measurement
            and self._experiment.measurement.model.labeling_measurement
        ):
            labeling = self._experiment.measurement.model.labeling_measurement
            data = (
                self._experiment.measurement.data.data
                if self._experiment.measurement.data
                else []
            )

            # Build a lookup from datum id to datum
            {d.id: d for d in data}

            for group in labeling.groups:
                # Parse group expression to get metabolite and fragment
                expression = (
                    group.expression.textual if group.expression else group.id
                )

                # Parse expression like "Ala[2,3]" or "AKG"
                metabolite, fragment = self._parse_ms_expression(expression)

                # Find all data points for this group
                # Group ID format varies, so we search for matching data
                for datum in data:
                    if datum.id.startswith(group.id):
                        # Extract isospecies from datum id
                        # (e.g., "group:M0" -> "M0")
                        parts = datum.id.split(":")
                        isospecies = parts[-1] if len(parts) > 1 else ""

                        time_str = (
                            str(datum.time) if datum.time is not None else ""
                        )

                        lines.append(
                            f"\t\t{metabolite}\t{fragment}\t{group.id}\t{isospecies}\t"
                            f"{datum.value}\t{datum.stddev}\t{time_str}"
                        )

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _parse_ms_expression(self, expression: str) -> Tuple[str, str]:
        """Parse MS expression like 'Ala[2,3]' or 'F[1,2,3]#M0,1,2,3' into (metabolite, fragment).

        Strips the optional '#M...' MDV-weight suffix used by x3cflux / 13CFlux2
        before parsing, so both 'F[1,2,3]' and 'F[1,2,3]#M0,1,2,3' are handled.
        """
        import re

        # Strip optional x3cflux MDV-weight suffix: 'F[1,2,3]#M0,1,2,3' -> 'F[1,2,3]'
        expression = expression.split("#")[0].strip()

        # Match patterns like "Ala[2,3]" or "Ala[2+3]" or "AKG"
        match = re.match(r"^(\w+)(?:\[([^\]]+)\])?$", expression)
        if match:
            metabolite = match.group(1)
            fragment_spec = match.group(2)

            if fragment_spec:
                # Expand range notation "2-4" -> "2,3,4"; normalise "2+3" -> "2,3"
                fragment = self._expand_fragment_spec(fragment_spec)
            else:
                # Full molecule — derive from atom count if available
                try:
                    metabolite_obj = self._model.model.metabolites[metabolite]
                    atom_count = metabolite_obj.atoms
                except (KeyError, IndexError):
                    atom_count = None
                if atom_count:
                    fragment = ",".join(str(i) for i in range(1, atom_count + 1))
                else:
                    fragment = ""

            return metabolite, fragment

        return expression, ""

    @staticmethod
    def _expand_fragment_spec(spec: str) -> str:
        """Expand a fragment spec into a comma-separated position list.

        Handles three notations:
        - ``"2,3,4"``  → ``"2,3,4"``  (passthrough)
        - ``"2+3+4"``  → ``"2,3,4"``  (+ → ,)
        - ``"2-4"``    → ``"2,3,4"``  (range expansion)
        - ``"1,3-5,7"``→ ``"1,3,4,5,7"``  (mixed)
        """
        import re

        positions: list[int] = []
        # Split on commas or '+' first
        for token in re.split(r"[,+]", spec):
            token = token.strip()
            if not token:
                continue
            range_match = re.match(r"^(\d+)-(\d+)$", token)
            if range_match:
                lo, hi = int(range_match.group(1)), int(range_match.group(2))
                positions.extend(range(lo, hi + 1))
            else:
                try:
                    positions.append(int(token))
                except ValueError:
                    pass

        return ",".join(str(p) for p in positions) if positions else spec

    def _write_mflux(self, filepath: Path) -> None:
        r"""Write .mflux file containing flux measurements.

        Format (TSV): Id\tComment\tFlux\tValue\tSD
        """
        lines = []
        lines.append("Id\tComment\tFlux\tValue\tSD")

        if (
            self._experiment
            and self._experiment.measurement
            and self._experiment.measurement.model.flux_measurement
        ):
            flux_meas = self._experiment.measurement.model.flux_measurement
            data = (
                self._experiment.measurement.data.data
                if self._experiment.measurement.data
                else []
            )

            # Build a lookup from datum id to datum
            datum_lookup = {d.id: d for d in data}

            for net_flux in flux_meas.net_fluxes:
                flux_id = net_flux.id
                expression = (
                    net_flux.expression.textual
                    if net_flux.expression
                    else flux_id
                )

                # Find matching datum
                datum = datum_lookup.get(flux_id)
                if datum:
                    lines.append(
                        f"\t\t{expression}\t{datum.value}\t{datum.stddev}"
                    )

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _write_cnstr(self, filepath: Path) -> None:
        r"""Write .cnstr file containing constraints.

        Format (TSV): Id\tComment\tKind\tFormula\tOperator\tValue
        """
        lines = []
        lines.append("Id\tComment\tKind\tFormula\tOperator\tValue")

        # Collect constraints from model level
        constraints = self._model.constraints

        # Also check experiment-level constraints
        if self._experiment and self._experiment.constraints:
            if constraints:
                # Merge constraints (experiment takes precedence)
                pass  # For now, just use model constraints
            else:
                constraints = self._experiment.constraints

        if constraints:
            # Write NET constraints
            if constraints.net and constraints.net.formulas:
                for formula in constraints.net.formulas:
                    expression = formula.expression
                    # Parse expression to extract operator and value
                    op, lhs, rhs = self._parse_constraint_expression(expression)
                    lines.append(f"\t\tNET\t{lhs}\t{op}\t{rhs}")

            # Write XCH constraints
            if constraints.xch and constraints.xch.formulas:
                for formula in constraints.xch.formulas:
                    expression = formula.expression
                    op, lhs, rhs = self._parse_constraint_expression(expression)
                    lines.append(f"\t\tXCH\t{lhs}\t{op}\t{rhs}")

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _parse_constraint_expression(
        self, expression: str
    ) -> Tuple[str, str, str]:
        """Parse constraint expression into (operator, lhs, rhs)."""
        import re

        # Try to match patterns like "A = B" or "A == B" or "A >= B" or "A <= B"
        match = re.match(r"^(.+?)\s*(==|>=|<=|=)\s*(.+)$", expression)
        if match:
            lhs = match.group(1).strip()
            op = match.group(2)
            rhs = match.group(3).strip()
            return op, lhs, rhs

        # No operator found - return as formula with ==0
        return "==", expression, "0"

    def _compute_dependent_fluxes(self) -> set:
        """Determine which net fluxes should be Dependent (D) in .tvar.

        Strategy: reactions with zero stored simulation value are D (they are
        computed from mass balance); reactions with non-zero stored values are F
        (they are the free parameters, matching the original model's intent).
        Drain reactions are always D.  If the zero-valued reactions do not form
        a valid (non-singular) basis, QR with column pivoting is used as a fallback.
        """
        import numpy as np
        from scipy.linalg import qr

        reactions = list(self._model.model.reactions)
        n_rxn = len(reactions)
        rxn_names = [r.id for r in reactions]
        metabolites = list(self._model.model.metabolites)
        met_idx = {m.id: i for i, m in enumerate(metabolites)}

        all_consumed: set = set()
        all_produced: set = set()
        for rxn in reactions:
            all_consumed.update(rxn.reactants)
            all_produced.update(rxn.products)
        input_pools = {
            t.metabolite
            for exp in self._model.experiments
            for t in exp.tracers
        }
        sink_mets = sorted(all_produced - all_consumed - input_pools)
        drain_ids = [f"{m}_out" for m in sink_mets]
        all_rxn_names = rxn_names + drain_ids
        n_all = len(all_rxn_names)

        # Build stoichiometric matrix (internal metabolites only)
        S = np.zeros((len(metabolites), n_all))
        for j, rxn in enumerate(reactions):
            for met_id in rxn.reactants:
                if met_id in met_idx:
                    S[met_idx[met_id], j] -= 1.0
            for met_id in rxn.products:
                if met_id in met_idx:
                    S[met_idx[met_id], j] += 1.0
        for k, met_id in enumerate(sink_mets):
            if met_id in met_idx:
                S[met_idx[met_id], n_rxn + k] -= 1.0

        int_row_indices = [
            i for i, m in enumerate(metabolites) if m.id not in input_pools
        ]
        S_int = S[int_row_indices, :]
        n_int = len(int_row_indices)

        # Stored simulation values: zero-valued reactions are candidates for D
        stored_vals: Dict[str, float] = {}
        if self._experiment and self._experiment.simulation:
            for fv in self._experiment.simulation.variables.flux_values:
                stored_vals[fv.flux] = abs(fv.value or 0.0)

        # D candidates: zero-valued reactions + drain reactions
        threshold = 1e-9
        dependent: set = set(drain_ids)
        for rxn_name in rxn_names:
            if stored_vals.get(rxn_name, 0.0) <= threshold:
                dependent.add(rxn_name)

        # Validate: the D columns must span rank(S_int) rows
        d_cols = [i for i, n in enumerate(all_rxn_names) if n in dependent]
        rank_needed = int(np.linalg.matrix_rank(S_int)) if n_int > 0 else 0

        if len(d_cols) == 0 or int(np.linalg.matrix_rank(S_int[:, d_cols])) < rank_needed:
            # Fallback: QR on uncovered rows to supplement D
            covered: set = set()
            for row_i in range(n_int):
                if any(S_int[row_i, c] != 0 for c in d_cols):
                    covered.add(row_i)
            uncovered = [i for i in range(n_int) if i not in covered]
            non_d_cols = [i for i in range(n_all) if i not in set(d_cols)]
            if uncovered and non_d_cols:
                S_uncov = S_int[np.ix_(uncovered, non_d_cols)]
                try:
                    _, _, perm = qr(S_uncov, pivoting=True)
                    r_uncov = int(np.linalg.matrix_rank(S_uncov))
                    for p in perm[:min(r_uncov, len(uncovered))]:
                        dependent.add(all_rxn_names[non_d_cols[p]])
                except Exception:
                    pass

        return dependent

    def _solve_balanced_starting_values(
        self,
        dependent_ids: set,
        flux_values: list,
        drain_ids: list,
    ) -> Dict[str, float]:
        """Compute mass-balanced D starting values from stored F values.

        Solves S[:,D]*v_D = -S[:,F]*v_F (S*v=0 with F values from simulation).
        Any near-zero fluxes at the starting point are handled by influx_si's
        --clownr option written to the .opt file.
        """
        import numpy as np

        reactions = list(self._model.model.reactions)
        n_rxn = len(reactions)
        rxn_names = [r.id for r in reactions]
        metabolites = list(self._model.model.metabolites)
        met_idx = {m.id: i for i, m in enumerate(metabolites)}

        input_pools = {
            t.metabolite
            for exp in self._model.experiments
            for t in exp.tracers
        }
        sink_mets = sorted(
            {m for rxn in reactions for m in rxn.products}
            - {m for rxn in reactions for m in rxn.reactants}
            - input_pools
        )

        all_rxn_names = rxn_names + drain_ids
        n_all = len(all_rxn_names)

        S = np.zeros((len(metabolites), n_all))
        for j, rxn in enumerate(reactions):
            for met_id in rxn.reactants:
                if met_id in met_idx:
                    S[met_idx[met_id], j] -= 1.0
            for met_id in rxn.products:
                if met_id in met_idx:
                    S[met_idx[met_id], j] += 1.0
        for k, met_id in enumerate(sink_mets):
            if met_id in met_idx:
                S[met_idx[met_id], n_rxn + k] -= 1.0

        int_rows = [i for i, m in enumerate(metabolites) if m.id not in input_pools]
        S_int = S[int_rows, :]

        stored = {fv.flux: (fv.value or 0.0) for fv in flux_values}
        d_cols = [i for i, n in enumerate(all_rxn_names) if n in dependent_ids]
        f_cols = [i for i, n in enumerate(all_rxn_names) if n not in dependent_ids]
        d_names = [all_rxn_names[i] for i in d_cols]
        f_names = [all_rxn_names[i] for i in f_cols]

        # Measured flux values — used as fallback anchors
        measured: Dict[str, float] = {}
        if (
            self._experiment
            and self._experiment.measurement
            and self._experiment.measurement.model.flux_measurement
            and self._experiment.measurement.data
        ):
            data_map = {d.id: d.value for d in self._experiment.measurement.data.data}
            for nf in self._experiment.measurement.model.flux_measurement.net_fluxes:
                v = data_map.get(nf.id)
                if v is not None:
                    measured[nf.id] = float(v)

        def _compute_d(f_vals: Dict[str, float]) -> Dict[str, float]:
            v_f = np.array([f_vals.get(n, 0.0) for n in f_names])
            try:
                v_d, *_ = np.linalg.lstsq(
                    S_int[:, d_cols], -S_int[:, f_cols] @ v_f, rcond=None
                )
            except Exception:
                v_d = np.zeros(len(d_cols))
            res = {n: float(f_vals.get(n, 0.0)) for n in f_names}
            res.update({d_names[i]: float(v_d[i]) for i in range(len(d_cols))})
            return res

        result = _compute_d({n: stored.get(n, 0.0) for n in f_names})

        # If any D drain flux is ≤ 0, fall back: measured F values stay,
        # non-measured free fluxes are set to 0 to get feasible D starting values.
        drain_d = [n for n in d_names if n in drain_ids]
        if any(result.get(n, 0.0) <= 1e-9 for n in drain_d):
            f_fallback = {n: measured.get(n, 0.0) for n in f_names}
            result2 = _compute_d(f_fallback)
            if all(result2.get(n, 0.0) > 1e-9 for n in drain_d):
                result = result2

        return result

    def _write_tvar(self, filepath: Path) -> None:
        r"""Write .tvar file containing variable types and starting values.

        Format (TSV): Id\tComment\tName\tKind\tType\tValue

        Type assignment (influx_si only supports F and D; C is reserved for
        fluxes without measurements — measured fluxes must be F or D):
          D — fluxes computed from F fluxes via S*v=0 mass balance
          F — free parameters the optimizer searches over
        """
        lines = []
        lines.append("Id\tComment\tName\tKind\tType\tValue")

        if self._experiment and self._experiment.simulation:
            sim = self._experiment.simulation

            if sim.variables:
                # Drain reactions (auto-generated sinks)
                all_consumed_: set = set()
                all_produced_: set = set()
                for rxn in self._model.model.reactions:
                    all_consumed_.update(rxn.reactants)
                    all_produced_.update(rxn.products)
                input_pools_ = {
                    t.metabolite
                    for exp in self._model.experiments
                    for t in exp.tracers
                }
                sink_ids_ = sorted(all_produced_ - all_consumed_ - input_pools_)
                drain_ids_ = [f"{m}_out" for m in sink_ids_]

                # QR-based F/D assignment from stoichiometric null space
                dependent_ids = self._compute_dependent_fluxes()

                # Compute mass-balanced starting values (solve S[:,D]*v_D = -S[:,F]*v_F)
                balanced = self._solve_balanced_starting_values(
                    dependent_ids,
                    sim.variables.flux_values,
                    drain_ids_,
                )

                # Write flux rows (F or D only — C is not valid for any flux here)
                for flux_val in sim.variables.flux_values:
                    flux_name = flux_val.flux
                    kind = "XCH" if flux_val.type == "xch" else "NET"
                    fdc = "D" if flux_name in dependent_ids else "F"
                    value = balanced.get(flux_name, flux_val.value or 0.0)
                    lines.append(f"\t\t{flux_name}\t{kind}\t{fdc}\t{value}")

                # Drain-reaction entries (always D, value from balanced solve)
                for sink_id in sink_ids_:
                    drain_name = f"{sink_id}_out"
                    value = balanced.get(drain_name, 0.0)
                    lines.append(f"\t\t{drain_name}\tNET\tD\t{value}")

                # Write metabolite size values (if any)
                for met_val in sim.variables.metabolitesize_values:
                    met_name = met_val.metabolite
                    var_type = met_val.type or "F"
                    value = met_val.value if met_val.value is not None else ""

                    lines.append(f"\t\t{met_name}\tMETAB\t{var_type}\t{value}")

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")


def write_mtf(
    model: FluxomicsData,
    base_path: str,
    experiment_name: Optional[str] = None,
) -> None:
    """Convenience function to write FluxomicsData to MTF files.

    Args:
        model: FluxomicsData to write
        base_path: Path to base filename (without extension)
        experiment_name: Name of experiment to write (optional)
    """
    writer = MTFWriter()
    writer.write(model, base_path, experiment_name)
