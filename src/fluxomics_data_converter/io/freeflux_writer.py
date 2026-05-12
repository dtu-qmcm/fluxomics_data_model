"""Freeflux tabular format writer for tsv/csv/xlsx files.

Writes FluxomicsData to Freeflux tabular format files:
    - reactions.tsv: Network definition with atom transitions
    - fluxes.tsv: Flux values (simulation/reference)
    - concentrations.tsv: Metabolite pool sizes
    - measured_MDVs.tsv: Steady-state mass distribution vector measurements
    - measured_fluxes.tsv: Measured flux values with uncertainties
    - measured_inst_MDVs.tsv: Time-course MDV measurements
"""

from pathlib import Path
from typing import Optional, Dict, List, Tuple
from collections import defaultdict

from ..core.core import FluxomicsData, LabelingExperiments
from ..model.atom_mapping import AtomTransition
from ..experiment.measurement import Datum


class FreefluxWriter:
    """Serialiser for Freeflux tabular format files (TSV).

    Writes a :class:`~fluxomics_data_converter.FluxomicsData` to the
    Freeflux directory-based format.  The following files are produced
    when the corresponding data is present:

    - ``reactions.tsv``       — always written
    - ``fluxes.tsv``          — when simulation variables exist
    - ``concentrations.tsv``  — when pool-size variables exist
    - ``measured_MDVs.tsv``   — when steady-state labelling data exists
    - ``measured_inst_MDVs.tsv`` — when time-course labelling data exists
    - ``measured_fluxes.tsv`` — when flux measurement data exists

    Net/exchange fluxes are converted back to the Freeflux ``_f``/``_b``
    (forward/backward) naming convention on write.

    Example::

        writer = FreefluxWriter()
        writer.write(model, "output/ecoli/", experiment_name="default")
        # or use the convenience function:
        from fluxomics_data_converter.io import write_freeflux
        write_freeflux(model, "output/ecoli/", experiment_name="default")
    """

    def write(
        self,
        model: FluxomicsData,
        output_dir: str,
        experiment_name: Optional[str] = None,
    ) -> None:
        """Write FluxomicsData to Freeflux tabular files.

        Args:
            model: FluxomicsData to write
            output_dir: Directory to write files to
            experiment_name: Name of experiment to write (required if
                multiple experiments exist)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        experiment = self._resolve_experiment(model, experiment_name)

        self._write_reactions(model, output_dir)

        if experiment:
            self._write_fluxes(experiment, model, output_dir)
            self._write_concentrations(experiment, output_dir)
            self._write_measurements(experiment, output_dir)

    def _resolve_experiment(
        self,
        model: FluxomicsData,
        experiment_name: Optional[str],
    ) -> Optional[LabelingExperiments]:
        """Resolve which experiment to write."""
        if not model.experiments:
            return None

        if experiment_name is not None:
            exp = model.get_experiments(experiment_name)
            if exp is None:
                raise ValueError(
                    f"Experiment '{experiment_name}' not found. "
                    f"Available: {model.experiments_names}"
                )
            return exp

        if len(model.experiments) == 1:
            return model.experiments[0]

        raise ValueError(
            f"Multiple experiments exist ({list(model.experiments_names)}). "
            f"Please specify experiment_name."
        )

    def _write_reactions(self, model: FluxomicsData, output_dir: Path) -> None:
        """Write reactions.tsv with atom transitions."""
        lines = [
            "#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility"
        ]

        for rxn in model.model.reactions:
            am = model.model.atom_mappings.get(rxn.id)
            if am is not None:
                reactants_str = self._format_side_from_mapping(am, "reactants")
                products_str = self._format_side_from_mapping(am, "products")
            else:
                reactants_str = "+".join(rxn.reactants)
                products_str = "+".join(rxn.products) if rxn.products else ""

            rev = 1 if rxn.reversibility else 0
            lines.append(f"{rxn.id}\t{reactants_str}\t{products_str}\t{rev}")

        filepath = output_dir / "reactions.tsv"
        filepath.write_text("\n".join(lines) + "\n")

    def _format_side_from_mapping(self, am: AtomTransition, side: str) -> str:
        """Format a reaction side from atom transition notation.

        Uses the letter notation directly from the mapping to preserve
        compound multiplicity (e.g., G3P(bdf)+G3P(cea) not
        G3P(bdf,cea)+G3P(bdf,cea)) and variant notation (e.g.,
        Suc(abcd,dcba) for symmetric compounds).

        Strategy:
        - Parse each map's side into ordered (compound, atoms) pairs
        - If all maps produce the same sequence of (compound, atoms) pairs,
          write them directly with multiplicities
        - If maps differ (variant reactions), group per-position variants
          using comma-separated atom strings
        """
        import re

        all_parsed: List[List[tuple]] = []
        for map_id in am.maps:
            try:
                notation = am.to_letter_notation(map_id)
                if notation is None:
                    continue
            except Exception:
                continue

            arrow = "->" if " -> " in notation else "=>"
            sides = notation.split(f" {arrow} ")
            if side == "reactants":
                target = sides[0]
            else:
                target = sides[1] if len(sides) > 1 else ""

            entries = re.findall(r"(\S+?)\(([^)]+)\)", target)
            if entries:
                all_parsed.append(entries)

        if not all_parsed:
            return ""

        n_positions = len(all_parsed[0])

        if all(len(p) == n_positions for p in all_parsed):
            result = []
            for pos in range(n_positions):
                cpd = all_parsed[0][pos][0]
                atom_variants = list(
                    dict.fromkeys(p[pos][1] for p in all_parsed)
                )
                if len(atom_variants) == 1:
                    result.append(f"{cpd}({atom_variants[0]})")
                else:
                    result.append(f"{cpd}({','.join(atom_variants)})")
            return "+".join(result)
        else:
            compound_atoms: Dict[str, List[str]] = defaultdict(list)
            for parsed in all_parsed:
                for cpd, atoms in parsed:
                    compound_atoms[cpd].append(atoms)

            parts = []
            for cpd, atom_strs in compound_atoms.items():
                unique_variants = list(dict.fromkeys(atom_strs))
                if len(unique_variants) == 1:
                    parts.append(f"{cpd}({unique_variants[0]})")
                else:
                    parts.append(f"{cpd}({','.join(unique_variants)})")
            return "+".join(parts)

    def _write_fluxes(
        self,
        experiment: LabelingExperiments,
        model: FluxomicsData,
        output_dir: Path,
    ) -> None:
        """Write fluxes.tsv from simulation variables.

        Converts 13CFlux2 net/xch convention to FreeFlux forward/backward:
            net = forward - backward
            xch = min(forward, backward)
        So:
            if net >= 0: forward = net + xch, backward = xch
            if net < 0:  forward = xch,        backward = xch - net
        """
        if not experiment.simulation or not experiment.simulation.variables:
            return

        flux_values = experiment.simulation.variables.flux_values
        if not flux_values:
            return

        computational_to_base = self._build_computational_to_base_map(model)

        grouped: Dict[str, Dict[str, float]] = defaultdict(dict)
        for fv in flux_values:
            flux_id = fv.flux
            value = fv.value
            if value is None:
                if fv.lo is not None and fv.hi is not None:
                    value = (fv.lo + fv.hi) / 2
                elif fv.lo is not None:
                    value = fv.lo
                else:
                    continue

            base_id = computational_to_base.get(flux_id, flux_id)
            grouped[base_id][fv.type] = value

        reversible_ids = set()
        for r in model.model.reactions:
            if r.reversibility:
                reversible_ids.add(r.id)

        lines = ["#flux_ID\tvalue"]
        for base_id in sorted(grouped.keys()):
            types = grouped[base_id]
            is_reversible = base_id in reversible_ids

            if is_reversible and "net" in types and "xch" in types:
                net = types["net"]
                xch = types["xch"]
                if net >= 0:
                    forward = net + xch
                    backward = xch
                else:
                    forward = xch
                    backward = xch - net
                lines.append(f"{base_id}_f\t{forward}")
                lines.append(f"{base_id}_b\t{backward}")
            elif is_reversible and "xch" in types:
                xch = types["xch"]
                net = types.get("net", 0.0)
                if net >= 0:
                    forward = net + xch
                    backward = xch
                else:
                    forward = xch
                    backward = xch - net
                lines.append(f"{base_id}_f\t{forward}")
                lines.append(f"{base_id}_b\t{backward}")
            elif is_reversible and "net" in types:
                net = types["net"]
                if net >= 0:
                    lines.append(f"{base_id}_f\t{net}")
                    lines.append(f"{base_id}_b\t0.0")
                else:
                    lines.append(f"{base_id}_f\t0.0")
                    lines.append(f"{base_id}_b\t{-net}")
            else:
                value = types.get("net", types.get("xch", 0.0))
                lines.append(f"{base_id}\t{value}")

        filepath = output_dir / "fluxes.tsv"
        filepath.write_text("\n".join(lines) + "\n")

    def _build_computational_to_base_map(
        self, model: FluxomicsData
    ) -> Dict[str, str]:
        """Map computational IDs (e.g., V4___1) back to base IDs (V4)."""
        mapping = {}
        for rxn in model.model.reactions:
            if rxn.atom_transition_ids:
                for comp_id in rxn.atom_transition_ids:
                    mapping[comp_id] = rxn.id
            else:
                mapping[rxn.id] = rxn.id
        return mapping

    def _write_concentrations(
        self, experiment: LabelingExperiments, output_dir: Path
    ) -> None:
        """Write concentrations.tsv from simulation variables."""
        if not experiment.simulation or not experiment.simulation.variables:
            return

        metab_values = experiment.simulation.variables.metabolitesize_values
        if not metab_values:
            return

        lines = ["#metab_ID\tvalue"]
        for mv in metab_values:
            value = mv.value
            if value is None:
                if mv.lo is not None and mv.hi is not None:
                    value = (mv.lo + mv.hi) / 2
                elif mv.lo is not None:
                    value = mv.lo
                else:
                    continue
            lines.append(f"{mv.metabolite}\t{value}")

        filepath = output_dir / "concentrations.tsv"
        filepath.write_text("\n".join(lines) + "\n")

    def _write_measurements(
        self, experiment: LabelingExperiments, output_dir: Path
    ) -> None:
        """Write measurement files (measured_MDVs, measured_fluxes, etc.)."""
        if not experiment.measurement:
            return

        measurement = experiment.measurement
        self._write_measured_mdvs(measurement, output_dir)
        self._write_measured_fluxes(measurement, output_dir)
        self._write_measured_inst_mdvs(measurement, output_dir)

    def _write_measured_mdvs(self, measurement, output_dir: Path) -> None:
        """Write measured_MDVs.tsv for steady-state MDV measurements."""
        if not measurement.model.labeling_measurement:
            return

        groups = measurement.model.labeling_measurement.groups
        data = measurement.data.data

        fragment_data: Dict[str, Tuple[list, list]] = {}
        for group in groups:
            fragment_id = group.id
            group_data = [
                d
                for d in data
                if d.id.startswith(fragment_id + ":") or d.id == fragment_id
            ]
            steady_data = [d for d in group_data if d.time is None]
            steady_data.sort(key=lambda d: (d.pos if d.pos is not None else 0))

            if not steady_data:
                continue

            means = [d.value for d in steady_data]
            sds = [d.stddev for d in steady_data]
            fragment_data[fragment_id] = (means, sds)

        if not fragment_data:
            return

        lines = ["#fragment_ID\tmean\tsd"]
        for frag_id, (means, sds) in fragment_data.items():
            mean_str = ",".join(f"{v:.10g}" for v in means)
            sd_str = ",".join(f"{v}" for v in sds)
            lines.append(f"{frag_id}\t{mean_str}\t{sd_str}")

        filepath = output_dir / "measured_MDVs.tsv"
        filepath.write_text("\n".join(lines) + "\n")

    def _write_measured_inst_mdvs(self, measurement, output_dir: Path) -> None:
        """Write measured_inst_MDVs.tsv for time-course MDV measurements."""
        if not measurement.model.labeling_measurement:
            return

        groups = measurement.model.labeling_measurement.groups
        data = measurement.data.data

        inst_groups = [g for g in groups if g.times is not None]
        if not inst_groups:
            return

        lines = ["#fragment_ID\ttime\tmean\tsd"]

        for group in inst_groups:
            fragment_id = group.id
            group_data = [
                d
                for d in data
                if d.id.startswith(fragment_id + ":") or d.id == fragment_id
            ]
            inst_data = [d for d in group_data if d.time is not None]

            time_groups: Dict[float, List[Datum]] = defaultdict(list)
            for d in inst_data:
                time_groups[d.time].append(d)

            for time_val in sorted(time_groups.keys()):
                time_data = sorted(
                    time_groups[time_val], key=lambda d: d.pos or 0
                )
                means = [d.value for d in time_data]
                sds = [d.stddev for d in time_data]
                mean_str = ",".join(f"{v:.10g}" for v in means)
                sd_str = ",".join(f"{v}" for v in sds)
                lines.append(
                    f'{fragment_id}\t{time_val}\t"{mean_str}"\t"{sd_str}"'
                )

        filepath = output_dir / "measured_inst_MDVs.tsv"
        filepath.write_text("\n".join(lines) + "\n")

    def _write_measured_fluxes(self, measurement, output_dir: Path) -> None:
        """Write measured_fluxes.tsv."""
        if not measurement.model.flux_measurement:
            return

        net_fluxes = measurement.model.flux_measurement.net_fluxes
        if not net_fluxes:
            return

        data = measurement.data.data

        lines = ["#reaction_ID\tmean\tsd"]
        for nf in net_fluxes:
            flux_data = [d for d in data if d.id == nf.id]
            if flux_data:
                d = flux_data[0]
                lines.append(f"{nf.id}\t{d.value}\t{d.stddev}")
            else:
                expr = (
                    nf.expression.textual
                    if nf.expression and nf.expression.textual
                    else nf.id
                )
                lines.append(f"{nf.id}\t{expr}\t0.01")

        filepath = output_dir / "measured_fluxes.tsv"
        filepath.write_text("\n".join(lines) + "\n")


def write_freeflux(
    model: FluxomicsData,
    output_dir: str,
    experiment_name: Optional[str] = None,
) -> None:
    """Convenience function to write FluxomicsData to Freeflux files.

    Args:
        model: FluxomicsData to write
        output_dir: Directory to write files to
        experiment_name: Name of experiment to write (optional)
    """
    writer = FreefluxWriter()
    writer.write(model, output_dir, experiment_name=experiment_name)
