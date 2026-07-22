"""FluxML XML writer for writing FluxML files.

FluxML is an XML-based format used by 13CFlux/13CFlux2 for 13C metabolic
flux analysis.
"""

import xml.etree.ElementTree as ET
from xml.dom import minidom
from typing import TYPE_CHECKING, Optional, Dict
from datetime import datetime
from pathlib import Path

from ..core.core import FluxomicsData, LabelingExperiments
from ..model.atom_mapping import AtomTransition

if TYPE_CHECKING:
    import numpy as np


class FluxMLWriter:
    """Serialiser for FluxML XML files (13CFlux2 format).

    Converts a :class:`~fluxomics_data_converter.FluxomicsData` into a
    pretty-printed UTF-8 XML file.

    Known limitations
    -----------------
    - Named constraints (``"name: expression"``) lose their name prefix on
      write; only the expression is emitted.
    - Experiment-level ``<constraints>`` sections are not written.
    - ``<variant>`` elements for scrambling reactions do not include their
      ``ratio`` / weight attributes.
    - Tracer ``type`` is always written as ``"isotopomer"``; the actual
      ``Tracers.type`` field value is ignored.

    Example::

        writer = FluxMLWriter()
        writer.write(model, "output/ecoli.fml")
        # or use the convenience function:
        from fluxomics_data_converter.io import write_fluxml
        write_fluxml(model, "output/ecoli.fml")
    """

    NAMESPACE = "http://www.13cflux.net/fluxml"

    def __init__(self):
        """Initialise the writer with empty state; call write() to serialise."""
        self._model: Optional[FluxomicsData] = None
        self._ns = self.NAMESPACE
        self._metabolite_atom_counts: Dict[str, int] = {}

    def _compute_metabolite_atom_counts(self) -> Dict[str, int]:
        """Compute atom counts for metabolites by scanning atom transitions.

        Returns a dict mapping metabolite id to atom count.
        """
        atom_counts: Dict[str, int] = {}

        for (
            reaction_id,
            atom_mapping,
        ) in self._model.model.atom_mappings.items():
            if not atom_mapping.maps:
                continue

            # Use the first atom map to get atom counts
            atom_map = next(iter(atom_mapping.maps.values()))

            # Scan reactant atoms (from source addresses)
            for src_addr in atom_map.mapping.values():
                mol_id = src_addr.mol
                atom_idx = src_addr.index
                if mol_id not in atom_counts:
                    atom_counts[mol_id] = atom_idx
                else:
                    atom_counts[mol_id] = max(atom_counts[mol_id], atom_idx)

            # Scan product atoms (from destination addresses)
            for prod_addr in atom_map.mapping.keys():
                mol_id = prod_addr.mol
                atom_idx = prod_addr.index
                if mol_id not in atom_counts:
                    atom_counts[mol_id] = atom_idx
                else:
                    atom_counts[mol_id] = max(atom_counts[mol_id], atom_idx)

        return atom_counts

    # Natural abundance of 13C per carbon position, as used by x3cflux internally.
    # Source: x3cflux/src/main/cpp/model/system/NaturalLabelingInitializer.cpp
    NATURAL_ABUNDANCE_CARBON: float = 0.01109

    def write(
        self,
        model: FluxomicsData,
        filepath: str,
        apply_na_correction: bool = False,
    ) -> None:
        """Write FluxomicsData to a FluxML file.

        Args:
            model: FluxomicsData to write
            filepath: Path to output file
            apply_na_correction: If True, apply the ¹³C natural abundance
                correction matrix M_C (p=0.01109/C) to all MS measurement
                data before writing.  Set this when the data is synthetic
                (e.g. from cmfa) and does NOT already contain natural
                abundance, so that x3cflux can recover unbiased flux
                estimates.  Leave False when the data comes from real MS
                experiments (natural abundance is already present in the
                raw signal).

                Note: M_C is a binomial approximation.  The exact correction
                is network-dependent; expect a residual of ~10⁻³ MID units.
        """
        self._model = model
        self._apply_na_correction = apply_na_correction
        filepath = Path(filepath)

        # Ensure output directory exists
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # Build XML tree
        root = self._build_fluxml()

        # Write to file with pretty formatting
        xml_str = ET.tostring(root, encoding="unicode")
        # Pretty print using minidom
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ", encoding="UTF-8")

        # Remove extra blank lines that minidom adds
        lines = pretty_xml.decode("UTF-8").split("\n")
        cleaned_lines = [line for line in lines if line.strip()]
        # Re-add proper XML declaration
        if cleaned_lines and cleaned_lines[0].startswith("<?xml"):
            cleaned_lines[0] = (
                '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
            )

        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(cleaned_lines))

    def _build_fluxml(self) -> ET.Element:
        """Build the root fluxml element."""
        # Register namespace - use empty prefix for default namespace
        ET.register_namespace("", self.NAMESPACE)

        # Create root element with namespace (don't use xmlns
        # attribute separately)
        root = ET.Element("fluxml")
        root.set("xmlns", self.NAMESPACE)

        # Add metadata section
        self._add_metadata(root)

        # Add reaction network
        self._add_reaction_network(root)

        # Add constraints
        self._add_constraints(root)

        # Add experiments/configurations
        self._add_experiments(root)

        return root

    def _add_metadata(self, root: ET.Element) -> None:
        """Add metadata element with metadata."""
        if not self._model.metadata:
            return

        metadata = ET.SubElement(root, "info")

        if self._model.metadata.date:
            date_elem = ET.SubElement(metadata, "date")
            date_elem.text = self._model.metadata.date.strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        else:
            # Add current timestamp
            date_elem = ET.SubElement(metadata, "date")
            date_elem.text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if self._model.metadata.comment:
            comment_elem = ET.SubElement(metadata, "comment")
            comment_elem.text = self._model.metadata.comment

        if self._model.metadata.modeler:
            modeler_elem = ET.SubElement(metadata, "modeler")
            modeler_elem.text = self._model.metadata.modeler

        if self._model.metadata.strain:
            strain_elem = ET.SubElement(metadata, "strain")
            strain_elem.text = self._model.metadata.strain

    def _add_reaction_network(self, root: ET.Element) -> None:
        """Add reactionnetwork element."""
        rn = ET.SubElement(root, "reactionnetwork")

        # Compute atom counts from atom transitions for metabolites
        # missing this info
        self._metabolite_atom_counts = self._compute_metabolite_atom_counts()

        # Add metabolite pools
        pools = ET.SubElement(rn, "metabolitepools")
        for metabolite in self._model.model.metabolites:
            pool = ET.SubElement(pools, "pool")
            pool.set("id", metabolite.id)

            # Use metabolite.atoms if available, otherwise infer from
            # atom transitions
            atom_count = metabolite.atoms
            if not atom_count:
                atom_count = self._metabolite_atom_counts.get(metabolite.id)
            if atom_count:
                pool.set("atoms", str(atom_count))

            # Add annotations
            if metabolite.annotations:
                for ann in metabolite.annotations:
                    ann_elem = ET.SubElement(pool, "annotation")
                    ann_elem.set("name", ann.name)
                    if ann.content:
                        ann_elem.text = ann.content

        # Add reactions
        for reaction in self._model.model.reactions:
            self._add_reaction(rn, reaction)

        # Add drain reactions for sink metabolites (produced but never consumed).
        # 13CFlux2 / x3cflux requires every output metabolite to have an explicit
        # sink reaction so the network can be validated as structurally consistent.
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
        for metab_id in sorted(
            m
            for m in (produced - consumed - input_pools)
            if not m.endswith("_ext")
        ):
            drain = ET.SubElement(rn, "reaction")
            drain.set("id", f"{metab_id}_out")
            drain.set("bidirectional", "false")
            reduct = ET.SubElement(drain, "reduct")
            reduct.set("id", metab_id)
            atom_count = self._metabolite_atom_counts.get(metab_id)
            if not atom_count:
                m = next(
                    (
                        m
                        for m in self._model.model.metabolites
                        if m.id == metab_id
                    ),
                    None,
                )
                if m:
                    atom_count = m.atoms
            if atom_count:
                cfg = " ".join(f"C#{i}@1" for i in range(1, atom_count + 1))
                reduct.set("cfg", cfg)

    def _add_reaction(self, rn: ET.Element, reaction) -> None:
        """Add a single reaction element."""
        rxn = ET.SubElement(rn, "reaction")

        # Handle variant reactions
        if (
            reaction.atom_transition_ids
            and len(reaction.atom_transition_ids) > 1
        ):
            # Multiple variants - use space-separated IDs
            rxn.set("id", " ".join(reaction.atom_transition_ids))
        else:
            rxn.set("id", reaction.id)

        # Reversibility (FluxML default is bidirectional=true)
        if not reaction.reversibility:
            rxn.set("bidirectional", "false")

        # Add annotations
        if hasattr(reaction, "annotations") and reaction.annotations:
            for ann in reaction.annotations:
                ann_elem = ET.SubElement(rxn, "annotation")
                ann_elem.set("name", ann.name)
                if ann.content:
                    ann_elem.text = ann.content

        # Get atom transition if exists
        atom_mapping = self._model.model.atom_mappings.get(reaction.id)

        # Add reactants
        for i, reactant_id in enumerate(reaction.reactants):
            reduct = ET.SubElement(rxn, "reduct")
            reduct.set("id", reactant_id)

            # Add cfg if atom transition exists
            if atom_mapping:
                cfg = self._build_reactant_cfg(atom_mapping, reactant_id, i)
                if cfg:
                    reduct.set("cfg", cfg)

        # Add products
        if reaction.products:
            # Check if this reaction has variants (multiple atom maps)
            has_variants = atom_mapping and len(atom_mapping.maps) > 1

            for i, product_id in enumerate(reaction.products):
                rproduct = ET.SubElement(rxn, "rproduct")
                rproduct.set("id", product_id)

                if atom_mapping:
                    if has_variants:
                        # Add variant elements for each mapping
                        self._add_product_variants(
                            rproduct, atom_mapping, product_id, i
                        )
                    else:
                        # Single mapping - add cfg attribute
                        cfg = self._build_product_cfg(
                            atom_mapping, product_id, i
                        )
                        if cfg:
                            rproduct.set("cfg", cfg)

    def _build_reactant_cfg(
        self,
        atom_mapping: AtomTransition,
        cpd_id: str,
        position: int,
    ) -> Optional[str]:
        """Build FluxML cfg string for a reactant."""
        if not atom_mapping.maps:
            return None

        # Get the first atom map (for reactants, all variants should
        # have same reactant atoms)
        atom_map = next(iter(atom_mapping.maps.values()))

        # Find all atoms that come from this reactant
        instance = 1  # Track instance (for repeated compounds)

        # Count how many times this compound appears before this position
        compound_count = 0
        for i, r_id in enumerate(atom_mapping.reactants):
            if r_id == cpd_id:
                compound_count += 1
                if i == position:
                    instance = compound_count
                    break

        # Collect source atoms for this reactant instance
        max_atom = 0
        for src_addr in atom_map.mapping.values():
            if src_addr.mol == cpd_id and src_addr.instance == instance:
                max_atom = max(max_atom, src_addr.index)

        if max_atom == 0:
            return None

        # Build cfg string like "C#1@1 C#2@1 C#3@1"
        cfg_parts = []
        for atom_idx in range(1, max_atom + 1):
            # Find this atom in the mapping values to get element type
            element = "C"  # Default to carbon
            for prod_addr, src_addr in atom_map.mapping.items():
                if (
                    src_addr.mol == cpd_id
                    and src_addr.instance == instance
                    and src_addr.index == atom_idx
                ):
                    element = src_addr.element or "C"
                    break
            cfg_parts.append(f"{element}#{atom_idx}@{position + 1}")

        return " ".join(cfg_parts) if cfg_parts else None

    def _build_product_cfg(
        self,
        atom_mapping: AtomTransition,
        cpd_id: str,
        position: int,
        atom_map_id: Optional[str] = None,
    ) -> Optional[str]:
        """Build FluxML cfg string for a product."""
        if not atom_mapping.maps:
            return None

        # Select which atom map to use
        if atom_map_id:
            atom_map = atom_mapping.maps.get(atom_map_id)
        else:
            atom_map = next(iter(atom_mapping.maps.values()))

        if not atom_map:
            return None

        # Count instance for this product
        instance = 1
        compound_count = 0
        for i, p_id in enumerate(atom_mapping.products):
            if p_id == cpd_id:
                compound_count += 1
                if i == position:
                    instance = compound_count
                    break

        # Build position lookup for reactants
        reactant_position = {}
        instance_counts: Dict[str, int] = {}
        for pos, r_id in enumerate(atom_mapping.reactants, start=1):
            instance_counts[r_id] = instance_counts.get(r_id, 0) + 1
            inst = instance_counts[r_id]
            reactant_position[(r_id, inst)] = pos

        # Collect product atoms for this compound instance
        atoms = []
        for prod_addr, src_addr in atom_map.mapping.items():
            if prod_addr.mol == cpd_id and prod_addr.instance == instance:
                src_pos = reactant_position.get(
                    (src_addr.mol, src_addr.instance)
                )
                if src_pos is not None:
                    element = prod_addr.element or "C"
                    atoms.append(
                        (prod_addr.index, element, src_addr.index, src_pos)
                    )

        if not atoms:
            return None

        # Sort by product atom index and build cfg string
        atoms.sort(key=lambda x: x[0])
        cfg_parts = [
            f"{elem}#{src_idx}@{src_pos}" for _, elem, src_idx, src_pos in atoms
        ]

        return " ".join(cfg_parts)

    def _add_product_variants(
        self,
        rproduct: ET.Element,
        atom_mapping: AtomTransition,
        cpd_id: str,
        position: int,
    ) -> None:
        """Add variant elements for a product with multiple mappings."""
        for map_id, atom_map in atom_mapping.maps.items():
            variant = ET.SubElement(rproduct, "variant")

            # Build cfg for this variant
            cfg = self._build_product_cfg(
                atom_mapping, cpd_id, position, atom_map_id=map_id
            )
            if cfg:
                variant.set("cfg", cfg)

    def _add_constraints(self, root: ET.Element) -> None:
        """Add constraints element."""
        constraints = self._model.constraints
        if not constraints:
            return

        # Check if there are any constraints to write
        has_net = constraints.net and constraints.net.formulas
        has_xch = constraints.xch and constraints.xch.formulas
        has_psize = (
            constraints.metabolitesize and constraints.metabolitesize.formulas
        )

        if not (has_net or has_xch or has_psize):
            return

        constr_elem = ET.SubElement(root, "constraints")

        # Add net constraints
        if has_net:
            net_elem = ET.SubElement(constr_elem, "net")
            textual = ET.SubElement(net_elem, "textual")
            formulas = []
            for formula in constraints.net.formulas:
                formulas.append(formula.expression)
            textual.text = ";\n".join(formulas)

        # Add xch constraints
        if has_xch:
            xch_elem = ET.SubElement(constr_elem, "xch")
            textual = ET.SubElement(xch_elem, "textual")
            formulas = []
            for formula in constraints.xch.formulas:
                formulas.append(formula.expression)
            textual.text = ";\n".join(formulas)

        # Add psize constraints
        if has_psize:
            psize_elem = ET.SubElement(constr_elem, "psize")
            textual = ET.SubElement(psize_elem, "textual")
            formulas = []
            for formula in constraints.metabolitesize.formulas:
                formulas.append(formula.expression)
            textual.text = ";\n".join(formulas)

    def _add_experiments(self, root: ET.Element) -> None:
        """Add configuration elements for experiments."""
        for experiment in self._model.experiments:
            self._add_experiment(root, experiment)

    def _add_experiment(
        self, root: ET.Element, experiment: LabelingExperiments
    ) -> None:
        """Add a single configuration/experiment element."""
        config = ET.SubElement(root, "configuration")
        config.set("name", experiment.name)
        config.set("stationary", "true" if experiment.stationary else "false")

        if experiment.time is not None:
            config.set("time", str(experiment.time))

        # Add tracer inputs
        for tracer in experiment.tracers:
            self._add_tracer_input(config, tracer)

        # Add measurement section
        if experiment.measurement:
            self._add_measurement(config, experiment)

        # Add simulation section
        if experiment.simulation:
            self._add_simulation(config, experiment)

    def _add_simulation(
        self, config: ET.Element, experiment: LabelingExperiments
    ) -> None:
        """Add simulation element with flux and pool-size initial values.

        Writes stored simulation values (F-type fluxes) then solves S·v=0
        to derive values for any reactions not yet listed — this covers:
        - D-type reactions from MTF parsing (E_out, F_out with no stored value)
        - Auto-generated drain reactions for sink metabolites (freeflux source)
        """
        import numpy as np

        sim = ET.SubElement(config, "simulation")
        sim.set("method", "auto")
        sim.set("type", "full")

        variables = ET.SubElement(sim, "variables")
        stored_ids: set = set()
        for fv in experiment.simulation.variables.flux_values:
            # Skip D-type fluxes with no stored value (influx_si MTF parsing
            # produces these; they must be derived from mass balance instead).
            if fv.value is None:
                continue
            fv_elem = ET.SubElement(variables, "fluxvalue")
            fv_elem.set("flux", fv.flux)
            fv_elem.set("type", fv.type or "net")
            fv_elem.text = str(fv.value)
            stored_ids.add(fv.flux)

        # Build the full reaction list: model reactions + auto-generated drains.
        # "Auto-generated drains" are for sink metabolites not already consumed
        # by any existing model reaction (and not _ext external boundaries).
        reactions = list(self._model.model.reactions)
        consumed: set = set()
        produced: set = set()
        for rxn in reactions:
            consumed.update(rxn.reactants)
            produced.update(rxn.products)
        input_pools = {
            t.metabolite for exp in self._model.experiments for t in exp.tracers
        }
        sink_mets = sorted(
            m
            for m in (produced - consumed - input_pools)
            if not m.endswith("_ext")
        )
        new_drain_ids = [f"{m}_out" for m in sink_mets]

        # Collect ALL reaction IDs (model + new drains) that need a value written.
        rxn_names = [r.id for r in reactions]
        all_names = rxn_names + new_drain_ids  # new drains appended at end
        missing = [n for n in all_names if n not in stored_ids]

        if missing:
            metabolites = list(self._model.model.metabolites)
            met_idx = {m.id: i for i, m in enumerate(metabolites)}
            n_rxn = len(rxn_names)

            S = np.zeros((len(metabolites), len(all_names)))
            for j, rxn in enumerate(reactions):
                for mid in rxn.reactants:
                    if mid in met_idx:
                        S[met_idx[mid], j] -= 1.0
                for mid in rxn.products:
                    if mid in met_idx:
                        S[met_idx[mid], j] += 1.0
            for k, mid in enumerate(sink_mets):
                if mid in met_idx:
                    S[met_idx[mid], n_rxn + k] -= 1.0

            int_rows = [
                i
                for i, m in enumerate(metabolites)
                if m.id not in input_pools and not m.id.endswith("_ext")
            ]
            S_int = S[int_rows, :]

            stored_vals = {
                fv.flux: (fv.value or 0.0)
                for fv in experiment.simulation.variables.flux_values
            }
            f_names = [n for n in all_names if n in stored_ids]
            d_names = [n for n in all_names if n not in stored_ids]
            f_cols = [i for i, n in enumerate(all_names) if n in stored_ids]
            d_cols = [i for i, n in enumerate(all_names) if n not in stored_ids]

            derived: dict = {}
            if d_cols:
                v_f = np.array([stored_vals.get(n, 0.0) for n in f_names])
                try:
                    v_d, *_ = np.linalg.lstsq(
                        S_int[:, d_cols], -S_int[:, f_cols] @ v_f, rcond=None
                    )
                    derived = {
                        d_names[i]: float(v_d[i]) for i in range(len(d_cols))
                    }
                except Exception:
                    derived = {n: 0.0 for n in d_names}

            for rxn_id in missing:
                fv_elem = ET.SubElement(variables, "fluxvalue")
                fv_elem.set("flux", rxn_id)
                fv_elem.set("type", "net")
                fv_elem.text = str(derived.get(rxn_id, 0.0))

        for pv in experiment.simulation.variables.metabolitesize_values:
            pv_elem = ET.SubElement(variables, "poolvalue")
            pv_elem.set("pool", pv.metabolite)
            pv_elem.text = str(pv.value)

    def _add_tracer_input(self, config: ET.Element, tracer) -> None:
        """Add input element for tracer specification."""
        input_elem = ET.SubElement(config, "input")
        input_elem.set("pool", tracer.metabolite)
        input_elem.set("type", "isotopomer")

        for label in tracer.labels:
            label_elem = ET.SubElement(input_elem, "label")
            label_elem.set("cfg", label.labeled_pattern)
            if hasattr(label, "purity") and label.purity:
                label_elem.set("purity", str(label.purity))
            label_elem.text = str(label.fraction)

    def _add_measurement(
        self, config: ET.Element, experiment: LabelingExperiments
    ) -> None:
        """Add measurement section to configuration."""
        meas_elem = ET.SubElement(config, "measurement")

        # Add mlabel with metadata
        mlabel = ET.SubElement(meas_elem, "mlabel")
        date_elem = ET.SubElement(mlabel, "date")
        date_elem.text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Check if this is time-course data
        if experiment.measurement.data and experiment.measurement.data.data:
            times = set()
            for datum in experiment.measurement.data.data:
                if datum.time is not None:
                    times.add(datum.time)
            if times:
                timeunit = ET.SubElement(mlabel, "timeunit")
                timeunit.text = "s"  # Default to seconds

        # Add model section
        model_elem = ET.SubElement(meas_elem, "model")

        # Add labeling measurements
        if experiment.measurement.model.labeling_measurement:
            labeling = ET.SubElement(model_elem, "labelingmeasurement")
            self._add_labeling_groups(
                labeling,
                experiment.measurement.model.labeling_measurement,
                experiment.measurement.data,
            )

        # Add flux measurements
        if experiment.measurement.model.flux_measurement:
            flux_meas = ET.SubElement(model_elem, "fluxmeasurement")
            self._add_flux_measurements(
                flux_meas, experiment.measurement.model.flux_measurement
            )

        # Add data section
        if experiment.measurement.data and experiment.measurement.data.data:
            data_elem = ET.SubElement(meas_elem, "data")

            # Add dlabel
            dlabel = ET.SubElement(data_elem, "dlabel")
            date_elem = ET.SubElement(dlabel, "date")
            date_elem.text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Add data points
            for datum in experiment.measurement.data.data:
                datum_elem = ET.SubElement(data_elem, "datum")
                datum_elem.set("id", datum.id.split(":")[0])  # Group ID
                datum_elem.set("stddev", str(datum.stddev))
                if datum.time is not None:
                    datum_elem.set("time", str(datum.time))
                # Weight is typically the isotopomer index
                if datum.weight is not None:
                    datum_elem.set("weight", str(datum.weight))
                elif "M" in datum.id:
                    weight = datum.id.split("M")[-1]
                    try:
                        datum_elem.set("weight", weight)
                    except (ValueError, TypeError):
                        pass
                datum_elem.text = str(datum.value)

            # Apply ¹³C natural abundance correction if requested
            if getattr(self, "_apply_na_correction", False):
                self._apply_na_to_data(data_elem)

    def _na_correction_matrix(self, n_C: int) -> "np.ndarray":
        """Binomial ¹³C natural-abundance correction matrix M_C of size (n_C+1)×(n_C+1).

        M_C[i,j] = Binom(n_C-j, i-j, p) for i>=j, else 0.
        Transforms a pure-tracer (no NA) MID vector into the expected
        measured (with NA) MID: mid_measured = M_C @ mid_tracer.
        p = NATURAL_ABUNDANCE_CARBON = 0.01109 (x3cflux hardcoded constant).
        """
        import numpy as np
        from scipy.stats import binom as sp_binom

        p = self.NATURAL_ABUNDANCE_CARBON
        n = n_C
        M = np.zeros((n + 1, n + 1))
        for j in range(n + 1):
            for extra in range(n - j + 1):
                i = j + extra
                M[i, j] = sp_binom.pmf(extra, n - j, p)
        return M

    def _apply_na_to_data(self, data_elem: ET.Element) -> None:
        """Re-write <datum> children of data_elem with ¹³C NA correction applied.

        Groups datums by group id, applies M_C to each group's (value, stddev)
        vectors, and rewrites the text / stddev attributes in place.
        Only groups with integer weights (MS isotopologue data) are corrected.
        """
        import numpy as np

        datums = list(data_elem)
        # Group datums: {group_id: [(weight_int, elem), ...]}
        from collections import defaultdict

        groups: dict = defaultdict(list)
        non_weight: list = []
        for elem in datums:
            if elem.tag == "datum":
                try:
                    w = int(elem.get("weight", ""))
                    groups[elem.get("id", "")].append((w, elem))
                except (ValueError, TypeError):
                    non_weight.append(elem)

        for gid, pairs in groups.items():
            pairs_sorted = sorted(pairs, key=lambda x: x[0])
            weights = [w for w, _ in pairs_sorted]
            if weights != list(range(len(weights))):
                continue  # Non-sequential weights — skip (NMR, etc.)
            n_C = len(weights) - 1
            M = self._na_correction_matrix(n_C)
            vals = np.array([float(e.text or "0") for _, e in pairs_sorted])
            sds = np.array(
                [float(e.get("stddev", "0")) for _, e in pairs_sorted]
            )
            new_vals = M @ vals
            new_sds = np.sqrt((M**2) @ (sds**2))
            for i, (_, elem) in enumerate(pairs_sorted):
                elem.text = f"{new_vals[i]:.8e}"
                elem.set("stddev", f"{new_sds[i]:.8e}")

    def _add_labeling_groups(
        self, labeling: ET.Element, labeling_meas, data
    ) -> None:
        """Add labeling measurement groups."""
        # Collect time points from data
        times = set()
        if data and data.data:
            for datum in data.data:
                if datum.time is not None:
                    times.add(datum.time)

        times_str = ",".join(str(t) for t in sorted(times)) if times else "0.0"

        for group in labeling_meas.groups:
            group_elem = ET.SubElement(labeling, "group")
            group_elem.set("id", group.id)
            group_elem.set("scale", group.scale or "auto")
            if times:
                group_elem.set("times", times_str)

            # Add textual expression
            textual = ET.SubElement(group_elem, "textual")
            if group.expression:
                textual.text = group.expression.textual
            else:
                textual.text = group.id

    def _add_flux_measurements(
        self, flux_meas: ET.Element, flux_measurement
    ) -> None:
        """Add flux measurement elements."""
        for net_flux in flux_measurement.net_fluxes:
            netflux_elem = ET.SubElement(flux_meas, "netflux")
            netflux_elem.set("id", net_flux.id)

            textual = ET.SubElement(netflux_elem, "textual")
            if net_flux.expression:
                textual.text = net_flux.expression.textual
            else:
                textual.text = net_flux.id


def write_fluxml(
    model: FluxomicsData,
    filepath: str,
    apply_na_correction: bool = False,
) -> None:
    """Convenience function to write FluxomicsData to a FluxML file.

    Args:
        model: FluxomicsData to write
        filepath: Path to output file
        apply_na_correction: Apply ¹³C natural abundance correction (p=0.01109/C)
            to all MS measurement data.  Use when data is synthetic (no NA present)
            and the target tool is x3cflux (which models NA internally).
    """
    writer = FluxMLWriter()
    writer.write(model, filepath, apply_na_correction=apply_na_correction)
