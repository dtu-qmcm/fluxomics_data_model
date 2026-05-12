"""Core Fluxomics Data Converter structures.

This module defines the four top-level classes that form the root of the
data hierarchy:

- :class:`Metadata`  — optional (modeler, strain, date …)
- :class:`MetabolicNetworkModel`     — required (metabolites, reactions, atom
  mappings)
- :class:`LabelingExperiments` — required (tracers, ,
  measurements, simulation variables)
- :class:`FluxomicsData` — the root ``<fluxml>`` element; holds one
  ``MetabolicNetworkModel`` and one or more ``LabelingExperiments``

All classes are immutable Pydantic models (``frozen=True``).  Numerical
(JAX) operations that require a stable ordering of metabolites and reactions
must use :attr:`MetabolicNetworkModel.ordered_metabolite_ids` and
:attr:`MetabolicNetworkModel.ordered_reaction_ids`.
"""

from typing import Optional, List, Dict, Any
from datetime import datetime
import re
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from ..model.metabolite import Metabolite
from .common import DictList, AtomTransitionsNetwork
from ..model.reaction import Reaction
from ..model.constraint import Constraints
from ..experiment.measurement import Measurement
from ..experiment.tracer import Tracers
from ..output.simulation import Simulation


class Metadata(BaseModel):
    """Optional metadata block.

    All fields are optional; a ``Metadata`` object with all ``None`` values
    is valid.  The ``date`` field accepts both ``datetime`` objects and
    FluxML-style timestamp strings (``"YYYY-MM-DD HH:MM:SS"``).

    Corresponds to ``fluxml/info``.
    """

    name: Optional[str] = Field(default=None, description="Model name")
    version: Optional[str] = Field(default=None, description="Model version")
    date: Optional[datetime] = Field(
        default=None, description="Creation/modification timestamp"
    )
    comment: Optional[str] = Field(
        default=None, description="Model description"
    )
    signature: Optional[bytes] = Field(
        default=None, description="Digital signature"
    )
    modeler: Optional[str] = Field(
        default=None, description="Modeler information"
    )
    strain: Optional[str] = Field(
        default=None, description="Strain information"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("date", mode="before")
    @classmethod
    def parse_date(cls, v):
        """Parse FluxML timestamp format: YYYY-MM-DD HH:MM:SS."""
        if isinstance(v, str):
            # FluxML timestamp pattern: YYYY-MM-DD HH:MM:SS
            if re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$", v):
                return datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
        return v


class MetabolicNetworkModel(BaseModel):
    """Metabolic network model.

    Holds all metabolites, reactions, and atom transitions.  Cross-references
    (every reactant/product ID must exist in ``metabolites``) are validated
    at construction time.

    Ordering note
    -------------
    ``metabolite_ids`` and ``reaction_ids`` are ``frozenset`` values and
    therefore have **no defined iteration order**.  Use
    :attr:`ordered_metabolite_ids` and :attr:`ordered_reaction_ids` whenever
    you need a deterministic sequence that corresponds to matrix rows/columns.

    """

    metabolites: DictList[Metabolite] = Field(
        default_factory=lambda: DictList[Metabolite](),
        description="Metabolite definitions",
    )
    reactions: DictList[Reaction] = Field(
        default_factory=lambda: DictList[Reaction](),
        description="Reaction definitions",
    )
    atom_mappings: AtomTransitionsNetwork = Field(
        default_factory=AtomTransitionsNetwork,
        description="Atom transitions keyed by reaction ID",
    )
    compartments: List[str] = Field(
        default_factory=list, description="List of compartments in the model"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
        """Initialise and validate metabolite/reaction cross-references."""
        super().__init__(**data)
        self._validate_cross_references()

    def _validate_cross_references(self):
        """Validate cross-references between metabolites and reactions."""
        metabolite_ids = self.metabolite_ids

        # Check all reactant and product IDs exist in metabolites
        for reaction in self.reactions:
            for reactant_id in reaction.reactants:
                if reactant_id not in metabolite_ids:
                    raise ValueError(
                        f"Reactant {reactant_id} in reaction {reaction.id} "
                        f"references unknown metabolite"
                    )
            for product_id in reaction.products:
                if product_id not in metabolite_ids:
                    raise ValueError(
                        f"Product {product_id} in reaction {reaction.id} "
                        f"references unknown metabolite"
                    )

    @property
    def metabolite_ids(self) -> frozenset[str]:
        """Get all metabolite IDs."""
        return frozenset(self.metabolites.ids)

    @property
    def reaction_ids(self) -> frozenset[str]:
        """Get all reaction IDs (biological representation)."""
        return frozenset(self.reactions.ids)

    @property
    def atom_transition_ids(self) -> List[str]:
        """All ids for atom transitions including variants."""
        ids = []
        for reaction in self.reactions:
            if reaction.atom_transition_ids is not None:
                ids.extend(reaction.atom_transition_ids)
        return ids

    @property
    def computational_reaction_ids(self) -> frozenset[str]:
        """Get all computational reaction IDs (both variants and base IDs)."""
        ids = []
        for reaction in self.reactions:
            if reaction.atom_transition_ids is not None:
                ids.append(reaction.id)
                ids.extend(reaction.atom_transition_ids)
            else:
                ids.append(reaction.id)
        return frozenset(ids)

    @property
    def ordered_metabolite_ids(self) -> List[str]:
        """Return metabolite IDs in stable insertion order.

        Use this (not ``metabolite_ids``) whenever the order of IDs must
        correspond to rows/columns of a numerical matrix, so that results
        are reproducible across interpreter runs.
        """
        return list(self.metabolites.ids)

    @property
    def ordered_reaction_ids(self) -> List[str]:
        """Return reaction IDs in stable insertion order.

        Use this (not ``reaction_ids``) whenever the order of IDs must
        correspond to columns of a numerical matrix.
        """
        return list(self.reactions.ids)

    def get_stoichiometric_matrix(self) -> jnp.ndarray:
        """Build a stoichiometric matrix of shape (n_metabolites, n_reactions).

        Rows correspond to metabolites and columns to reactions, both in
        stable insertion order (see :attr:`ordered_metabolite_ids` and
        :attr:`ordered_reaction_ids`).  Using ``DictList.ids`` instead of
        iterating over a ``frozenset`` guarantees that the row ordering is
        deterministic across interpreter runs regardless of Python's hash
        randomisation (``PYTHONHASHSEED``).

        Returns:
            JAX array of shape ``(n_metabolites, n_reactions)`` where
            ``S[i, j]`` is the stoichiometric coefficient of metabolite *i*
            in reaction *j* (negative for reactants, positive for products).
        """
        # Use insertion-ordered list, NOT self.metabolite_ids (a frozenset).
        metabolite_ids = self.ordered_metabolite_ids

        matrix = []
        for reaction in self.reactions:
            column = reaction.get_stoichiometric_vector(metabolite_ids)
            matrix.append(column)

        return jnp.stack(matrix, axis=1)


class LabelingExperiments(BaseModel):
    """Isotope labeling experiment data.

    A specific set of tracers, local constraints,
    measurement data, and simulation variables.  The ``name`` attribute must
    be unique within a :class:`FluxomicsData`; 

    """

    name: str = Field(description="Experiment name")
    stationary: bool = Field(default=True, description="Stationary assumption")
    time: Optional[float] = Field(default=None, description="Time point")
    comment: Optional[str] = Field(
        default=None, description="Experiment comment"
    )
    tracers: List[Tracers] = Field(
        default_factory=list, description="Tracer specifications"
    )
    constraints: Optional[Constraints] = Field(
        default=None, description="Experiment-specific constraints"
    )
    measurement: Optional[Measurement] = Field(
        default=None, description="Measurement data"
    )
    simulation: Optional[Simulation] = Field(
        default=None, description="Experiment-specific simulation settings"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def traced_metabolites(self) -> frozenset[str]:
        """Get all traced metabolite IDs."""
        return frozenset(tracer_spec.metabolite for tracer_spec in self.tracers)

    def get_tracer_for_metabolite(
        self, metabolite_id: str
    ) -> Optional[Tracers]:
        """Get tracer specification for a metabolite."""
        for tracer_spec in self.tracers:
            if tracer_spec.metabolite == metabolite_id:
                return tracer_spec
        return None

    def get_tracer_composition_matrix(
        self, metabolite_ids: List[str]
    ) -> jnp.ndarray:
        """Get tracer composition matrix for JAX computations.

        Returns:
            JAX array of shape (n_metabolites, n_isotopomers)
        """
        compositions = []
        for metabolite_id in metabolite_ids:
            tracer_spec = self.get_tracer_for_metabolite(metabolite_id)
            if tracer_spec:
                compositions.append(tracer_spec.composition_vector)
            else:
                compositions.append(jnp.array([1.0]))  # Unlabeled

        # Pad to same length
        max_len = max(len(comp) for comp in compositions)
        padded = [
            jnp.pad(comp, (0, max_len - len(comp))) for comp in compositions
        ]

        return jnp.stack(padded)


class FluxomicsData(BaseModel):
    """Root data model for fluxomics data.

    This is the primary entry point for all fluxomics data.  It combines one
    :class:`MetabolicNetworkModel` (the reaction network) with one or more
    :class:`LabelingExperiments` (experimental configurations), optional global
    :class:`~fluxomics_data_converter.model.constraint.Constraints`, and
    optional :class:`Metadata`.

    Cross-references between experiments and the model (tracer metabolites,
    simulation flux variables) are validated at construction time.

    Numerical usage
    ---------------
    Call :meth:`to_jax_representation` to obtain a dictionary of JAX arrays
    suitable for gradient-based optimisation.  When the model contains more
    than one experiment, pass ``experiment_name`` to select which
    experiment's bounds to use::

        jax_data = fdm.to_jax_representation(experiment_name="default")
        S = jax_data["stoichiometric_matrix"]   # shape (n_met, n_rxn)

    """

    model: MetabolicNetworkModel = Field(
        description="Metabolic network model definition"
    )
    metadata: Optional[Metadata] = Field(default=None, description="Metadata")
    constraints: Optional[Constraints] = Field(
        default=None, description="Model constraints"
    )
    experiments: List[LabelingExperiments] = Field(
        default_factory=list, description="Experimental setups"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
        """Initialise and validate experiment cross-references."""
        super().__init__(**data)
        self._validate_experiments_references()

    def _validate_experiments_references(self):
        """Validate experiment references to metabolites and reactions."""
        metabolite_ids = self.metabolite_ids
        # Use computational reaction IDs for validation (includes variant IDs)
        computational_reaction_ids = self.model.computational_reaction_ids

        # Also allow drain reactions auto-generated for sink metabolites
        # (produced but never consumed and not an input pool). These are not
        # stored in the model but are created during FML/MTF writing.
        all_produced: set = set()
        all_consumed: set = set()
        for rxn in self.model.reactions:
            all_produced.update(rxn.products)
            all_consumed.update(rxn.reactants)
        input_pools = frozenset(
            t.metabolite
            for exp in self.experiments
            for t in exp.tracers
        )
        drain_ids = frozenset(
            f"{m}_out" for m in (all_produced - all_consumed - input_pools)
        )
        valid_reaction_ids = computational_reaction_ids | drain_ids

        # Check experiment names are unique
        if len(self.experiments) > 1:
            experiments_names = [e.name for e in self.experiments]
            if len(set(experiments_names)) != len(experiments_names):
                raise ValueError("Experiment names must be unique")

        # Check metabolite and reaction references
        for experiment in self.experiments:
            # Check tracer metabolite references
            for tracer_spec in experiment.tracers:
                if tracer_spec.metabolite not in metabolite_ids:
                    raise ValueError(
                        f"Tracer metabolite {tracer_spec.metabolite} in "
                        f"experiment {experiment.name} references unknown "
                        f"metabolite"
                    )

            # Check simulation variable references
            if experiment.simulation and experiment.simulation.variables:
                for flux_var in experiment.simulation.variables.flux_values:
                    # Check against computational IDs (includes variants) plus
                    # auto-generated drain reactions for sink metabolites
                    if flux_var.flux not in valid_reaction_ids:
                        raise ValueError(
                            f"Flux variable {flux_var.flux} in configuration "
                            f"{experiment.name} references unknown reaction"
                        )

                for (
                    m_var
                ) in experiment.simulation.variables.metabolitesize_values:
                    if m_var.metabolite not in metabolite_ids:
                        raise ValueError(
                            f"Metabolite size variable {m_var.metabolite} in "
                            f"experiment {experiment.name} references "
                            f"unknown metabolite"
                        )

    @property
    def metabolite_ids(self) -> frozenset[str]:
        """Get all metabolite IDs in the model."""
        return self.model.metabolite_ids

    @property
    def reaction_ids(self) -> frozenset[str]:
        """Get all reaction IDs in the model."""
        return self.model.reaction_ids

    @property
    def experiments_names(self) -> frozenset[str]:
        """Get all experiment names."""
        return frozenset(exp.name for exp in self.experiments)

    def get_experiments(self, name: str) -> Optional[LabelingExperiments]:
        """Get experiments by name."""
        for exp in self.experiments:
            if exp.name == name:
                return exp
        return None

    def to_jax_representation(
        self, experiment_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """Convert to a JAX-compatible dictionary for numerical computation.

        Metabolite and reaction IDs are returned in stable insertion order
        (matching the rows/columns of the stoichiometric matrix).  When the
        model contains more than one experiment, ``experiment_name`` must be
        provided so that the correct flux and pool-size bounds are used.

        Args:
            experiment_name: Name of the experiment whose simulation bounds
                should be included.  Required when the model contains more
                than one experiment.  If the model has exactly one experiment
                and this argument is ``None``, that experiment is used
                automatically.  If the model has no experiments, bounds
                default to ``[-inf, inf]`` for fluxes and ``[0, inf]`` for
                pool sizes.

        Returns:
            Dictionary with the following keys:

            - ``"stoichiometric_matrix"``: JAX array of shape
              ``(n_metabolites, n_reactions)``
            - ``"flux_bounds"``: JAX array of shape ``(n_reactions, 2)``
            - ``"metabolitesize_bounds"``: JAX array of shape
              ``(n_metabolites, 2)``
            - ``"metabolite_ids"``: ordered list of metabolite ID strings
            - ``"reaction_ids"``: ordered list of reaction ID strings
            - ``"n_metabolites"``: int
            - ``"n_reactions"``: int
            - ``"n_experiments"``: int

        Raises:
            ValueError: If the model has more than one experiment and
                ``experiment_name`` is not provided, or if
                ``experiment_name`` does not match any experiment.
        """
        # Use insertion-ordered lists for deterministic matrix indexing.
        metabolite_ids = self.model.ordered_metabolite_ids
        reaction_ids = self.model.ordered_reaction_ids

        # Get stoichiometric matrix
        S = self.model.get_stoichiometric_matrix()

        # Resolve which experiment supplies the simulation bounds.
        experiment = None
        if experiment_name is not None:
            experiment = self.get_experiments(experiment_name)
            if experiment is None:
                raise ValueError(
                    f"Experiment '{experiment_name}' not found. "
                    f"Available: {sorted(self.experiments_names)}"
                )
        elif len(self.experiments) == 1:
            experiment = self.experiments[0]
        elif len(self.experiments) > 1:
            raise ValueError(
                f"Model has {len(self.experiments)} experiments. "
                f"Provide experiment_name to select one. "
                f"Available: {sorted(self.experiments_names)}"
            )

        flux_bounds = None
        metabolitesize_bounds = None

        if experiment and experiment.simulation:
            flux_bounds, metabolitesize_bounds = (
                experiment.simulation.get_optimization_bounds(
                    reaction_ids, metabolite_ids
                )
            )

        if flux_bounds is None:
            flux_bounds = jnp.array([[-jnp.inf, jnp.inf]] * len(reaction_ids))
        if metabolitesize_bounds is None:
            metabolitesize_bounds = jnp.array(
                [[0.0, jnp.inf]] * len(metabolite_ids)
            )

        return {
            "stoichiometric_matrix": S,
            "flux_bounds": flux_bounds,
            "metabolitesize_bounds": metabolitesize_bounds,
            "metabolite_ids": metabolite_ids,
            "reaction_ids": reaction_ids,
            "n_metabolites": len(metabolite_ids),
            "n_reactions": len(reaction_ids),
            "n_experiments": len(self.experiments),
        }

    def get_tracer_experiment_data(
        self, experiments_name: str
    ) -> Optional[Dict[str, Any]]:
        """Get tracer experiment data for a specific experiment.

        Returns:
            Dictionary with JAX arrays for tracer experiment analysis
        """
        experiment = self.get_experiments(experiments_name)
        if not experiment:
            return None

        # Use insertion-ordered list for deterministic matrix indexing.
        metabolite_ids = self.model.ordered_metabolite_ids

        # Get tracer composition matrix
        tracer_matrix = experiment.get_tracer_composition_matrix(metabolite_ids)

        # Get measurement data
        measurement_data = None
        if experiment.measurement:
            values = experiment.measurement.data.values
            errors = experiment.measurement.data.errors
            times = experiment.measurement.data.times

            measurement_data = {
                "values": values,
                "errors": errors,
                "times": times,
            }

        return {
            "tracer_composition": tracer_matrix,
            "measurement_data": measurement_data,
            "stationary": experiment.stationary,
            "time_point": experiment.time,
        }

    def __repr__(self) -> str:
        """Return a summary of the data model including Metadata and Stats."""
        lines = ["Fluxomics Data Converter Summary", "=" * 30, ""]

        # Metadata table
        if self.metadata:
            lines.append("Model Information:")
            lines.append("-" * 18)
            info_items = [
                ("Name", self.metadata.name),
                ("Version", self.metadata.version),
                (
                    "Date",
                    self.metadata.date.strftime("%Y-%m-%d %H:%M:%S")
                    if self.metadata.date
                    else None,
                ),
                ("Comment", self.metadata.comment),
                ("Modeler", self.metadata.modeler),
                ("Strain", self.metadata.strain),
            ]

            max_key_len = max(len(key) for key, _ in info_items)
            for key, value in info_items:
                if value is not None:
                    lines.append(f"{key:<{max_key_len}} : {value}")
        else:
            lines.append("Model Information: Not available")

        lines.append("")

        # Counts
        lines.append("Model Components:")
        lines.append("-" * 17)
        lines.append(f"Reactions        : {len(self.model.reactions)}")  # noqa: E501
        lines.append(f"Metabolites      : {len(self.model.metabolites)}")
        lines.append(f"Experiments      : {len(self.experiments)}")

        return "\n".join(lines)
