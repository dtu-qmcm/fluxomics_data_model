"""
Core FluxML data structures.
"""

from typing import Optional, List, Dict, Any
from datetime import datetime
import re
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from ..model.metabolite import Metabolite
from .common import DictList, AtomMappingsDict
from ..model.reaction import Reaction
from ..model.constraint import Constraints
from ..experiment.measurement import Measurement
from ..experiment.tracer import Tracers
from ..output.simulation import Simulation


class Metadata(BaseModel):
    """
    FluxML info section containing metadata.

    Corresponds to fluxml/info
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
        """Parse FluxML timestamp format: YYYY-MM-DD HH:MM:SS"""
        if isinstance(v, str):
            # FluxML timestamp pattern: YYYY-MM-DD HH:MM:SS
            if re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$", v):
                return datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
        return v


class Model(BaseModel):
    """
    FluxML model containing metabolites and reactions.

    Corresponds to fluxml/reactionnetwork
    """

    metabolites: DictList[Metabolite] = Field(
        default_factory=lambda: DictList[Metabolite](),
        description="Metabolite definitions",
    )
    reactions: DictList[Reaction] = Field(
        default_factory=lambda: DictList[Reaction](),
        description="Reaction definitions",
    )
    atom_mappings: AtomMappingsDict = Field(
        default_factory=AtomMappingsDict,
        description="Atom mappings keyed by reaction ID",
    )
    compartments: List[str] = Field(
        default_factory=list, description="List of compartments in the model"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
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
    def atom_mapping_ids(self) -> List[str]:
        """All ids for atom mappings including variants."""
        ids = []
        for reaction in self.reactions:
            if reaction.atom_mapping_ids is not None:
                ids.extend(reaction.atom_mapping_ids)
        return ids

    @property
    def computational_reaction_ids(self) -> frozenset[str]:
        """Get all computational reaction IDs (includes variants and base IDs)."""
        ids = []
        for reaction in self.reactions:
            if reaction.atom_mapping_ids is not None:
                ids.append(reaction.id)
                ids.extend(reaction.atom_mapping_ids)
            else:
                ids.append(reaction.id)
        return frozenset(ids)

    def get_stoichiometric_matrix(self) -> jnp.ndarray:
        """
        Get stoichiometric matrix for JAX computations.

        Returns:
            JAX array of shape (n_metabolites, n_reactions)
        """
        metabolite_ids = list(self.metabolite_ids)
        # reaction_ids = list(self.reaction_ids)  # Unused variable

        matrix = []
        for reaction in self.reactions:
            column = reaction.get_stoichiometric_vector(metabolite_ids)
            matrix.append(column)

        return jnp.stack(matrix, axis=1)


class Experiments(BaseModel):
    """
    FluxML experimental setup.

    Corresponds to fluxml/experiments
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
        """
        Get tracer composition matrix for JAX computations.

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


class FluxomicsDataModel(BaseModel):
    """
    Root FluxML object containing complete model specification.

    This is the main entry point for FluxML data, designed for JAX compatibility
    with immutable data structures and validation.

    Corresponds to fluxml root element
    """

    model: Model = Field(description="Model definition")
    info: Optional[Metadata] = Field(default=None, description="Model metadata")
    constraints: Optional[Constraints] = Field(
        default=None, description="Model constraints"
    )
    experiments: List[Experiments] = Field(
        default_factory=list, description="Experimental setups"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
        super().__init__(**data)
        self._validate_experiments_references()

    def _validate_experiments_references(self):
        """Validate experiment references to metabolites and reactions."""
        metabolite_ids = self.metabolite_ids
        # Use computational reaction IDs for validation (includes variant IDs)
        computational_reaction_ids = self.model.computational_reaction_ids

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
                    # Check against computational IDs (includes variants)
                    if flux_var.flux not in computational_reaction_ids:
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

    def get_experiments(self, name: str) -> Optional[Experiments]:
        """Get experiments by name."""
        for exp in self.experiments:
            if exp.name == name:
                return exp
        return None

    def to_jax_representation(self) -> Dict[str, Any]:
        """
        Convert to JAX-compatible dictionary representation.

        Returns:
            Dictionary with JAX arrays for numerical computations
        """
        metabolite_ids = list(self.metabolite_ids)
        reaction_ids = list(self.reaction_ids)

        # Get stoichiometric matrix
        S = self.model.get_stoichiometric_matrix()

        # Get bounds matrices from first experiment (if available)
        flux_bounds = None
        metabolitesize_bounds = None

        if self.experiments:
            experiment = self.experiments[0]
            if experiment.simulation:
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
        """
        Get tracer experiment data for a specific experiment.

        Returns:
            Dictionary with JAX arrays for tracer experiment analysis
        """
        experiment = self.get_experiments(experiments_name)
        if not experiment:
            return None

        metabolite_ids = list(self.metabolite_ids)

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
        """
        Return a summary of the data model including Metadata and Stats.
        """
        lines = ["Fluxomics Data Model Summary", "=" * 30, ""]

        # Metadata table
        if self.info:
            lines.append("Model Information:")
            lines.append("-" * 18)
            info_items = [
                ("Name", self.info.name),
                ("Version", self.info.version),
                (
                    "Date",
                    self.info.date.strftime("%Y-%m-%d %H:%M:%S")
                    if self.info.date
                    else None,
                ),
                ("Comment", self.info.comment),
                ("Modeler", self.info.modeler),
                ("Strain", self.info.strain),
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
