"""
Core FluxML data structures.
"""

from typing import Optional, List, Dict, Any
from datetime import datetime
import re
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .pools import MetabolitePools
from .reactions import Reaction
from .constraints import Constraints
from .configuration import Configuration


class Info(BaseModel):
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


class ReactionNetwork(BaseModel):
    """
    FluxML reaction network containing metabolite pools and reactions.

    Corresponds to fluxml/reactionnetwork
    """

    metabolitepools: MetabolitePools = Field(
        description="Metabolite pool definitions"
    )
    reactions: List[Reaction] = Field(
        default_factory=list, description="Reaction definitions"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
        super().__init__(**data)
        self._validate_cross_references()

    def _validate_cross_references(self):
        """Validate cross-references between pools and reactions."""
        pool_ids = self.metabolitepools.pool_ids

        # Check all reduct and rproduct IDs exist in pools
        for reaction in self.reactions:
            for reduct in reaction.reducts:
                if reduct.id not in pool_ids:
                    raise ValueError(
                        f"Reduct {reduct.id} in reaction {reaction.id} "
                        f"references unknown pool"
                    )
            for rproduct in reaction.rproducts:
                if rproduct.id not in pool_ids:
                    raise ValueError(
                        f"RProduct {rproduct.id} in reaction {reaction.id} "
                        f"references unknown pool"
                    )

    @field_validator("reactions")
    @classmethod
    def validate_unique_reaction_ids(cls, v: List[Reaction]) -> List[Reaction]:
        """Validate reaction IDs are unique."""
        ids = [reaction.id for reaction in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Reaction IDs must be unique")
        return v

    @property
    def reaction_dict(self) -> Dict[str, Reaction]:
        """Get reactions as dictionary keyed by ID."""
        return {reaction.id: reaction for reaction in self.reactions}

    @property
    def reaction_ids(self) -> frozenset[str]:
        """Get all reaction IDs."""
        return frozenset(reaction.id for reaction in self.reactions)

    def get_stoichiometric_matrix(self) -> jnp.ndarray:
        """
        Get stoichiometric matrix for JAX computations.

        Returns:
            JAX array of shape (n_pools, n_reactions)
        """
        pool_ids = list(self.metabolitepools.pool_ids)
        # reaction_ids = list(self.reaction_ids)  # Unused variable

        matrix = []
        for reaction in self.reactions:
            column = reaction.get_stoichiometric_vector(pool_ids)
            matrix.append(column)

        return jnp.stack(matrix, axis=1)


class FluxML(BaseModel):
    """
    Root FluxML object containing complete model specification.

    This is the main entry point for FluxML data, designed for JAX compatibility
    with immutable data structures and validation.

    Corresponds to fluxml root element
    """

    reactionnetwork: ReactionNetwork = Field(
        description="Reaction network definition"
    )
    info: Optional[Info] = Field(default=None, description="Model metadata")
    constraints: Optional[Constraints] = Field(
        default=None, description="Model constraints"
    )
    configurations: List[Configuration] = Field(
        default_factory=list, description="Experimental configurations"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def __init__(self, **data):
        super().__init__(**data)
        self._validate_configuration_references()

    def _validate_configuration_references(self):
        """Validate configuration references to pools and reactions."""
        pool_ids = self.pool_ids
        reaction_ids = self.reaction_ids

        # Check configuration names are unique
        if len(self.configurations) > 1:
            config_names = [c.name for c in self.configurations]
            if len(set(config_names)) != len(config_names):
                raise ValueError("Configuration names must be unique")

        # Check pool and reaction references
        for config in self.configurations:
            # Check input pool references
            for input_spec in config.inputs:
                if input_spec.pool not in pool_ids:
                    raise ValueError(
                        f"Input pool {input_spec.pool} in configuration "
                        f"{config.name} references unknown pool"
                    )

            # Check simulation variable references
            if config.simulation and config.simulation.variables:
                for flux_var in config.simulation.variables.flux_values:
                    if flux_var.flux not in reaction_ids:
                        raise ValueError(
                            f"Flux variable {flux_var.flux} in configuration "
                            f"{config.name} references unknown reaction"
                        )

                for pool_var in config.simulation.variables.poolsize_values:
                    if pool_var.pool not in pool_ids:
                        raise ValueError(
                            f"Pool size variable {pool_var.pool} in "
                            f"configuration {config.name} references "
                            f"unknown pool"
                        )

    @property
    def pool_ids(self) -> frozenset[str]:
        """Get all pool IDs in the model."""
        return self.reactionnetwork.metabolitepools.pool_ids

    @property
    def reaction_ids(self) -> frozenset[str]:
        """Get all reaction IDs in the model."""
        return self.reactionnetwork.reaction_ids

    @property
    def configuration_names(self) -> frozenset[str]:
        """Get all configuration names."""
        return frozenset(config.name for config in self.configurations)

    def get_configuration(self, name: str) -> Optional[Configuration]:
        """Get configuration by name."""
        for config in self.configurations:
            if config.name == name:
                return config
        return None

    def to_jax_representation(self) -> Dict[str, Any]:
        """
        Convert to JAX-compatible dictionary representation.

        Returns:
            Dictionary with JAX arrays for numerical computations
        """
        pool_ids = list(self.pool_ids)
        reaction_ids = list(self.reaction_ids)

        # Get stoichiometric matrix
        S = self.reactionnetwork.get_stoichiometric_matrix()

        # Get bounds matrices from first configuration (if available)
        flux_bounds = None
        poolsize_bounds = None

        if self.configurations:
            config = self.configurations[0]
            if config.simulation:
                flux_bounds, poolsize_bounds = (
                    config.simulation.get_optimization_bounds(
                        reaction_ids, pool_ids
                    )
                )

        if flux_bounds is None:
            flux_bounds = jnp.array([[-jnp.inf, jnp.inf]] * len(reaction_ids))
        if poolsize_bounds is None:
            poolsize_bounds = jnp.array([[0.0, jnp.inf]] * len(pool_ids))

        return {
            "stoichiometric_matrix": S,
            "flux_bounds": flux_bounds,
            "poolsize_bounds": poolsize_bounds,
            "pool_ids": pool_ids,
            "reaction_ids": reaction_ids,
            "n_pools": len(pool_ids),
            "n_reactions": len(reaction_ids),
            "n_configurations": len(self.configurations),
        }

    def get_tracer_experiment_data(
        self, config_name: str
    ) -> Optional[Dict[str, Any]]:
        """
        Get tracer experiment data for a specific configuration.

        Returns:
            Dictionary with JAX arrays for tracer experiment analysis
        """
        config = self.get_configuration(config_name)
        if not config:
            return None

        pool_ids = list(self.pool_ids)

        # Get tracer composition matrix
        tracer_matrix = config.get_tracer_composition_matrix(pool_ids)

        # Get measurement data
        measurement_data = None
        if config.measurement:
            values = config.measurement.data.values
            errors = config.measurement.data.errors
            times = config.measurement.data.times

            measurement_data = {
                "values": values,
                "errors": errors,
                "times": times,
            }

        return {
            "tracer_composition": tracer_matrix,
            "measurement_data": measurement_data,
            "stationary": config.stationary,
            "time_point": config.time,
        }
