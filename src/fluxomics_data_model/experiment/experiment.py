"""
FluxML experimental setup definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field
import jax.numpy as jnp
from .tracer import Tracers
from ..model.constraint import Constraints
from .measurement import Measurement
from .simulation import Simulation


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
