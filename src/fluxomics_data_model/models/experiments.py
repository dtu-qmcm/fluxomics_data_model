"""
FluxML experimental setup definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .common import TextualOrMath, JAXArray, TimeSeries
from .constraints import Constraints
from .measurements import Measurement
from .simulation import Simulation


class Label(BaseModel):
    """
    FluxML isotope label specification.

    Corresponds to fluxml/experiments/tracers/label
    """

    cfg: str = Field(
        description="Label configuration pattern", pattern=r"[01xX]+"
    )
    purity: Optional[str] = Field(default=None, description="Label purity")
    cost: Optional[float] = Field(default=None, description="Label cost")
    content: Optional[str] = Field(default=None, description="Label content")
    expression: Optional[TextualOrMath] = Field(
        default=None, description="Mathematical expression"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("purity")
    @classmethod
    def validate_purity(cls, v: Optional[str]) -> Optional[str]:
        """Validate purity string format."""
        if v is None:
            return v
        # Basic validation - full implementation would parse purity values
        return v


class Tracers(BaseModel):
    """
    FluxML tracer specification for tracer experiments.

    Corresponds to fluxml/experiments/tracers
    """

    id: Optional[str] = Field(default=None, description="Tracer identifier")
    metabolite: str = Field(description="Metabolite identifier")
    type: str = Field(default="isotopomer", description="Tracer type")
    profile: Optional[str] = Field(default=None, description="Time profile")
    labels: List[Label] = Field(  # noqa: B008
        default_factory=list, description="Isotope labels"
    )

    # JAX-compatible numerical representation
    composition_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Composition vector"
    )
    time_profile_data: Optional[TimeSeries] = Field(
        default=None, exclude=True, description="Time profile data"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def composition_vector(self) -> jnp.ndarray:
        """Get JAX array representation of tracer composition."""
        if self.composition_array is None:
            # Default: create composition from labels
            if not self.labels:
                return jnp.array([1.0])  # Unlabeled
            # Simple implementation - would need full label parsing
            return jnp.array([1.0])
        return self.composition_array.to_jax_array()

    def with_composition(self, composition: jnp.ndarray) -> "Tracers":
        """Create new tracer with specified composition."""
        return Tracers(
            id=self.id,
            metabolite=self.metabolite,
            type=self.type,
            profile=self.profile,
            labels=self.labels,
            composition_array=JAXArray.from_jax_array(composition),
            time_profile_data=self.time_profile_data,
        )

    def with_time_profile(
        self, times: jnp.ndarray, values: jnp.ndarray
    ) -> "Tracers":
        """Create new tracer with time profile."""
        return Tracers(
            id=self.id,
            metabolite=self.metabolite,
            type=self.type,
            profile=self.profile,
            labels=self.labels,
            composition_array=self.composition_array,
            time_profile_data=TimeSeries.from_arrays(times, values),
        )


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
