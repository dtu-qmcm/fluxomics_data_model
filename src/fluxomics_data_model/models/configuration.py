"""
FluxML configuration definitions.
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

    Corresponds to fluxml/configuration/input/label
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


class Input(BaseModel):
    """
    FluxML input specification for tracer experiments.

    Corresponds to fluxml/configuration/input
    """

    id: Optional[str] = Field(default=None, description="Input identifier")
    pool: str = Field(description="Pool identifier")
    type: str = Field(default="isotopomer", description="Input type")
    profile: Optional[str] = Field(default=None, description="Time profile")
    labels: List[Label] = Field(
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

    def with_composition(self, composition: jnp.ndarray) -> "Input":
        """Create new input with specified composition."""
        return Input(
            id=self.id,
            pool=self.pool,
            type=self.type,
            profile=self.profile,
            labels=self.labels,
            composition_array=JAXArray.from_jax_array(composition),
            time_profile_data=self.time_profile_data,
        )

    def with_time_profile(
        self, times: jnp.ndarray, values: jnp.ndarray
    ) -> "Input":
        """Create new input with time profile."""
        return Input(
            id=self.id,
            pool=self.pool,
            type=self.type,
            profile=self.profile,
            labels=self.labels,
            composition_array=self.composition_array,
            time_profile_data=TimeSeries.from_arrays(times, values),
        )


class Configuration(BaseModel):
    """
    FluxML experimental configuration.

    Corresponds to fluxml/configuration
    """

    name: str = Field(description="Configuration name")
    stationary: bool = Field(default=True, description="Stationary assumption")
    time: Optional[float] = Field(default=None, description="Time point")
    comment: Optional[str] = Field(
        default=None, description="Configuration comment"
    )
    inputs: List[Input] = Field(
        default_factory=list, description="Input specifications"
    )
    constraints: Optional[Constraints] = Field(
        default=None, description="Configuration constraints"
    )
    measurement: Optional[Measurement] = Field(
        default=None, description="Measurement data"
    )
    simulation: Optional[Simulation] = Field(
        default=None, description="Simulation settings"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def input_pools(self) -> frozenset[str]:
        """Get all input pool IDs."""
        return frozenset(input_spec.pool for input_spec in self.inputs)

    def get_input_for_pool(self, pool_id: str) -> Optional[Input]:
        """Get input specification for a pool."""
        for input_spec in self.inputs:
            if input_spec.pool == pool_id:
                return input_spec
        return None

    def get_tracer_composition_matrix(self, pool_ids: List[str]) -> jnp.ndarray:
        """
        Get tracer composition matrix for JAX computations.

        Returns:
            JAX array of shape (n_pools, n_isotopomers)
        """
        compositions = []
        for pool_id in pool_ids:
            input_spec = self.get_input_for_pool(pool_id)
            if input_spec:
                compositions.append(input_spec.composition_vector)
            else:
                compositions.append(jnp.array([1.0]))  # Unlabeled

        # Pad to same length
        max_len = max(len(comp) for comp in compositions)
        padded = [
            jnp.pad(comp, (0, max_len - len(comp))) for comp in compositions
        ]

        return jnp.stack(padded)
