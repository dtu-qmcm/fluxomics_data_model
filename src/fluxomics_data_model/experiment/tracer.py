"""
FluxML tracer and label composition definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator, model_validator
import jax.numpy as jnp
from ..core.common import TextualOrMath, JAXArray, TimeSeries


class LabelComposition(BaseModel):
    """
    FluxML isotope label composition specification.

    Corresponds to fluxml/experiments/tracers/label
    """

    labeled_pattern: str = Field(
        description="Label configuration pattern", pattern=r"[01xX]+"
    )
    purity: Optional[float] = Field(
        default=None,
        description="Label purity (between 0 and 1)",
        ge=0.0,
        le=1.0,
    )
    cost: Optional[float] = Field(
        default=None, description="Label cost (must be positive)", gt=0.0
    )
    fraction: Optional[float] = Field(
        default=None,
        description="Label fraction (for float fractions, must sum to 1)",
    )
    expression: Optional[TextualOrMath] = Field(
        default=None, description="Mathematical expression"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("purity")
    @classmethod
    def validate_purity(cls, v: Optional[float]) -> Optional[float]:
        """Validate purity is between 0 and 1."""
        if v is None:
            return v
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"Purity must be between 0 and 1, got {v}")
        return v

    @field_validator("cost")
    @classmethod
    def validate_cost(cls, v: Optional[float]) -> Optional[float]:
        """Validate cost is positive."""
        if v is None:
            return v
        if v <= 0.0:
            raise ValueError(f"Cost must be positive, got {v}")
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
    labels: List[LabelComposition] = Field(
        default_factory=list, description="Isotope label compositions"
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

    @model_validator(mode="after")
    def validate_fractions(self) -> "Tracers":
        """Validate that float fractions sum to approximately 1.0.

        Allows fractions < 1.0 when natural abundance is not explicitly
        specified (common in influx_si .linp files where the unlabeled
        remainder is implied). Raises if sum exceeds 1.0.
        """
        import warnings

        float_fractions = []
        has_float_fractions = False

        for label in self.labels:
            if label.fraction is not None and isinstance(label.fraction, float):
                has_float_fractions = True
                float_fractions.append(label.fraction)
            elif label.fraction is not None:
                return self

        if has_float_fractions and float_fractions:
            total = sum(float_fractions)
            if total > 1.0 + 1e-6:
                raise ValueError(
                    f"Label fractions sum to {total}, which exceeds 1.0"
                )
            if abs(total - 1.0) > 1e-6:
                warnings.warn(
                    f"Label fractions sum to {total:.6f} (expected 1.0). "
                    f"The remainder is assumed to be natural abundance.",
                    stacklevel=2,
                )

        return self

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
