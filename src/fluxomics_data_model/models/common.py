"""
Common data structures used across FluxML models.
"""

from typing import Optional
from pydantic import BaseModel, Field
import jax.numpy as jnp


class Annotation(BaseModel):
    """
    FluxML annotation element for additional metadata.

    Corresponds to fluxml/reactionnetwork/metabolitepools/pool/annotation
    and fluxml/reactionnetwork/reaction/annotation
    """

    name: str = Field(description="Annotation name")
    content: Optional[str] = Field(
        default=None, description="Annotation content"
    )

    class Config:
        frozen = True


class TextualOrMath(BaseModel):
    """
    FluxML textual or MathML content.

    Used for constraints and mathematical expressions.
    """

    textual: Optional[str] = Field(
        default=None, description="Textual representation"
    )
    mathml: Optional[str] = Field(
        default=None, description="MathML representation"
    )

    class Config:
        frozen = True

    def __init__(self, **data):
        super().__init__(**data)
        if not self.textual and not self.mathml:
            raise ValueError("Either textual or mathml must be provided")


class ErrorModel(BaseModel):
    """
    FluxML error model for measurement uncertainties.

    Corresponds to fluxml/configuration/measurement/model/*/errormodel
    """

    expression: TextualOrMath = Field(description="Error model expression")

    class Config:
        frozen = True


class JAXArray(BaseModel):
    """
    JAX array wrapper for Pydantic serialization.
    """

    shape: tuple[int, ...] = Field(description="Array shape")
    dtype: str = Field(description="Array dtype")
    data: list = Field(description="Array data as nested lists")

    class Config:
        frozen = True

    @classmethod
    def from_jax_array(cls, arr: jnp.ndarray) -> "JAXArray":
        """Create from JAX array."""
        return cls(shape=arr.shape, dtype=str(arr.dtype), data=arr.tolist())

    def to_jax_array(self) -> jnp.ndarray:
        """Convert to JAX array."""
        return jnp.array(self.data, dtype=self.dtype).reshape(self.shape)


class TimeSeries(BaseModel):
    """
    Time series data structure for non-stationary measurements.
    """

    times: JAXArray = Field(description="Time points")
    values: JAXArray = Field(description="Measurement values")
    errors: Optional[JAXArray] = Field(
        default=None, description="Measurement errors"
    )

    class Config:
        frozen = True

    @classmethod
    def from_arrays(
        cls,
        times: jnp.ndarray,
        values: jnp.ndarray,
        errors: Optional[jnp.ndarray] = None,
    ) -> "TimeSeries":
        """Create from JAX arrays."""
        return cls(
            times=JAXArray.from_jax_array(times),
            values=JAXArray.from_jax_array(values),
            errors=(
                JAXArray.from_jax_array(errors) if errors is not None else None
            ),
        )

    def to_arrays(
        self,
    ) -> tuple[jnp.ndarray, jnp.ndarray, Optional[jnp.ndarray]]:
        """Convert to JAX arrays."""
        return (
            self.times.to_jax_array(),
            self.values.to_jax_array(),
            self.errors.to_jax_array() if self.errors else None,
        )
