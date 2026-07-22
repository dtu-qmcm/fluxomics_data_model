"""FluxML tracer and isotope label composition definitions.

A *tracer* specifies the isotopically labelled substrate fed to the cell in
a ¹³C MFA experiment.  The two classes here map to the FluxML
``<input>`` element:

- :class:`LabelComposition` — one ``<label>`` entry: a specific isotopomer
  pattern (e.g. ``"111111"`` for uniformly labelled glucose), its fractional
  abundance in the tracer mixture, and optional purity and cost metadata.
- :class:`Tracers` — the full ``<input>`` block for one substrate pool:
  a list of :class:`LabelComposition` entries plus the metabolite ID,
  tracer type (``"isotopomer"``, ``"cumomer"``, or ``"emu"``), and an
  optional time-varying labelling profile.

JAX integration
---------------
Call :meth:`Tracers.with_composition` to attach a pre-computed composition
vector (a JAX array) and :meth:`Tracers.with_time_profile` to attach
time-course data.  These methods return new immutable instances; the
original object is unchanged.

Note: :attr:`Tracers.composition_vector` currently returns a placeholder
``jnp.array([1.0])`` when no composition array has been explicitly attached
via :meth:`with_composition`.  Numerical code that needs the actual isotopomer
fractions must call :meth:`with_composition` first.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator, model_validator
import jax.numpy as jnp
from ..core.common import TextualOrMath, JAXArray, TimeSeries


class LabelComposition(BaseModel):
    """FluxML isotope label composition specification.

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
    """FluxML tracer specification for tracer experiments.

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
        """Return the tracer composition as a JAX array.

        Resolution order:

        1. If a composition array has been explicitly attached via
           :meth:`with_composition`, that array is returned unchanged.
        2. Otherwise, the vector is derived from ``self.labels``:
           each :class:`LabelComposition` whose ``fraction`` field is a
           plain ``float`` contributes one element.  The resulting vector
           contains those fractions in the order the labels appear.
        3. If no labels are defined, or none carry a numeric ``fraction``
           (e.g. all labels use time-varying ``expression`` values), a
           fallback of ``jnp.array([1.0])`` is returned, representing a
           fully unlabelled substrate.

        Examples::

            # Labels with fractions → vector built from labels
            t = Tracers(
                metabolite="Glc",
                labels=[
                    LabelComposition(labeled_pattern="111111", fraction=0.2),
                    LabelComposition(labeled_pattern="000000", fraction=0.8),
                ],
            )
            t.composition_vector  # jnp.array([0.2, 0.8])

            # Explicitly set vector takes priority
            t2 = t.with_composition(jnp.array([0.5, 0.5]))
            t2.composition_vector  # jnp.array([0.5, 0.5])

            # No labels → fallback
            Tracers(metabolite="Glc").composition_vector  # jnp.array([1.0])

        Returns:
            1-D JAX array of fractional abundances.
        """
        if self.composition_array is not None:
            return self.composition_array.to_jax_array()

        # Derive from labels: collect numeric fractions in order.
        fractions = [
            label.fraction
            for label in self.labels
            if label.fraction is not None and isinstance(label.fraction, float)
        ]

        if fractions:
            return jnp.array(fractions)

        # Fallback: no labels, or all labels use expression-only profiles.
        return jnp.array([1.0])

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
