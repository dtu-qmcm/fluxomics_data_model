"""
FluxML simulation definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .measurements import MeasurementModel


class FluxValue(BaseModel):
    """
    FluxML flux value specification for simulation.

    Corresponds to fluxml/experiments/simulation/variables/fluxvalue
    """

    flux: str = Field(description="Flux identifier")
    type: str = Field(description="Flux type")
    value: Optional[float] = Field(default=None, description="Flux value")
    lo: Optional[float] = Field(default=None, description="Lower bound")
    hi: Optional[float] = Field(default=None, description="Upper bound")
    inc: Optional[float] = Field(default=None, description="Increment")
    edweight: float = Field(
        default=1.0, ge=0.0, le=1.0, description="ED weight"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def bounds(self) -> tuple[Optional[float], Optional[float]]:
        """Get flux bounds."""
        return (self.lo, self.hi)

    @property
    def bounds_array(self) -> jnp.ndarray:
        """Get bounds as JAX array."""
        lo = self.lo if self.lo is not None else -jnp.inf
        hi = self.hi if self.hi is not None else jnp.inf
        return jnp.array([lo, hi])


class MetaboliteSizeValue(BaseModel):
    """
    FluxML metabolite size value specification for simulation.

    Corresponds to fluxml/experiments/simulation/variables/metabolitesizevalue
    """

    metabolite: str = Field(description="Metabolite identifier")
    value: Optional[float] = Field(
        default=None, description="Metabolite size value"
    )
    lo: Optional[float] = Field(default=None, description="Lower bound")
    hi: Optional[float] = Field(default=None, description="Upper bound")
    inc: Optional[float] = Field(default=None, description="Increment")
    edweight: float = Field(
        default=1.0, ge=0.0, le=1.0, description="ED weight"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def bounds(self) -> tuple[Optional[float], Optional[float]]:
        """Get metabolite size bounds."""
        return (self.lo, self.hi)

    @property
    def bounds_array(self) -> jnp.ndarray:
        """Get bounds as JAX array."""
        lo = self.lo if self.lo is not None else 0.0
        hi = self.hi if self.hi is not None else jnp.inf
        return jnp.array([lo, hi])


class Variables(BaseModel):
    """
    FluxML simulation variables.

    Corresponds to fluxml/experiments/simulation/variables
    """

    flux_values: List[FluxValue] = Field(
        default_factory=list, description="Flux variables"
    )
    metabolitesize_values: List[MetaboliteSizeValue] = Field(
        default_factory=list, description="Metabolite size variables"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("flux_values")
    @classmethod
    def validate_unique_flux_ids(cls, v: List[FluxValue]) -> List[FluxValue]:
        """Validate flux IDs are unique."""
        ids = [(fv.flux, fv.type) for fv in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Flux variable (flux, type) pairs must be unique")
        return v

    @field_validator("metabolitesize_values")
    @classmethod
    def validate_unique_metabolite_ids(
        cls, v: List[MetaboliteSizeValue]
    ) -> List[MetaboliteSizeValue]:
        """Validate metabolite IDs are unique."""
        ids = [pv.metabolite for pv in v]
        if len(set(ids)) != len(ids):
            raise ValueError(
                "Metabolite size variable metabolite IDs must be unique"
            )
        return v

    def get_flux_bounds_matrix(self, flux_ids: List[str]) -> jnp.ndarray:
        """
        Get flux bounds matrix for JAX computations.

        Returns:
            JAX array of shape (n_fluxes, 2) with [lower, upper] bounds
        """
        flux_dict = {(fv.flux, fv.type): fv for fv in self.flux_values}

        bounds = []
        for flux_id in flux_ids:
            # Check both net and xch types
            net_var = flux_dict.get((flux_id, "net"))
            xch_var = flux_dict.get((flux_id, "xch"))

            if net_var:
                bounds.append(net_var.bounds_array)
            elif xch_var:
                bounds.append(xch_var.bounds_array)
            else:
                # Default bounds
                bounds.append(jnp.array([-jnp.inf, jnp.inf]))

        return jnp.stack(bounds)

    def get_metabolitesize_bounds_matrix(
        self, metabolite_ids: List[str]
    ) -> jnp.ndarray:
        """
        Get metabolite size bounds matrix for JAX computations.

        Returns:
            JAX array of shape (n_metabolites, 2) with [lower, upper] bounds
        """
        metabolite_dict = {
            pv.metabolite: pv for pv in self.metabolitesize_values
        }

        bounds = []
        for metabolite_id in metabolite_ids:
            metabolite_var = metabolite_dict.get(metabolite_id)
            if metabolite_var:
                bounds.append(metabolite_var.bounds_array)
            else:
                # Default bounds for metabolite sizes
                bounds.append(jnp.array([0.0, jnp.inf]))

        return jnp.stack(bounds)


class Simulation(BaseModel):
    """
    FluxML simulation specification.

    Corresponds to fluxml/experiments/simulation
    """

    type: str = Field(default="auto", description="Simulation type")
    method: str = Field(default="auto", description="Simulation method")
    model: Optional[MeasurementModel] = Field(
        default=None, description="Simulation model"
    )
    variables: Optional[Variables] = Field(
        default=None, description="Simulation variables"
    )

    class Config:
        frozen = True
        extra = "forbid"

    def get_optimization_bounds(
        self, flux_ids: List[str], metabolite_ids: List[str]
    ) -> tuple[jnp.ndarray, jnp.ndarray]:
        """
        Get optimization bounds for JAX optimization.

        Returns:
            Tuple of (flux_bounds, metabolitesize_bounds) as JAX arrays
        """
        if self.variables is None:
            # Default bounds
            flux_bounds = jnp.array([[-jnp.inf, jnp.inf]] * len(flux_ids))
            metabolitesize_bounds = jnp.array(
                [[0.0, jnp.inf]] * len(metabolite_ids)
            )
        else:
            flux_bounds = self.variables.get_flux_bounds_matrix(flux_ids)
            metabolitesize_bounds = (
                self.variables.get_metabolitesize_bounds_matrix(metabolite_ids)
            )

        return flux_bounds, metabolitesize_bounds
