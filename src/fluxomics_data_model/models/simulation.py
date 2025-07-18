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

    Corresponds to fluxml/configuration/simulation/variables/fluxvalue
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


class PoolSizeValue(BaseModel):
    """
    FluxML pool size value specification for simulation.

    Corresponds to fluxml/configuration/simulation/variables/poolsizevalue
    """

    pool: str = Field(description="Pool identifier")
    value: Optional[float] = Field(default=None, description="Pool size value")
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
        """Get pool size bounds."""
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

    Corresponds to fluxml/configuration/simulation/variables
    """

    flux_values: List[FluxValue] = Field(
        default_factory=list, description="Flux variables"
    )
    poolsize_values: List[PoolSizeValue] = Field(
        default_factory=list, description="Pool size variables"
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

    @field_validator("poolsize_values")
    @classmethod
    def validate_unique_pool_ids(
        cls, v: List[PoolSizeValue]
    ) -> List[PoolSizeValue]:
        """Validate pool IDs are unique."""
        ids = [pv.pool for pv in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Pool size variable pool IDs must be unique")
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

    def get_poolsize_bounds_matrix(self, pool_ids: List[str]) -> jnp.ndarray:
        """
        Get pool size bounds matrix for JAX computations.

        Returns:
            JAX array of shape (n_pools, 2) with [lower, upper] bounds
        """
        pool_dict = {pv.pool: pv for pv in self.poolsize_values}

        bounds = []
        for pool_id in pool_ids:
            pool_var = pool_dict.get(pool_id)
            if pool_var:
                bounds.append(pool_var.bounds_array)
            else:
                # Default bounds for pool sizes
                bounds.append(jnp.array([0.0, jnp.inf]))

        return jnp.stack(bounds)


class Simulation(BaseModel):
    """
    FluxML simulation specification.

    Corresponds to fluxml/configuration/simulation
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
        self, flux_ids: List[str], pool_ids: List[str]
    ) -> tuple[jnp.ndarray, jnp.ndarray]:
        """
        Get optimization bounds for JAX optimization.

        Returns:
            Tuple of (flux_bounds, poolsize_bounds) as JAX arrays
        """
        if self.variables is None:
            # Default bounds
            flux_bounds = jnp.array([[-jnp.inf, jnp.inf]] * len(flux_ids))
            poolsize_bounds = jnp.array([[0.0, jnp.inf]] * len(pool_ids))
        else:
            flux_bounds = self.variables.get_flux_bounds_matrix(flux_ids)
            poolsize_bounds = self.variables.get_poolsize_bounds_matrix(
                pool_ids
            )

        return flux_bounds, poolsize_bounds
