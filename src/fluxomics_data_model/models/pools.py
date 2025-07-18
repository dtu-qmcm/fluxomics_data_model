"""
FluxML metabolite pool definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .common import Annotation, JAXArray


class Pool(BaseModel):
    """
    FluxML metabolite pool definition.

    Corresponds to fluxml/reactionnetwork/metabolitepools/pool
    """

    id: str = Field(description="Pool identifier")
    atoms: int = Field(default=0, ge=0, le=1024, description="Number of atoms")
    size: float = Field(default=1.0, description="Pool size")
    cfg: str = Field(default="0", description="Atom configuration")
    annotations: List[Annotation] = Field(
        default_factory=list, description="Annotations"
    )

    # JAX-compatible numerical representation
    jax_atoms: Optional[JAXArray] = Field(
        default=None, exclude=True, description="JAX atom array"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("cfg")
    @classmethod
    def validate_cfg(cls, v: str) -> str:
        """Validate atom configuration string."""
        if not v or v == "0":
            return v
        # Basic validation - full implementation would parse atom mapping
        return v

    @property
    def atom_vector(self) -> jnp.ndarray:
        """Get JAX array representation of atoms."""
        if self.jax_atoms is None:
            # Create atom vector based on configuration
            if self.atoms == 0:
                return jnp.array([])
            return jnp.ones(self.atoms, dtype=jnp.float32)
        return self.jax_atoms.to_jax_array()

    def with_jax_atoms(self, atoms: jnp.ndarray) -> "Pool":
        """Create new pool with JAX atom array."""
        return Pool(
            id=self.id,
            atoms=self.atoms,
            size=self.size,
            cfg=self.cfg,
            annotations=self.annotations,
            jax_atoms=JAXArray.from_jax_array(atoms),
        )


class MetabolitePools(BaseModel):
    """
    FluxML metabolite pools collection.

    Corresponds to fluxml/reactionnetwork/metabolitepools
    """

    pools: List[Pool] = Field(
        min_length=1, description="List of metabolite pools"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("pools")
    @classmethod
    def validate_unique_ids(cls, v: List[Pool]) -> List[Pool]:
        """Validate pool IDs are unique."""
        ids = [pool.id for pool in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Pool IDs must be unique")
        return v

    @property
    def pool_dict(self) -> dict[str, Pool]:
        """Get pools as dictionary keyed by ID."""
        return {pool.id: pool for pool in self.pools}

    @property
    def pool_ids(self) -> frozenset[str]:
        """Get all pool IDs."""
        return frozenset(pool.id for pool in self.pools)

    def get_pool(self, pool_id: str) -> Optional[Pool]:
        """Get pool by ID."""
        return self.pool_dict.get(pool_id)

    def to_stoichiometric_matrix(self, reaction_ids: List[str]) -> jnp.ndarray:
        """
        Create stoichiometric matrix for JAX computations.

        Returns:
            JAX array of shape (n_pools, n_reactions)
        """
        # This would be implemented with actual stoichiometry data
        n_pools = len(self.pools)
        n_reactions = len(reaction_ids)
        return jnp.zeros((n_pools, n_reactions), dtype=jnp.float32)
