"""
FluxML metabolite definitions.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .common import Annotation, JAXArray


class Metabolite(BaseModel):
    """
    FluxML metabolite definition.

    Corresponds to fluxml/reactionnetwork/metabolites/metabolite
    """

    id: str = Field(description="Metabolite identifier")
    atoms: int = Field(default=0, ge=0, le=1024, description="Number of atoms")
    size: float = Field(default=1.0, description="Metabolite size")
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

    def with_jax_atoms(self, atoms: jnp.ndarray) -> "Metabolite":
        """Create new metabolite with JAX atom array."""
        return Metabolite(
            id=self.id,
            atoms=self.atoms,
            size=self.size,
            cfg=self.cfg,
            annotations=self.annotations,
            jax_atoms=JAXArray.from_jax_array(atoms),
        )


class Metabolites(BaseModel):
    """
    FluxML metabolites collection.

    Corresponds to fluxml/reactionnetwork/metabolites
    """

    metabolites: List[Metabolite] = Field(
        min_length=1, description="List of metabolites"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("metabolites")
    @classmethod
    def validate_unique_ids(cls, v: List[Metabolite]) -> List[Metabolite]:
        """Validate metabolite IDs are unique."""
        ids = [metabolite.id for metabolite in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Metabolite IDs must be unique")
        return v

    @property
    def metabolite_dict(self) -> dict[str, Metabolite]:
        """Get metabolites as dictionary keyed by ID."""
        return {metabolite.id: metabolite for metabolite in self.metabolites}

    @property
    def metabolite_ids(self) -> frozenset[str]:
        """Get all metabolite IDs."""
        return frozenset(metabolite.id for metabolite in self.metabolites)

    def get_metabolite(self, metabolite_id: str) -> Optional[Metabolite]:
        """Get metabolite by ID."""
        return self.metabolite_dict.get(metabolite_id)

    def to_stoichiometric_matrix(self, reaction_ids: List[str]) -> jnp.ndarray:
        """
        Create stoichiometric matrix for JAX computations.

        Returns:
            JAX array of shape (n_metabolites, n_reactions)
        """
        # This would be implemented with actual stoichiometry data
        n_metabolites = len(self.metabolites)
        n_reactions = len(reaction_ids)
        return jnp.zeros((n_metabolites, n_reactions), dtype=jnp.float32)
