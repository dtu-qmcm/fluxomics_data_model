"""
FluxML reaction definitions.
"""

from typing import Optional, List, Dict
from pydantic import BaseModel, Field
import jax.numpy as jnp
from .common import Annotation, JAXArray


class Variant(BaseModel):
    """
    FluxML reaction variant for alternative atom mappings.

    Corresponds to fluxml/reactionnetwork/reaction/(reduct|rproduct)/variant
    """

    cfg: str = Field(description="Atom configuration")
    ratio: Optional[float] = Field(default=None, description="Variant ratio")

    class Config:
        frozen = True
        extra = "forbid"


class Reduct(BaseModel):
    """
    FluxML reaction reactant (substrate).

    Corresponds to fluxml/reactionnetwork/reaction/reduct
    """

    id: str = Field(description="Pool ID reference")
    cfg: Optional[str] = Field(default=None, description="Atom configuration")
    variants: List[Variant] = Field(
        default_factory=list, description="Alternative mappings"
    )

    class Config:
        frozen = True
        extra = "forbid"


class RProduct(BaseModel):
    """
    FluxML reaction product.

    Corresponds to fluxml/reactionnetwork/reaction/rproduct
    """

    id: str = Field(description="Pool ID reference")
    cfg: Optional[str] = Field(default=None, description="Atom configuration")
    variants: List[Variant] = Field(
        default_factory=list, description="Alternative mappings"
    )

    class Config:
        frozen = True
        extra = "forbid"


class Reaction(BaseModel):
    """
    FluxML reaction definition.

    Corresponds to fluxml/reactionnetwork/reaction
    """

    id: str = Field(description="Reaction identifier")
    bidirectional: bool = Field(
        default=True, description="Reaction reversibility"
    )
    annotations: List[Annotation] = Field(
        default_factory=list, description="Annotations"
    )
    reducts: List[Reduct] = Field(default_factory=list, description="Reactants")
    rproducts: List[RProduct] = Field(
        default_factory=list, description="Products"
    )

    # JAX-compatible numerical representation
    flux_bounds_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Flux bounds"
    )
    stoichiometry_dict: Optional[Dict[str, float]] = Field(
        default=None, exclude=True, description="Stoichiometry"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def reactant_ids(self) -> frozenset[str]:
        """Get all reactant pool IDs."""
        return frozenset(reduct.id for reduct in self.reducts)

    @property
    def product_ids(self) -> frozenset[str]:
        """Get all product pool IDs."""
        return frozenset(rproduct.id for rproduct in self.rproducts)

    @property
    def participating_pools(self) -> frozenset[str]:
        """Get all participating pool IDs."""
        return self.reactant_ids | self.product_ids

    @property
    def flux_bounds(self) -> jnp.ndarray:
        """Get flux bounds as JAX array [lower, upper]."""
        if self.flux_bounds_array is None:
            if self.bidirectional:
                return jnp.array([-jnp.inf, jnp.inf])
            else:
                return jnp.array([0.0, jnp.inf])
        return self.flux_bounds_array.to_jax_array()

    def with_flux_bounds(self, lower: float, upper: float) -> "Reaction":
        """Create new reaction with specified flux bounds."""
        bounds_array = jnp.array([lower, upper])
        return Reaction(
            id=self.id,
            bidirectional=self.bidirectional,
            annotations=self.annotations,
            reducts=self.reducts,
            rproducts=self.rproducts,
            flux_bounds_array=JAXArray.from_jax_array(bounds_array),
            stoichiometry_dict=self.stoichiometry_dict,
        )

    def with_stoichiometry(self, stoichiometry: Dict[str, float]) -> "Reaction":
        """Create new reaction with specified stoichiometry."""
        return Reaction(
            id=self.id,
            bidirectional=self.bidirectional,
            annotations=self.annotations,
            reducts=self.reducts,
            rproducts=self.rproducts,
            flux_bounds_array=self.flux_bounds_array,
            stoichiometry_dict=stoichiometry,
        )

    def get_stoichiometric_vector(self, pool_ids: List[str]) -> jnp.ndarray:
        """
        Get stoichiometric vector for this reaction.

        Args:
            pool_ids: Ordered list of pool IDs

        Returns:
            JAX array with stoichiometric coefficients
        """
        if self.stoichiometry_dict is None:
            # Default stoichiometry: -1 for reactants, +1 for products
            coefficients = []
            for pool_id in pool_ids:
                if pool_id in self.reactant_ids:
                    coefficients.append(-1.0)
                elif pool_id in self.product_ids:
                    coefficients.append(1.0)
                else:
                    coefficients.append(0.0)
            return jnp.array(coefficients)
        else:
            return jnp.array(
                [
                    self.stoichiometry_dict.get(pool_id, 0.0)
                    for pool_id in pool_ids
                ]
            )

    @property
    def equation(self) -> str:
        """Get reaction equation string."""
        reactants = " + ".join(reduct.id for reduct in self.reducts)
        products = " + ".join(rproduct.id for rproduct in self.rproducts)
        arrow = " <=> " if self.bidirectional else " => "
        return f"{reactants}{arrow}{products}"
