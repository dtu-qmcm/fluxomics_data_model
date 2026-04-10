"""
FluxML reaction definitions.
"""

from typing import Optional, List, Dict
from pydantic import BaseModel, Field
import jax.numpy as jnp
from ..core.common import Annotation, JAXArray


class Reaction(BaseModel):
    """
    Represents a biochemical reaction in the metabolic network.

    A reaction describes the transformation of reactants to products.
    Reactions can have multiple atom mapping variants when symmetric
    compounds create ambiguity in carbon atom transitions.

    Attributes:
        id: Base reaction identifier (e.g., "SCS" for variants
            "SCS___1", "SCS___2")
        name: Optional human-readable name
        reversibility: Whether the reaction can proceed in both directions
        annotations: Additional metadata annotations
        reactants: List of reactant metabolite IDs (in order)
        products: List of product metabolite IDs (in order)
        atom_mapping_ids: Computational flux variable IDs for variants.
            None for single-map reactions, list of IDs for variant reactions.
            Example: ["bsDAP___1", "bsDAP___2", "bsDAP___3", "bsDAP___4"]

    Examples:
        Simple reaction without variants:
            Reaction(id="PGI", reactants=["G6P"], products=["F6P"])

        Reaction with 4 variants due to symmetric compounds:
            Reaction(id="bsDAP", reactants=["ASA", "PYR"],
                    products=["DAP", "H2O"],
                    atom_mapping_ids=["bsDAP___1", "bsDAP___2",
                                     "bsDAP___3", "bsDAP___4"])
    """

    id: str = Field(description="Reaction identifier (base name for variants)")
    name: Optional[str] = Field(default=None, description="Reaction name")
    reversibility: bool = Field(
        default=True, description="Reaction reversibility"
    )
    annotations: List[Annotation] = Field(
        default_factory=list, description="Annotations"
    )
    reactants: List[str] = Field(
        default_factory=list, description="Reactant metabolite IDs"
    )
    products: List[str] = Field(
        default_factory=list, description="Product metabolite IDs"
    )
    atom_mapping_ids: Optional[List[str]] = Field(
        default=None,
        description="Computational flux variable IDs for variants "
        "(e.g., ['SCS___1', 'SCS___2'])",
    )

    # JAX-compatible numerical representation
    stoichiometry_dict: Optional[Dict[str, float]] = Field(
        default=None, exclude=True, description="Stoichiometry"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def reactant_ids(self) -> frozenset[str]:
        """Get all reactant metabolite IDs."""
        return frozenset(self.reactants)

    @property
    def product_ids(self) -> frozenset[str]:
        """Get all product metabolite IDs."""
        return frozenset(self.products)

    @property
    def metabolites(self) -> frozenset[str]:
        """Get all participating metabolite IDs."""
        return self.reactant_ids | self.product_ids

    @property
    def is_variant_reaction(self) -> bool:
        """Check if this reaction has multiple atom map variants."""
        return (
            self.atom_mapping_ids is not None and len(self.atom_mapping_ids) > 1
        )

    @property
    def n_variants(self) -> int:
        """Number of atom map variants for this reaction."""
        return len(self.atom_mapping_ids) if self.atom_mapping_ids else 1

    @property
    def flux_bounds(self) -> jnp.ndarray:
        """Get flux bounds as JAX array [lower, upper]."""
        if self.flux_bounds_array is None:
            if self.reversibility:
                return jnp.array([-1000.0, 1000.0])
            else:
                return jnp.array([0.0, 1000.0])
        return self.flux_bounds_array.to_jax_array()

    @property
    def equation(self) -> str:
        """Get reaction equation string."""
        reactants = " + ".join(self.reactants)
        products = " + ".join(self.products)
        arrow = " <=> " if self.reversibility else " => "
        return f"{reactants}{arrow}{products}"

    def with_flux_bounds(self, lower: float, upper: float) -> "Reaction":
        """Create new reaction with specified flux bounds."""
        bounds_array = jnp.array([lower, upper])
        return Reaction(
            id=self.id,
            reversibility=self.reversibility,
            annotations=self.annotations,
            reactants=self.reactants,
            products=self.products,
            flux_bounds_array=JAXArray.from_jax_array(bounds_array),
            stoichiometry_dict=self.stoichiometry_dict,
        )

    def with_stoichiometry(self, stoichiometry: Dict[str, float]) -> "Reaction":
        """Create new reaction with specified stoichiometry."""
        return Reaction(
            id=self.id,
            reversibility=self.reversibility,
            annotations=self.annotations,
            reactants=self.reactants,
            products=self.products,
            flux_bounds_array=self.flux_bounds_array,
            stoichiometry_dict=stoichiometry,
        )

    def get_stoichiometric_vector(
        self, metabolite_ids: List[str]
    ) -> jnp.ndarray:
        """
        Get stoichiometric vector for this reaction.

        Args:
            metabolite_ids: Ordered list of metabolite IDs

        Returns:
            JAX array with stoichiometric coefficients
        """
        if self.stoichiometry_dict is None:
            # Default stoichiometry: -1 for reactants, +1 for products
            coefficients = []
            for metabolite_id in metabolite_ids:
                if metabolite_id in self.reactant_ids:
                    coefficients.append(-1.0)
                elif metabolite_id in self.product_ids:
                    coefficients.append(1.0)
                else:
                    coefficients.append(0.0)
            return jnp.array(coefficients)
        else:
            return jnp.array(
                [
                    self.stoichiometry_dict.get(metabolite_id, 0.0)
                    for metabolite_id in metabolite_ids
                ]
            )
