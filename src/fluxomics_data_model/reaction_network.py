"""Format for a reaction network."""

from pydantic import BaseModel, NonNegativeInt


class Species(BaseModel):
    """A species"""

    id: str
    name: str | None = None
    n_labellable_atom: NonNegativeInt


class Reactant(BaseModel):
    """A reactant"""

    id: str
    is_product: bool
    atom_pattern: str
    name: str | None = None

    @property
    def is_substrate(self) -> bool:
        return not self.is_product


class Reaction(BaseModel):
    """A reaction"""

    id: str
    bidirectional: bool
    reactants: list[Reactant]
    name: str
    comment: str | None = None

    @property
    def substrates(self) -> list[Reactant]:
        return [r for r in self.reactants if r.is_substrate]

    @property
    def products(self) -> list[Reactant]:
        return [r for r in self.reactants if r.is_product]


class ReactionNetwork(BaseModel):
    """A class representing the reactions in a metabolic network."""

    reactions: list[Reaction]
    species: list[Species]
