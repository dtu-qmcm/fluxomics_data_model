"""Unit tests for the reaction_network module"""

from fluxomics_data_model.reaction_network import (
    Species,
    Reactant,
    Reaction,
    ReactionNetwork,
)


def test_species_creation():
    """Test creation of Species objects."""
    species = Species(id="glc", name="glucose", n_labellable_atom=6)
    assert species.id == "glc"
    assert species.name == "glucose"
    assert species.n_labellable_atom == 6


def test_reactant_creation():
    """Test creation of Reactant objects."""
    reactant = Reactant(
        id="glc", is_product=False, atom_pattern="ABCDEF", name="glucose"
    )
    assert reactant.id == "glc"
    assert reactant.is_product is False
    assert reactant.atom_pattern == "ABCDEF"
    assert reactant.name == "glucose"
    assert reactant.is_substrate is True


def test_reaction_creation():
    """Test creation of Reaction objects."""
    substrate = Reactant(
        id="g6p",
        is_product=False,
        atom_pattern="ABCDEF",
        name="glucose 6-phosphate",
    )
    product = Reactant(
        id="f6p",
        is_product=True,
        atom_pattern="ABCDEF",
        name="fructose 6-phosphate",
    )
    reaction = Reaction(
        id="r1",
        bidirectional=True,
        reactants=[substrate, product],
        name="glucose to pyruvate",
        comment="glycolysis reaction",
    )

    assert reaction.id == "r1"
    assert reaction.bidirectional is True
    assert len(reaction.reactants) == 2
    assert len(reaction.substrates) == 1
    assert len(reaction.products) == 1
    assert reaction.name == "glucose to pyruvate"
    assert reaction.comment == "glycolysis reaction"


def test_reaction_network_creation():
    """Test creation of ReactionNetwork objects."""
    glc = Species(id="glc", name="glucose", n_labellable_atom=6)
    pyr = Species(id="pyr", name="pyruvate", n_labellable_atom=3)

    r1_glc = Reactant(
        id="glc",
        is_product=False,
        atom_pattern="abcdef",
        name="glucose in reaction r1",
    )
    r1_pyr = Reactant(
        id="pyr",
        is_product=True,
        atom_pattern="abc",
        name="pyruvate in reaction r1",
    )

    r1 = Reaction(
        id="r1",
        bidirectional=True,
        reactants=[r1_glc, r1_pyr],
        name="glucose to pyruvate",
    )

    network = ReactionNetwork(reactions=[r1], species=[glc, pyr])

    assert len(network.reactions) == 1
    assert len(network.species) == 2


def test_optional_fields():
    """Test handling of optional fields."""
    species = Species(id="glc", n_labellable_atom=6)
    assert species.name is None

    reactant = Reactant(id="glc", is_product=False, atom_pattern="ABCDEF")
    assert reactant.name is None

    reaction = Reaction(
        id="r1",
        bidirectional=True,
        reactants=[],
        name="test reaction",
    )
    assert reaction.comment is None
