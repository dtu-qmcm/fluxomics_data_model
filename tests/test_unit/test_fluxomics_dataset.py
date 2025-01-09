"""Unit tests for the fluxomics_dataset module."""

from fluxomics_data_model.fluxomics_dataset import (
    FluxomicsDataset,
    FluxomicsDatasetMetadata,
)
from fluxomics_data_model.reaction_network import (
    Species,
    Reactant,
    Reaction,
    ReactionNetwork,
)


def test_fluxomics_dataset_metadata_creation():
    """Test creation of FluxomicsDatasetMetadata objects."""
    metadata = FluxomicsDatasetMetadata(
        fluxomics_data_model_version="1.0.0", description="Test dataset"
    )
    assert metadata.fluxomics_data_model_version == "1.0.0"
    assert metadata.description == "Test dataset"


def test_fluxomics_dataset_creation():
    """Test creation of FluxomicsDataset objects."""
    # Create metadata
    metadata = FluxomicsDatasetMetadata(
        fluxomics_data_model_version="1.0.0", description="Test dataset"
    )

    # Create a simple reaction network
    g6p = Species(id="g6p", name="glucose 6-phosphate", n_labellable_atom=6)
    f6p = Species(id="f6p", name="fructose 6-phosphate", n_labellable_atom=6)
    pgi_g6p = Reactant(
        id="g6p",
        is_product=False,
        atom_pattern="abcdef",
    )
    pgi_f6p = Reactant(
        id="f6p",
        is_product=True,
        atom_pattern="abcdef",
    )
    pgi = Reaction(
        id="pgi",
        bidirectional=True,
        reactants=[pgi_g6p, pgi_f6p],
        name="Glucose-6-phosphate isomerase",
    )
    network = ReactionNetwork(reactions=[pgi], species=[g6p, f6p])

    # Create dataset
    dataset = FluxomicsDataset(metadata=metadata, reaction_network=network)

    assert dataset.metadata == metadata
    assert dataset.reaction_network == network
    assert dataset.metadata.fluxomics_data_model_version == "1.0.0"
    assert len(dataset.reaction_network.reactions) == 1
    assert len(dataset.reaction_network.species) == 2
