"""Provides the class FluxomicsDataset."""

from pydantic import BaseModel
from fluxomics_data_model.reaction_network import ReactionNetwork


class FluxomicsDatasetMetadata(BaseModel):
    """Metadata for a fluxomics dataset."""

    fluxomics_data_model_version: str
    description: str


class FluxomicsDataset(BaseModel):
    """A fluxomics dataset"""

    metadata: FluxomicsDatasetMetadata
    reaction_network: ReactionNetwork
