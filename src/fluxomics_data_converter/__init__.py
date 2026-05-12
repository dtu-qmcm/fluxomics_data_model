"""FluxML Data Model.

A JAX-compatible Python library for handling FluxML data structures with:

- Pydantic validation
- JAX compatibility for numerical operations
- Immutable data structures
- Bidirectional XML conversion
"""

from .model import (
    FluxomicsData,
    Metadata,
    MetabolicNetworkModel,
    LabelingExperiments,
    Metabolite,
    Reaction,
    AtomTransition,
    AtomMap,
    AtomAddress,
    Measurement,
    MeasurementData,
    Group,
    Constraints,
    NetConstraints,
    ExchangeConstraints,
    MetaboliteSizeConstraints,
    Simulation,
    FluxValue,
    MetaboliteSizeValue,
    Tracers,
    LabelComposition,
    Annotation,
    DictList,
)
from .io import parse_fluxml_file

__version__ = "0.1.0"
__all__ = [
    "FluxomicsData",
    "Metadata",
    "MetabolicNetworkModel",
    "LabelingExperiments",
    "Metabolite",
    "Reaction",
    "AtomTransition",
    "AtomMap",
    "AtomAddress",
    "Measurement",
    "MeasurementData",
    "Group",
    "Constraints",
    "NetConstraints",
    "ExchangeConstraints",
    "MetaboliteSizeConstraints",
    "Simulation",
    "FluxValue",
    "MetaboliteSizeValue",
    "Tracers",
    "LabelComposition",
    "Annotation",
    "DictList",
    "parse_fluxml_file",
]
