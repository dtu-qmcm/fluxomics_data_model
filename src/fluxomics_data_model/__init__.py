"""
FluxML Data Model

A JAX-compatible Python library for handling FluxML data structures with:
- Pydantic validation
- JAX compatibility for numerical operations
- Immutable data structures
- Bidirectional XML conversion
"""

from .models import (
    FluxML,
    Info,
    ReactionNetwork,
    Configuration,
    Pool,
    MetabolitePools,
    Reaction,
    Reduct,
    RProduct,
    Variant,
    Measurement,
    MeasurementData,
    Group,
    Constraints,
    NetConstraints,
    XchConstraints,
    PsizeConstraints,
    Simulation,
    FluxValue,
    PoolSizeValue,
    Input,
    Label,
    Annotation,
)
from .parser import parse_fluxml_file

__version__ = "0.1.0"
__all__ = [
    "FluxML",
    "Info",
    "ReactionNetwork",
    "Configuration",
    "Pool",
    "MetabolitePools",
    "Reaction",
    "Reduct",
    "RProduct",
    "Variant",
    "Measurement",
    "MeasurementData",
    "Group",
    "Constraints",
    "NetConstraints",
    "XchConstraints",
    "PsizeConstraints",
    "Simulation",
    "FluxValue",
    "PoolSizeValue",
    "Input",
    "Label",
    "Annotation",
    "parse_fluxml_file",
]
