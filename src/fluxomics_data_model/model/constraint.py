"""
FluxML constraint definitions.
"""

from typing import Optional
from pydantic import BaseModel, Field
from ..core.common import TextualOrMath


class NetConstraints(BaseModel):
    """
    FluxML net flux constraints.

    Corresponds to fluxml/constraints/net
    """

    expression: TextualOrMath = Field(
        description="Net flux constraint expression"
    )

    class Config:
        frozen = True
        extra = "forbid"


class ExchangeConstraints(BaseModel):
    """
    FluxML exchange flux constraints.

    Corresponds to fluxml/constraints/xch
    """

    expression: TextualOrMath = Field(
        description="Exchange flux constraint expression"
    )

    class Config:
        frozen = True
        extra = "forbid"


class MetaboliteSizeConstraints(BaseModel):
    """
    FluxML metabolite size constraints.

    Corresponds to fluxml/constraints/metabolitesize
    """

    expression: TextualOrMath = Field(
        description="Metabolite size constraint expression"
    )

    class Config:
        frozen = True
        extra = "forbid"


class Constraints(BaseModel):
    """
    FluxML constraints collection.

    Corresponds to fluxml/constraints
    """

    net: Optional[NetConstraints] = Field(
        default=None, description="Net flux constraints"
    )
    xch: Optional[ExchangeConstraints] = Field(
        default=None, description="Exchange flux constraints"
    )
    metabolitesize: Optional[MetaboliteSizeConstraints] = Field(
        default=None, description="Metabolite size constraints"
    )

    class Config:
        frozen = True
        extra = "forbid"
