"""
FluxML constraint definitions.
"""

from typing import Optional
from pydantic import BaseModel, Field
from .common import TextualOrMath


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


class XchConstraints(BaseModel):
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


class PsizeConstraints(BaseModel):
    """
    FluxML pool size constraints.

    Corresponds to fluxml/constraints/psize
    """

    expression: TextualOrMath = Field(
        description="Pool size constraint expression"
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
    xch: Optional[XchConstraints] = Field(
        default=None, description="Exchange flux constraints"
    )
    psize: Optional[PsizeConstraints] = Field(
        default=None, description="Pool size constraints"
    )

    class Config:
        frozen = True
        extra = "forbid"
