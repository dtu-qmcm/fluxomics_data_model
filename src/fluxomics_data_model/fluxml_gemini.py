from datetime import datetime
from typing import Annotated, Any, List, Literal, Optional, Union

from pydantic import BaseModel, Field

# --- Type Aliases and Constrained Types ---

# Placeholder for SymPy expressions. Actual integration would require SymPy.
# For now, we allow string representations or Any
# (to hold a sympy object later).
MathExpression = Union[str, Any]

# Datetime format from schema: '%Y-%m-%d %H:%M:%S'
# Pydantic can validate datetime objects, conversion handled separately
# if needed.
TsType = datetime

# Integer constrained by schema
AtomType = Annotated[int, Field(ge=0, le=1024)]
WeightType = Annotated[int, Field(ge=0, le=16)]
RowType = Annotated[int, Field(ge=1, le=256)]
EdWeightType = Annotated[float, Field(ge=0.0, le=1.0)]

# String types with complex patterns - represented as basic strings for now.
# Parsing/validation logic can be added later if needed.
ComposedWeightType = str  # Pattern: \d+(\s*,\s*\d+)? (or .*)

CfgType = str  # Pattern: ([A-Za-z]+|(\s*[CHNOS]#[0-9]+@[A-Za-z0-9]
# +(\s+[CHNOS]#[0-9]+@[A-Za-z0-9]+)*\s*)) (or .*)

AtomCfgType = str  # Pattern: ([A-Za-z]+|(\s*[CHNOS]#[0-9]
# +(\s+[CHNOS]#[0-9]+)*\s*)) (or .*)

LabelCfgType = str  # Pattern: [01xX]+
TimesType = str  # Comma/space separated list of floats or 'inf'

# Literal types from schema enumerations
GroupScaleType = Literal["auto", "one"]
InputTypeEnum = Literal["isotopomer", "cumomer", "emu"]
DatumTypeEnum = Literal["S", "DL", "DR", "DD", "T"]
SimulationTypeEnum = Literal["auto", "explicit", "full"]
SimulationMethodEnum = Literal["auto", "cumomer", "emu"]
FluxValueTypeEnum = Literal["net", "xch"]


# --- Model Definitions ---


class Info(BaseModel):
    """
    Contains metadata about the FluxML document.
    Corresponds to fluxml/info
    """

    name: Optional[str] = None
    version: Optional[str] = None
    date: Optional[TsType] = None
    comment: Optional[str] = None
    signature: Optional[bytes] = None  # Corresponds to xs:base64Binary
    modeler: Optional[str] = None
    strain: Optional[str] = None


class Annotation(BaseModel):
    """
    Generic annotation element.
    Corresponds to fluxml/.../annotation
    """

    name: str  # The type of annotation (e.g., 'KEGGID', 'Description')
    value: str  # The content of the annotation


class Pool(BaseModel):
    """
    Represents a metabolite pool.
    Corresponds to fluxml/reactionnetwork/metabolitepools/pool
    """

    id: str = Field(
        ..., description="Unique identifier for the metabolite pool."
    )
    atoms: AtomType = Field(
        0,
        description="Number of atoms tracked for labeling (e.g. carbon atoms).",
    )
    size: float = Field(
        1.0, description="Relative or absolute size/concentration of the pool."
    )
    cfg: AtomCfgType = Field(
        "0",
        description="Atom configuration/composition string e.g. element counts"
        ". Default '0' might imply undefined or carbon-only tracking.",
    )
    annotations: List[Annotation] = Field(default_factory=list)


class MetabolitePools(BaseModel):
    """
    Container for metabolite pool definitions.
    Corresponds to fluxml/reactionnetwork/metabolitepools
    """

    pools: Annotated[List[Pool], Field(min_length=2)]


class Variant(BaseModel):
    """
    Represents alternative atom mappings or stoichiometries within a reaction.
    Corresponds to fluxml/reactionnetwork/reaction/(reduct|rproduct)/variant
    """

    cfg: CfgType = Field(
        ..., description="Atom mapping configuration string for this variant."
    )
    ratio: Optional[float] = Field(
        None,
        description="Relative ratio/probability of this variant occurring.",
    )


class Reduct(BaseModel):
    """
    Represents a reactant (educt) in a reaction.
    Corresponds to fluxml/reactionnetwork/reaction/reduct
    """

    pool_id: str = Field(
        ...,
        description="Identifier of the metabolite pool acting as reactant"
        " (references Pool.id).",
    )
    cfg: Optional[CfgType] = Field(
        None,
        description="Default atom mapping configuration if no variants are"
        " specified.",
    )
    variants: List[Variant] = Field(
        default_factory=list,
        description="Alternative atom mappings for this reactant.",
    )


class RProduct(BaseModel):
    """
    Represents a product in a reaction.
    Corresponds to fluxml/reactionnetwork/reaction/rproduct
    """

    pool_id: str = Field(
        ...,
        description="Identifier of the metabolite pool acting as product"
        " (references Pool.id).",
    )
    cfg: Optional[CfgType] = Field(
        None,
        description="Default atom mapping configuration if no variants are "
        "specified.",
    )
    variants: List[Variant] = Field(
        default_factory=list,
        description="Alternative atom mappings for this product.",
    )


class Reaction(BaseModel):
    """
    Represents a metabolic reaction.
    Corresponds to fluxml/reactionnetwork/reaction
    """

    id: str = Field(..., description="Unique identifier for the reaction.")
    bidirectional: bool = Field(
        True, description="Indicates if the reaction is reversible."
    )
    annotations: List[Annotation] = Field(default_factory=list)
    # Schema allows 0..unbounded, but functionally >0 reducts usually needed.
    # Sticking to schema.
    reducts: List[Reduct] = Field(default_factory=list)
    # Schema allows 0..unbounded products (e.g., sink reactions).
    rproducts: List[RProduct] = Field(default_factory=list)


class ReactionNetwork(BaseModel):
    """
    Defines the metabolic network structure.
    Corresponds to fluxml/reactionnetwork
    """

    metabolite_pools: MetabolitePools
    # Schema allows 1..unbounded reactions.
    reactions: Annotated[List[Reaction], Field(min_length=1)]


# --- Constraints Elements ---


class Textual(BaseModel):
    """
    Holds a string-based constraint expression.
    Corresponds to fluxml/constraints/.../textual
    """

    value: str


# Classes representing elements that can contain either textual
# or MathML constraints
class NetConstraint(BaseModel):
    """Corresponds to fluxml/constraints/net"""

    # Holds either string from <textual>
    # or SymPy expr from <mml:math>
    content: MathExpression


class XchConstraint(BaseModel):
    """Corresponds to fluxml/constraints/xch"""

    content: MathExpression


class PsizeConstraint(BaseModel):
    """Corresponds to fluxml/constraints/psize"""

    content: MathExpression


class Constraints(BaseModel):
    """
    Container for flux balance or pool size constraints.
    Corresponds to fluxml/constraints
    """

    net: Optional[NetConstraint] = (
        None  # Net flux constraints (usually stoichiometric balance)
    )
    xch: Optional[XchConstraint] = None  # Exchange flux constraints
    psize: Optional[PsizeConstraint] = (
        None  # Pool size constraints (for non-stationary)
    )


# --- Configuration Elements ---


class Sep(BaseModel):
    """
    Represents a separator marker used in EMU definitions within input labels.
    Corresponds to fluxml/configuration/input/label/sep
    """

    pass  # Empty element, acts as a marker


class Label(BaseModel):
    """
    Specifies the isotopic labeling state of a specific isotopomer/cumomer/EMU.
    Corresponds to fluxml/configuration/input/label
    """

    cfg: LabelCfgType = Field(
        ...,
        description="Binary string indicating labeled (1) or unlabeled (0)"
        " positions, 'x' for unknown/variable.",
    )
    purity: Optional[str] = Field(
        None, description="Isotopic purity specification (format might vary)."
    )
    cost: Optional[float] = Field(
        None, description="Cost associated with this labeled substrate."
    )
    # Represents the value associated with this label configuration.
    # Can be a fraction (isotopomer), cumomer value, or EMU value.
    # The choice between value/separators might need custom logic
    #  during parsing.
    value: Optional[MathExpression] = Field(
        None,
        description="Number (fraction/abundance) or MathML expression.",
    )
    separators: Optional[List[Sep]] = Field(
        None, description="Used for EMU definitions, indicates fragmentation."
    )

    # TODO: Add validator to ensure either value or separators are
    #  relevant based on context?


class Input(BaseModel):
    """
    Defines a labeled substrate input to the system.
    Corresponds to fluxml/configuration/input and mixture/input
    """

    pool: str = Field(
        ...,
        description="Identifier of the metabolite pool being input "
        "(references Pool.id).",
    )
    input_type: InputTypeEnum = Field(
        "isotopomer", description="Type of labeling representation used."
    )
    profile: Optional[TimesType] = Field(
        None,
        description="Time profile for substrate input"
        " (comma/space separated times/values).",
    )
    labels: List[Label] = Field(
        default_factory=list,
        description="List of labeling specifications for this input.",
    )
    # Optional 'id' attribute, primarily for mixture/input
    id: Optional[str] = None


# --- Measurement Model Elements ---


class ErrorModel(BaseModel):
    """
    Defines an error model for measurements.
    Corresponds to fluxml/configuration/measurement/model/.../errormodel
    """

    content: MathExpression  # Textual or MathML description of the error model


class Group(BaseModel):
    """
    Defines a group of labeling measurements (e.g., mass isotopomers of a fragment).
    Corresponds to fluxml/configuration/measurement/model/labelingmeasurement/group
    """

    id: str = Field(
        ..., description="Unique identifier for this measurement group."
    )
    times: Optional[TimesType] = Field(
        None,
        description="Comma/space separated time points for non-stationary"
        " measurements.",
    )
    scale: GroupScaleType = Field(
        "auto",
        description="Scaling factor for the measurements in this group"
        " ('one' or 'auto').",
    )
    error_model: Optional[ErrorModel] = None

    # Textual or MathML defining the fragment or measured entity
    content: MathExpression


class LabelingMeasurement(BaseModel):
    """
    Container for labeling measurement group definitions.
    Corresponds to fluxml/configuration/measurement/model/labelingmeasurement
    """

    groups: List[Group] = Field(default_factory=list)


class NetFlux(BaseModel):
    """
    Definition of a measured net flux.
    Corresponds to fluxml/configuration/measurement/model/fluxmeasurement/netflux
    """

    id: str = Field(
        ...,
        description="Identifier for the measured net flux"
        "(can be reaction ID or custom name).",
    )
    error_model: Optional[ErrorModel] = None

    # Usually textual ID, potentially MathML for complex definitions
    content: MathExpression


class XchFlux(BaseModel):
    """
    Definition of a measured exchange flux.
    Corresponds to fluxml/configuration/measurement/model/fluxmeasurement/xchflux.
    """  # noqa: E501

    id: str = Field(
        ..., description="Identifier for the measured exchange flux."
    )
    error_model: Optional[ErrorModel] = None
    content: MathExpression


class FluxMeasurement(BaseModel):
    """
    Container for flux measurement definitions.
    Corresponds to fluxml/configuration/measurement/model/fluxmeasurement
    """

    net_fluxes: List[NetFlux] = Field(default_factory=list)
    xch_fluxes: List[XchFlux] = Field(default_factory=list)


class PoolSize(BaseModel):
    """
    Definition of a measured pool size.
    Corresponds to fluxml/configuration/measurement/model/poolsizemeasurement/poolsize.
    """  # noqa: E501

    id: str = Field(
        ...,
        description="Identifier for the measured pool (references Pool.id).",
    )
    error_model: Optional[ErrorModel] = None
    content: MathExpression  # Usually textual ID (Pool.id)


class PoolSizeMeasurement(BaseModel):
    """
    Container for pool size measurement definitions.
    Corresponds to fluxml/configuration/measurement/model/poolsizemeasurement
    """

    pool_sizes: List[PoolSize] = Field(default_factory=list)


class NetRatio(BaseModel):
    """
    Definition of a measured net flux ratio.
    Corresponds to fluxml/configuration/measurement/model/fluxratios/netratio
    """

    id: str = Field(..., description="Identifier for the net flux ratio.")
    content: (
        MathExpression  # Textual or MathML defining the ratio (e.g., 'R1/R2')
    )


class XchRatio(BaseModel):
    """
    Definition of a measured exchange flux ratio.
    Corresponds to fluxml/configuration/measurement/model/fluxratios/xchratio
    """

    id: str = Field(..., description="Identifier for the exchange flux ratio.")
    content: MathExpression


class FluxRatios(BaseModel):
    """
    Container for flux ratio definitions.
    Corresponds to fluxml/configuration/measurement/model/fluxratios
    """

    net_ratios: List[NetRatio] = Field(default_factory=list)
    xch_ratios: List[XchRatio] = Field(default_factory=list)


class PoolSizeRatio(BaseModel):
    """
    Definition of a measured pool size ratio.
    Corresponds to fluxml/configuration/measurement/model/poolsizeratios/poolsizeratio
    """  # noqa: E501

    id: str = Field(..., description="Identifier for the pool size ratio.")
    content: (
        MathExpression  # Textual or MathML defining the ratio (e.g., 'P1/P2')
    )


class PoolSizeRatios(BaseModel):
    """
    Container for pool size ratio definitions.
    Corresponds to fluxml/configuration/measurement/model/poolsizeratios
    """

    pool_size_ratios: List[PoolSizeRatio] = Field(default_factory=list)


class MeasurementModel(BaseModel):
    """
    Defines all types of measurements included in the experiment.
    Corresponds to fluxml/configuration/measurement/model
    """

    labeling_measurement: Optional[LabelingMeasurement] = None
    flux_measurement: Optional[FluxMeasurement] = None
    pool_size_measurement: Optional[PoolSizeMeasurement] = None
    flux_ratios: Optional[FluxRatios] = None
    pool_size_ratios: Optional[PoolSizeRatios] = None


# --- Measurement Data Elements ---


class Experiment(BaseModel):
    """Metadata about the experimental procedure."""

    operator: str
    description: str


class Analytics(BaseModel):
    """Metadata about the sample analytics (e.g., MS method)."""

    operator: str
    description: str


class Analysis(BaseModel):
    """Metadata about the data analysis/processing."""

    operator: str
    description: str


class DLabel(BaseModel):
    """
    Metadata associated with the measurement data set.
    Corresponds to fluxml/configuration/measurement/data/dlabel
    """

    strain: Optional[str] = None
    date: Optional[TsType] = None
    start: Optional[TsType] = None  # Compatibility alias for date?
    finish: Optional[TsType] = None
    experiment: Optional[Experiment] = None
    analytics: Optional[Analytics] = None
    analysis: Optional[Analysis] = None
    comment: Optional[str] = None


class Datum(BaseModel):
    """
    A single measurement data point.
    Corresponds to fluxml/configuration/measurement/data/datum
    """

    id: str = Field(
        ...,
        description="Identifier linking this datum to a definition in"
        " MeasurementModel (e.g., Group.id, NetFlux.id).",
    )
    stddev: float = Field(
        ..., description="Standard deviation of the measurement."
    )
    value: float  # The measured numerical value (extracted from mixed content)
    row: Optional[RowType] = Field(
        None,
        description="Row index, often M+0, M+1 etc. for mass isotopomers"
        " (1-based).",
    )
    time: Optional[float] = Field(
        None, description="Time point for non-stationary data."
    )
    weight: Optional[ComposedWeightType] = Field(
        None, description="Weighting factor(s) for parameter estimation."
    )
    pos: Optional[AtomType] = Field(
        None, description="Atom position information (context-dependent)."
    )
    type: Optional[DatumTypeEnum] = Field(
        None,
        description="Type specifier for MS/MS data "
        "(Singlet, Doublet Left/Right/Double, Triplet).",
    )


class MeasurementData(BaseModel):
    """
    Container for the actual measurement data points.
    Corresponds to fluxml/configuration/measurement/data
    """

    dlabel: Optional[DLabel] = None
    # Schema allows 0..unbounded datum elements.
    datum: List[Datum] = Field(default_factory=list)


# --- Top-Level Measurement Container ---


class MLabel(BaseModel):
    """
    Metadata associated specifically with the measurement section.
    Corresponds to fluxml/configuration/measurement/mlabel
    """

    date: Optional[TsType] = None
    version: Optional[str] = None
    comment: Optional[str] = None
    fluxunit: Optional[str] = None
    poolsizeunit: Optional[str] = None
    timeunit: Optional[str] = None


class Measurement(BaseModel):
    """
    Encapsulates the definition (model) and the actual data of measurements.
    Corresponds to fluxml/configuration/measurement
    """

    mlabel: Optional[MLabel] = None
    model: MeasurementModel
    data: MeasurementData


# --- Simulation Elements ---


class FluxValue(BaseModel):
    """
    Represents a flux value, potentially estimated or fixed, with bounds.
    Corresponds to fluxml/configuration/simulation/variables/fluxvalue
    """

    flux: str = Field(
        ...,
        description="Identifier of the flux "
        "(references Reaction.id or custom name).",
    )
    type: FluxValueTypeEnum = Field(
        ...,
        description="Specifies whether it's a 'net' or 'xch' (exchange) flux.",
    )
    value: float  # The actual flux value (estimated or fixed)
    lo: Optional[float] = Field(
        None, description="Lower bound for the flux value."
    )
    hi: Optional[float] = Field(
        None, description="Upper bound for the flux value."
    )
    inc: Optional[float] = Field(
        None,
        description="Increment step"
        " (usage depends on context, e.g., exploration).",
    )
    edweight: EdWeightType = Field(
        1.0,
        description="Weight used in experimental design optimization "
        "(0.0 to 1.0).",
    )


class PoolSizeValue(BaseModel):
    """
    Represents a pool size value, potentially estimated or fixed, with bounds.
    Corresponds to fluxml/configuration/simulation/variables/poolsizevalue
    """

    pool: str = Field(
        ..., description="Identifier of the pool (references Pool.id)."
    )
    value: float  # The actual pool size value
    lo: Optional[float] = Field(
        None, description="Lower bound for the pool size."
    )
    hi: Optional[float] = Field(
        None, description="Upper bound for the pool size."
    )
    inc: Optional[float] = Field(None, description="Increment step.")
    edweight: EdWeightType = Field(
        1.0, description="Weight used in experimental design optimization."
    )


class SimulationVariables(BaseModel):
    """
    Container for simulation variables (fluxes and pool sizes).
    Corresponds to fluxml/configuration/simulation/variables
    """

    flux_values: List[FluxValue] = Field(default_factory=list)
    pool_size_values: List[PoolSizeValue] = Field(default_factory=list)


class Simulation(BaseModel):
    """
    Contains simulation settings and results (variables).
    Corresponds to fluxml/configuration/simulation
    """

    type: SimulationTypeEnum = Field(
        "auto", description="Type of simulation balancing."
    )
    method: SimulationMethodEnum = Field(
        "auto", description="Computational method used for labeling simulation."
    )
    # Schema allows reusing the MeasurementModel structure here,
    # semantic meaning might differ.
    model: Optional[MeasurementModel] = None
    variables: Optional[SimulationVariables] = None


# --- Configuration Container ---


class Configuration(BaseModel):
    """
    Represents a specific experimental setup or simulation scenario.
    Corresponds to fluxml/configuration
    """

    name: str = Field(..., description="Unique name for this configuration.")
    stationary: bool = Field(
        True,
        description="Indicates if the system is assumed to be at "
        "isotopic steady state.",
    )
    time: Optional[float] = Field(
        None,
        description="Specific time point associated with this configuration"
        " (e.g., for steady-state snapshot).",
    )
    comment: Optional[str] = None
    # Schema requires at least one input element.
    inputs: Annotated[List[Input], Field(min_length=1)]
    constraints: Optional[Constraints] = None
    measurement: Optional[Measurement] = None
    simulation: Optional[Simulation] = None


# --- Root Elements ---


class FluxML(BaseModel):
    """
    Root element for FluxML documents defining metabolic models and MFA results.
    Corresponds to fluxml
    """

    info: Optional[Info] = None
    reaction_network: ReactionNetwork
    constraints: Optional[Constraints] = (
        None  # Global constraints applying to all configurations
    )
    configurations: List[Configuration] = Field(default_factory=list)

    # TODO: Implement model_validator to check key/keyref integrity
    #  if strict validation is needed.
    # Example: Check if all Pool.id referenced in Reaction, Input,
    #  PoolSizeValue etc. exist.


class Mixture(BaseModel):
    """
    Root element for defining substrate mixtures, often used in experimental
     design. Corresponds to mixture
    """

    # Schema requires at least one input.
    inputs: Annotated[List[Input], Field(min_length=1)]
    objvalue: Optional[float] = Field(
        None,
        description="Objective function value associated with this mixture"
        " (e.g., cost, information content).",
    )


# --- Example Usage (Optional) ---
