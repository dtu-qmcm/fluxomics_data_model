# --- START OF FILE fluxml_gemini.py ---

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
ExperimentDesignWeightType = Annotated[float, Field(ge=0.0, le=1.0)]

# String types with complex patterns - represented as basic strings for now.
# Parsing/validation logic can be added later if needed.
ComposedWeightType = str  # Pattern: \d+(\s*,\s*\d+)? (or .*)

CfgType = str  # Pattern: ([A-Za-z]+|(\s*[CHNOS]#[0-9]+@[A-Za-z0-9]
# +(\s+[CHNOS]#[0-9]+@[A-Za-z0-9]+)*\s*)) (or .*)

AtomCfgType = str  # Pattern: ([A-Za-z]+|(\s*[CHNOS]#[0-9]
# +(\s+[CHNOS]#[0-9]+)*\s*)) (or .*)

LabelCfgType = str  # Pattern: [01xX]+
TimeListType = str  # Comma/space separated list of floats or 'inf'

# Literal types from schema enumerations
GroupScaleTypeEnum = Literal["auto", "one"]
InputTypeEnum = Literal["isotopomer", "cumomer", "emu"]
DatumTypeEnum = Literal["S", "DL", "DR", "DD", "T"]
SimulationTypeEnum = Literal["auto", "explicit", "full"]
SimulationMethodEnum = Literal["auto", "cumomer", "emu"]
FluxValueTypeEnum = Literal["net", "xch"]


# --- Model Definitions ---


class Info(BaseModel):
    """Contains metadata about the FluxML document.

    This class represents the 'fluxml/info' element in the FluxML specification.
    It stores various metadata attributes about the model, such as name, version,
    creation date, and other identifying information.

    Attributes:
        name: Optional name of the FluxML model.
        version: Optional version identifier of the model.
        date: Optional timestamp indicating when the model was created or last modified.
        comment: Optional additional notes or descriptions about the model.
        signature: Optional Base64-encoded binary signature for model authentication.
        modeler: Optional name or identifier of the person or software that created the model.
        strain: Optional biological strain information relevant to the model.
    """  # noqa: E501

    name: Optional[str] = None
    version: Optional[str] = None
    date: Optional[TsType] = None
    comment: Optional[str] = None
    signature: Optional[bytes] = None  # Corresponds to xs:base64Binary
    modeler: Optional[str] = None
    strain: Optional[str] = None


class Annotation(BaseModel):
    """Represents a generic annotation element for metadata.

    This class models annotation elements commonly used in biological data
    exchange formats to associate metadata with biological entities. It corresponds
    to the 'annotation' element found within various FluxML structures
    (e.g., `fluxml/reactionnetwork/metabolitepools/pool/annotation`).

    Attributes:
        annotation_type: The type or category of the annotation (e.g., 'KEGGID', 'Description', 'EC Number').
        value: The actual content or value of the annotation.

    Examples:
        >>> kegg_annotation = Annotation(annotation_type="KEGGID", value="K00001")
        >>> description = Annotation(annotation_type="Description", value="Glucose transporter")
    """  # noqa: E501

    annotation_type: str = Field(
        ...,
        description="The type of annotation "
        "(e.g., 'KEGGID', 'Description', 'EC Number')",
    )
    value: str = Field(..., description="The content of the annotation")


class Pool(BaseModel):
    """Represents a metabolite pool within a metabolic network.

    A metabolite pool defines a chemical species involved in reactions, along
    with properties relevant for isotopic labeling simulations. This model corresponds
    to the 'pool' element in the FluxML specification
    (`fluxml/reactionnetwork/metabolitepools/pool`).

    Attributes:
        id: Unique identifier for the metabolite pool.
        atoms: Number of atoms tracked for isotopic labeling studies (typically carbon atoms). Defaults to 0.
        size: Relative or absolute size/concentration of the pool. Defaults to 1.0.
        cfg: Atom configuration string specifying elemental composition or labeling status. Defaults to '0'.
        annotations: Optional list of annotations providing additional metadata for this pool.
    """  # noqa: E501

    id: str = Field(
        ..., description="Unique identifier for the metabolite pool."
    )
    atoms: AtomType = Field(
        0,
        description="Number of atoms tracked for labeling"
        " (e.g., carbon atoms).",
    )
    size: float = Field(
        1.0, description="Relative or absolute size/concentration of the pool."
    )
    cfg: AtomCfgType = Field(
        "0",
        description="Atom configuration/composition string"
        " (e.g., element counts). Default '0' might imply undefined or "
        "carbon-only tracking.",
    )
    annotations: List[Annotation] = Field(
        default_factory=list,
        description="Optional list of annotations for this pool.",
    )


class MetabolitePools(BaseModel):
    """Container for metabolite pool definitions in the reaction network.

    Corresponds to the 'metabolitepools' element in the FluxML specification
    (`fluxml/reactionnetwork/metabolitepools`).

    Attributes:
        pools: A list containing at least two `Pool` objects defining the metabolites in the network.
    """  # noqa: E501

    pools: Annotated[List[Pool], Field(min_length=2)] = Field(
        ...,
        description="List of metabolite pools defined in the network "
        "(minimum of 2).",
    )


class Variant(BaseModel):
    """Represents an alternative atom mapping or stoichiometry for a reaction component.

    Corresponds to the 'variant' element used within 'reduct' and 'rproduct'
    elements in the FluxML specification
    (`fluxml/reactionnetwork/reaction/(reduct|rproduct)/variant`).

    Attributes:
        cfg: Atom mapping configuration string defining this specific variant.
        ratio: Optional relative ratio or probability of this variant occurring compared to other variants for the same reactant/product.
    """  # noqa: E501

    cfg: CfgType = Field(
        ..., description="Atom mapping configuration string for this variant."
    )
    ratio: Optional[float] = Field(
        None,
        description="Optional relative ratio/probability of this variant"
        " occurring.",
    )


class Reduct(BaseModel):
    """Represents a reactant (educt) participating in a reaction.

    Corresponds to the 'reduct' element within a 'reaction' in the FluxML
    specification (`fluxml/reactionnetwork/reaction/reduct`).

    Attributes:
        pool_id: Identifier of the metabolite pool acting as the reactant (references `Pool.id`).
        cfg: Optional default atom mapping configuration used if no variants are specified. Defines how atoms from this reactant are mapped in the reaction.
        variants: Optional list of `Variant` objects representing alternative atom mappings for this reactant.
    """  # noqa: E501

    pool_id: str = Field(
        ...,
        description="Identifier of the metabolite pool acting as reactant "
        "(references Pool.id).",
    )
    cfg: Optional[CfgType] = Field(
        None,
        description="Optional default atom mapping configuration if no variants"
        " are specified.",
    )
    variants: List[Variant] = Field(
        default_factory=list,
        description="Optional list of alternative atom mappings for this"
        " reactant.",
    )


class RProduct(BaseModel):
    """Represents a product generated by a reaction.

    Corresponds to the 'rproduct' element within a 'reaction' in the FluxML
    specification (`fluxml/reactionnetwork/reaction/rproduct`).

    Attributes:
        pool_id: Identifier of the metabolite pool acting as the product (references `Pool.id`).
        cfg: Optional default atom mapping configuration used if no variants are specified. Defines how atoms map to this product.
        variants: Optional list of `Variant` objects representing alternative atom mappings for this product.
    """  # noqa: E501

    pool_id: str = Field(
        ...,
        description="Identifier of the metabolite pool acting as product "
        "(references Pool.id).",
    )
    cfg: Optional[CfgType] = Field(
        None,
        description="Optional default atom mapping configuration if no"
        " variants are specified.",
    )
    variants: List[Variant] = Field(
        default_factory=list,
        description="Optional list of alternative atom mappings for this "
        " product.",
    )


class Reaction(BaseModel):
    """Represents a metabolic reaction in the network.

    Corresponds to the 'reaction' element in the FluxML specification
    (`fluxml/reactionnetwork/reaction`).

    Attributes:
        id: Unique identifier for the reaction.
        bidirectional: Boolean indicating if the reaction is reversible. Defaults to True.
        annotations: Optional list of annotations providing additional metadata for this reaction.
        reducts: List of `Reduct` objects representing the reactants of the reaction.
        rproducts: List of `RProduct` objects representing the products of the reaction.
    """  # noqa: E501

    id: str = Field(..., description="Unique identifier for the reaction.")
    bidirectional: bool = Field(
        True, description="Indicates if the reaction is reversible."
    )
    annotations: List[Annotation] = Field(
        default_factory=list,
        description="Optional list of annotations for this reaction.",
    )
    # Schema allows 0..unbounded, functionally usually >0 needed.
    reducts: List[Reduct] = Field(
        default_factory=list, description="List of reactants (reducts)."
    )
    # Schema allows 0..unbounded (e.g., sink reactions).
    rproducts: List[RProduct] = Field(
        default_factory=list, description="List of products."
    )


class ReactionNetwork(BaseModel):
    """Defines the complete structure of the metabolic reaction network.

    Corresponds to the 'reactionnetwork' element in the FluxML specification
    (`fluxml/reactionnetwork`).

    Attributes:
        metabolite_pools: A `MetabolitePools` object containing all metabolite pool definitions.
        reactions: A list containing at least one `Reaction` object defining the reactions in the network.
    """  # noqa: E501

    metabolite_pools: MetabolitePools = Field(
        ..., description="Definitions of all metabolite pools in the network."
    )
    # Schema requires 1..unbounded reactions.
    reactions: Annotated[List[Reaction], Field(min_length=1)] = Field(
        ...,
        description="List of reactions defined in the network (minimum of 1).",
    )


# --- Constraints Elements ---


class Textual(BaseModel):
    """Holds a constraint expression represented as a plain string.

    Corresponds to the 'textual' element used within various constraint
    definitions in the FluxML specification (e.g., `fluxml/constraints/net/textual`).

    Attributes:
        value: The string representation of the constraint.
    """  # noqa: E501

    value: str = Field(
        ..., description="The constraint expression as a string."
    )


class NetConstraint(BaseModel):
    """Represents a constraint on net fluxes, typically for stoichiometric balance.

    Corresponds to the 'net' element in the FluxML specification
    (`fluxml/constraints/net`). This constraint usually defines relationships
    that must hold for the net fluxes in the system (e.g., steady-state mass balance).

    Attributes:
        content: The constraint expression, either as a string (from `<textual>`) or potentially a mathematical expression object (e.g., SymPy from `<mml:math>`).
    """  # noqa: E501

    # Holds either string from <textual> or SymPy expr from <mml:math>
    content: MathExpression = Field(
        ..., description="Constraint expression (textual or mathematical)."
    )


class XchConstraint(BaseModel):
    """Represents a constraint specifically on exchange fluxes.

    Corresponds to the 'xch' element in the FluxML specification
    (`fluxml/constraints/xch`). These constraints define relationships or bounds
    involving exchange fluxes.

    Attributes:
        content: The constraint expression, either as a string (from `<textual>`) or potentially a mathematical expression object (e.g., SymPy from `<mml:math>`).
    """  # noqa: E501

    content: MathExpression = Field(
        ..., description="Constraint expression (textual or mathematical)."
    )


class PsizeConstraint(BaseModel):
    """Represents a constraint on metabolite pool sizes.

    Corresponds to the 'psize' element in the FluxML specification
    (`fluxml/constraints/psize`). These constraints are relevant for
    non-stationary models, defining relationships or bounds for pool sizes.

    Attributes:
        content: The constraint expression, either as a string (from `<textual>`) or potentially a mathematical expression object (e.g., SymPy from `<mml:math>`).
    """  # noqa: E501

    content: MathExpression = Field(
        ..., description="Constraint expression (textual or mathematical)."
    )


class Constraints(BaseModel):
    """Container for defining flux balance or pool size constraints.

    Corresponds to the 'constraints' element in the FluxML specification
    (`fluxml/constraints`). These constraints can apply globally or within
    specific configurations.

    Attributes:
        net: Optional constraint on net fluxes (e.g., steady-state balance).
        xch: Optional constraint on exchange fluxes.
        psize: Optional constraint on pool sizes (for non-stationary scenarios).
    """  # noqa: E501

    net: Optional[NetConstraint] = Field(
        None,
        description="Optional net flux constraint"
        " (e.g., stoichiometric balance).",
    )
    xch: Optional[XchConstraint] = Field(
        None, description="Optional exchange flux constraint."
    )
    psize: Optional[PsizeConstraint] = Field(
        None,
        description="Optional pool size constraint"
        " (for non-stationary models).",
    )


# --- Configuration Elements ---


class Sep(BaseModel):
    """Represents a separator marker used in EMU definitions within input labels.

    Corresponds to the empty 'sep' element in the FluxML specification
    (`fluxml/configuration/input/label/sep`). Its presence acts as a marker
    to indicate fragmentation points in EMU definitions.
    """  # noqa: E501

    pass  # Empty element, acts as a marker


class Label(BaseModel):
    """Specifies the isotopic labeling state of an input isotopomer, cumomer, or EMU.

    Corresponds to the 'label' element within an 'input' definition in the
    FluxML specification (`fluxml/configuration/input/label`).

    Attributes:
        cfg: Binary-like string indicating labeled ('1'), unlabeled ('0'), or unknown/variable ('x'/'X') positions.
        purity: Optional string specifying the isotopic purity of the labeled substrate (format may vary).
        cost: Optional cost associated with using this specific labeled substrate.
        value: Optional numerical value or MathML expression representing the fraction (isotopomer), cumomer value, or EMU value associated with this label configuration.
        separators: Optional list of `Sep` markers, used for EMU definitions to indicate fragmentation points. Typically used instead of `value` for EMUs.
    """  # noqa: E501

    cfg: LabelCfgType = Field(
        ...,
        description="Binary string indicating labeled (1), unlabeled (0),"
        " or variable (x/X) positions.",
    )
    purity: Optional[str] = Field(
        None, description="Optional isotopic purity specification."
    )
    cost: Optional[float] = Field(
        None,
        description="Optional cost associated with this labeled substrate.",
    )
    # Represents the value associated with this label configuration.
    # Can be a fraction (isotopomer), cumomer value, or EMU value.
    value: Optional[MathExpression] = Field(
        None,
        description="Optional number (fraction/abundance) "
        "or MathML expression.",
    )
    # Used for EMU definitions.
    separators: Optional[List[Sep]] = Field(
        None, description="Optional list of separators for EMU definitions."
    )

    # TODO: Add validator to ensure either value or separators are relevant?


class Input(BaseModel):
    """Defines a labeled substrate input to the metabolic system.

    Corresponds to the 'input' element in the FluxML specification, found
    within 'configuration' (`fluxml/configuration/input`) or 'mixture'
    (`mixture/input`).

    Attributes:
        pool: Identifier of the metabolite pool being input (references `Pool.id`).
        input_type: The type of labeling representation used ('isotopomer', 'cumomer', or 'emu'). Defaults to 'isotopomer'.
        profile: Optional time profile for substrate input dynamics (comma/space separated times/values or 'inf'). Relevant for non-stationary models.
        labels: List of `Label` objects specifying the different labeling states and their abundances/values for this input pool.
        id: Optional identifier, primarily used when this input is part of a `Mixture` definition.
    """  # noqa: E501

    pool: str = Field(
        ...,
        description="Identifier of the metabolite pool being input"
        " (references Pool.id).",
    )
    input_type: InputTypeEnum = Field(
        "isotopomer", description="Type of labeling representation used."
    )
    profile: Optional[TimeListType] = Field(
        None,
        description="Optional time profile for substrate input "
        "(comma/space separated times/values or 'inf').",
    )
    labels: List[Label] = Field(
        default_factory=list,
        description="List of labeling specifications for this input.",
    )
    # Optional 'id' attribute, primarily for mixture/input
    id: Optional[str] = Field(
        None, description="Optional identifier, mainly used within mixtures."
    )


# --- Measurement Model Elements ---


class ErrorModel(BaseModel):
    """Defines an error model associated with a measurement type or group.

    Corresponds to the 'errormodel' element found within various measurement
    model components in the FluxML specification (e.g.,
    `fluxml/configuration/measurement/model/labelingmeasurement/group/errormodel`).

    Attributes:
        content: The error model definition, either as a string (from `<textual>`) or potentially a mathematical expression object (e.g., SymPy from `<mml:math>`).
    """  # noqa: E501

    content: MathExpression = Field(
        ..., description="Error model definition (textual or mathematical)."
    )


class Group(BaseModel):
    """Defines a group of related labeling measurements, like mass isotopomers of a fragment.

    Corresponds to the 'group' element within 'labelingmeasurement' in the
    FluxML specification
    (`fluxml/configuration/measurement/model/labelingmeasurement/group`).

    Attributes:
        id: Unique identifier for this measurement group (e.g., fragment name).
        times: Optional comma/space separated list of time points for non-stationary measurements within this group.
        scale: Scaling factor type for measurements in this group ('auto' or 'one'). Defaults to 'auto'.
        error_model: Optional `ErrorModel` specific to this group.
        content: Definition of the measured entity (e.g., fragment), either as a string (from `<textual>`) or potentially a mathematical expression object (e.g., SymPy from `<mml:math>`).
    """  # noqa: E501

    id: str = Field(
        ..., description="Unique identifier for this measurement group."
    )
    times: Optional[TimeListType] = Field(
        None,
        description="Optional comma/space separated time points for "
        "non-stationary measurements.",
    )
    scale: GroupScaleTypeEnum = Field(
        "auto",
        description="Scaling factor for measurements in this group "
        "('one' or 'auto').",
    )
    error_model: Optional[ErrorModel] = Field(
        None, description="Optional error model specific to this group."
    )
    # Textual or MathML defining the fragment or measured entity
    content: MathExpression = Field(
        ...,
        description="Definition of the measured entity "
        "(textual or mathematical).",
    )


class LabelingMeasurement(BaseModel):
    """Container for definitions of labeling measurement groups.

    Corresponds to the 'labelingmeasurement' element in the FluxML specification
    (`fluxml/configuration/measurement/model/labelingmeasurement`).

    Attributes:
        groups: List of `Group` objects, each defining a set of related labeling measurements.
    """  # noqa: E501

    groups: List[Group] = Field(
        default_factory=list, description="List of labeling measurement groups."
    )


class NetFlux(BaseModel):
    """Defines a measured net flux value within the measurement model.

    Corresponds to the 'netflux' element within 'fluxmeasurement' in the
    FluxML specification
    (`fluxml/configuration/measurement/model/fluxmeasurement/netflux`).

    Attributes:
        id: Identifier for the measured net flux (can be a reaction ID or a custom name).
        error_model: Optional `ErrorModel` specific to this flux measurement.
        content: Definition or identifier of the measured flux, usually a textual ID (matching `id` or `Reaction.id`) but potentially a mathematical expression.
    """  # noqa: E501

    id: str = Field(
        ...,
        description="Identifier for the measured net flux "
        "(e.g., reaction ID or custom name).",
    )
    error_model: Optional[ErrorModel] = Field(
        None, description="Optional error model for this flux measurement."
    )
    # Usually textual ID, potentially MathML for complex definitions
    content: MathExpression = Field(
        ..., description="Identifier or definition of the measured net flux."
    )


class XchFlux(BaseModel):
    """Defines a measured exchange flux value within the measurement model.

    Corresponds to the 'xchflux' element within 'fluxmeasurement' in the FluxML
    specification (`fluxml/configuration/measurement/model/fluxmeasurement/xchflux`).

    Attributes:
        id: Identifier for the measured exchange flux.
        error_model: Optional `ErrorModel` specific to this flux measurement.
        content: Definition or identifier of the measured flux, usually a textual ID (matching `id` or `Reaction.id`) but potentially a mathematical expression.
    """  # noqa: E501

    id: str = Field(
        ..., description="Identifier for the measured exchange flux."
    )
    error_model: Optional[ErrorModel] = Field(
        None, description="Optional error model for this flux measurement."
    )
    content: MathExpression = Field(
        ...,
        description="Identifier or definition of the measured exchange flux.",
    )


class FluxMeasurement(BaseModel):
    """Container for definitions of flux measurements (net and exchange).

    Corresponds to the 'fluxmeasurement' element in the FluxML specification
    (`fluxml/configuration/measurement/model/fluxmeasurement`).

    Attributes:
        net_fluxes: List of `NetFlux` objects defining measured net fluxes.
        xch_fluxes: List of `XchFlux` objects defining measured exchange fluxes.
    """  # noqa: E501

    net_fluxes: List[NetFlux] = Field(
        default_factory=list,
        description="List of measured net flux definitions.",
    )
    xch_fluxes: List[XchFlux] = Field(
        default_factory=list,
        description="List of measured exchange flux definitions.",
    )


class PoolSize(BaseModel):
    """Defines a measured metabolite pool size within the measurement model.

    Corresponds to the 'poolsize' element within 'poolsizemeasurement' in the
    FluxML specification
    (`fluxml/configuration/measurement/model/poolsizemeasurement/poolsize`).

    Attributes:
        id: Identifier for the measured pool (references `Pool.id`).
        error_model: Optional `ErrorModel` specific to this pool size measurement.
        content: Identifier of the measured pool, typically the `Pool.id`.
    """  # noqa: E501

    id: str = Field(
        ...,
        description="Identifier for the measured pool (references Pool.id).",
    )
    error_model: Optional[ErrorModel] = Field(
        None, description="Optional error model for this pool size measurement."
    )
    content: MathExpression = Field(
        ..., description="Identifier of the measured pool (usually Pool.id)."
    )


class PoolSizeMeasurement(BaseModel):
    """Container for definitions of pool size measurements.

    Corresponds to the 'poolsizemeasurement' element in the FluxML specification
    (`fluxml/configuration/measurement/model/poolsizemeasurement`).

    Attributes:
        pool_sizes: List of `PoolSize` objects defining measured pool sizes.
    """  # noqa: E501

    pool_sizes: List[PoolSize] = Field(
        default_factory=list,
        description="List of measured pool size definitions.",
    )


class NetRatio(BaseModel):
    """Defines a measured ratio between net fluxes within the measurement model.

    Corresponds to the 'netratio' element within 'fluxratios' in the FluxML
    specification (`fluxml/configuration/measurement/model/fluxratios/netratio`).

    Attributes:
        id: Identifier for the measured net flux ratio.
        content: Definition of the ratio, typically a string like 'R1/R2' or a mathematical expression.
    """  # noqa: E501

    id: str = Field(..., description="Identifier for the net flux ratio.")
    content: MathExpression = Field(
        ...,
        description="Definition of the ratio "
        "(e.g., 'R1/R2', textual or mathematical).",
    )


class XchRatio(BaseModel):
    """Defines a measured ratio between exchange fluxes within the measurement model.

    Corresponds to the 'xchratio' element within 'fluxratios' in the FluxML
    specification (`fluxml/configuration/measurement/model/fluxratios/xchratio`).

    Attributes:
        id: Identifier for the measured exchange flux ratio.
        content: Definition of the ratio, typically a string or a mathematical expression.
    """  # noqa: E501

    id: str = Field(..., description="Identifier for the exchange flux ratio.")
    content: MathExpression = Field(
        ..., description="Definition of the ratio (textual or mathematical)."
    )


class FluxRatios(BaseModel):
    """Container for definitions of flux ratio measurements (net and exchange).

    Corresponds to the 'fluxratios' element in the FluxML specification
    (`fluxml/configuration/measurement/model/fluxratios`).

    Attributes:
        net_ratios: List of `NetRatio` objects defining measured net flux ratios.
        xch_ratios: List of `XchRatio` objects defining measured exchange flux ratios.
    """  # noqa: E501

    net_ratios: List[NetRatio] = Field(
        default_factory=list,
        description="List of measured net flux ratio definitions.",
    )
    xch_ratios: List[XchRatio] = Field(
        default_factory=list,
        description="List of measured exchange flux ratio definitions.",
    )


class PoolSizeRatio(BaseModel):
    """Defines a measured ratio between pool sizes within the measurement model.

    Corresponds to the 'poolsizeratio' element within 'poolsizeratios' in the
    FluxML specification
    (`fluxml/configuration/measurement/model/poolsizeratios/poolsizeratio`).

    Attributes:
        id: Identifier for the measured pool size ratio.
        content: Definition of the ratio, typically a string like 'P1/P2' or a mathematical expression.
    """  # noqa: E501

    id: str = Field(..., description="Identifier for the pool size ratio.")
    content: MathExpression = Field(
        ...,
        description="Definition of the ratio "
        "(e.g., 'P1/P2', textual or mathematical).",
    )


class PoolSizeRatios(BaseModel):
    """Container for definitions of pool size ratio measurements.

    Corresponds to the 'poolsizeratios' element in the FluxML specification
    (`fluxml/configuration/measurement/model/poolsizeratios`).

    Attributes:
        pool_size_ratios: List of `PoolSizeRatio` objects defining measured pool size ratios.
    """  # noqa: E501

    pool_size_ratios: List[PoolSizeRatio] = Field(
        default_factory=list,
        description="List of measured pool size ratio definitions.",
    )


class MeasurementModel(BaseModel):
    """Defines all types of measurements expected or included in the experiment(s).

    Corresponds to the 'model' element within 'measurement' in the FluxML
    specification (`fluxml/configuration/measurement/model`). It acts as a
    blueprint for the actual measurement data.

    Attributes:
        labeling_measurement: Optional definitions for labeling measurements (e.g., mass isotopomers).
        flux_measurement: Optional definitions for absolute flux measurements.
        pool_size_measurement: Optional definitions for absolute pool size measurements.
        flux_ratios: Optional definitions for flux ratio measurements.
        pool_size_ratios: Optional definitions for pool size ratio measurements.
    """  # noqa: E501

    labeling_measurement: Optional[LabelingMeasurement] = Field(
        None, description="Optional definition of labeling measurements."
    )
    flux_measurement: Optional[FluxMeasurement] = Field(
        None, description="Optional definition of flux measurements."
    )
    pool_size_measurement: Optional[PoolSizeMeasurement] = Field(
        None, description="Optional definition of pool size measurements."
    )
    flux_ratios: Optional[FluxRatios] = Field(
        None, description="Optional definition of flux ratios."
    )
    pool_size_ratios: Optional[PoolSizeRatios] = Field(
        None, description="Optional definition of pool size ratios."
    )


# --- Measurement Data Elements ---


class Experiment(BaseModel):
    """Contains metadata about the experimental procedure.

    Corresponds to the 'experiment' element within 'dlabel' in the FluxML
    specification (`fluxml/configuration/measurement/data/dlabel/experiment`).

    Attributes:
        operator: Name or identifier of the person who performed the experiment.
        description: Description of the experimental setup or protocol.
    """  # noqa: E501

    operator: str = Field(
        ..., description="Operator who performed the experiment."
    )
    description: str = Field(
        ..., description="Description of the experimental procedure."
    )


class Analytics(BaseModel):
    """Contains metadata about the sample analytics (e.g., MS method).

    Corresponds to the 'analytics' element within 'dlabel' in the FluxML
    specification (`fluxml/configuration/measurement/data/dlabel/analytics`).

    Attributes:
        operator: Name or identifier of the person who performed the analytical measurements.
        description: Description of the analytical method used (e.g., GC-MS, LC-MS details).
    """  # noqa: E501

    operator: str = Field(
        ..., description="Operator who performed the analytics."
    )
    description: str = Field(
        ..., description="Description of the analytical method."
    )


class Analysis(BaseModel):
    """Contains metadata about the data analysis and processing steps.

    Corresponds to the 'analysis' element within 'dlabel' in the FluxML
    specification (`fluxml/configuration/measurement/data/dlabel/analysis`).

    Attributes:
        operator: Name or identifier of the person who performed the data analysis.
        description: Description of the data processing and analysis steps.
    """  # noqa: E501

    operator: str = Field(
        ..., description="Operator who performed the analysis."
    )
    description: str = Field(
        ..., description="Description of the data analysis."
    )


class DLabel(BaseModel):
    """Contains metadata associated specifically with a measurement data set.

    Corresponds to the 'dlabel' element within 'data' in the FluxML
    specification (`fluxml/configuration/measurement/data/dlabel`).

    Attributes:
        strain: Optional biological strain information relevant to this specific dataset.
        date: Optional primary date associated with the dataset (e.g., measurement date).
        start: Optional start timestamp for the experiment or data acquisition.
        finish: Optional finish timestamp for the experiment or data acquisition.
        experiment: Optional detailed metadata about the experimental procedure.
        analytics: Optional detailed metadata about the sample analytics.
        analysis: Optional detailed metadata about the data analysis.
        comment: Optional comments specific to this dataset.
    """  # noqa: E501

    strain: Optional[str] = None
    date: Optional[TsType] = None
    start: Optional[TsType] = None  # Compatibility alias for date?
    finish: Optional[TsType] = None
    experiment: Optional[Experiment] = None
    analytics: Optional[Analytics] = None
    analysis: Optional[Analysis] = None
    comment: Optional[str] = None


class Datum(BaseModel):
    """Represents a single measurement data point.

    Corresponds to the 'datum' element within 'data' in the FluxML
    specification (`fluxml/configuration/measurement/data/datum`).

    Attributes:
        id: Identifier linking this datum to its definition in the `MeasurementModel` (e.g., `Group.id`, `NetFlux.id`).
        stddev: Standard deviation associated with the measured value.
        value: The measured numerical value.
        row: Optional row index, often used for mass isotopomers (M+0, M+1, etc., typically 1-based).
        time: Optional time point associated with this measurement, for non-stationary data.
        weight: Optional weighting factor(s) used in parameter estimation (string format, e.g., '1.0' or '1.0, 0.5').
        pos: Optional atom position information (context-dependent, e.g., for specific atom labeling).
        type: Optional type specifier, primarily for MS/MS fragment data ('S', 'DL', 'DR', 'DD', 'T').
    """  # noqa: E501

    id: str = Field(
        ...,
        description="Identifier linking datum to MeasurementModel definition"
        " (e.g., Group.id, NetFlux.id).",
    )
    stddev: float = Field(
        ..., description="Standard deviation of the measurement."
    )
    value: float = Field(..., description="The measured numerical value.")
    row: Optional[RowType] = Field(
        None,
        description="Optional row index "
        "(e.g., M+0, M+1 for mass isotopomers, 1-based).",
    )
    time: Optional[float] = Field(
        None, description="Optional time point for non-stationary data."
    )
    weight: Optional[ComposedWeightType] = Field(
        None,
        description="Optional weighting factor(s) for parameter estimation.",
    )
    pos: Optional[AtomType] = Field(
        None,
        description="Optional atom position information (context-dependent).",
    )
    type: Optional[DatumTypeEnum] = Field(
        None,
        description="Optional type specifier for MS/MS data "
        "(Singlet, Doublet Left/Right/Double, Triplet).",
    )


class MeasurementData(BaseModel):
    """Container for the actual measurement data points for a configuration.

    Corresponds to the 'data' element within 'measurement' in the FluxML
    specification (`fluxml/configuration/measurement/data`).

    Attributes:
        dlabel: Optional metadata specific to this dataset.
        datum: List of `Datum` objects representing the individual measurement points.
    """  # noqa: E501

    dlabel: Optional[DLabel] = Field(
        None, description="Optional metadata specific to this dataset."
    )
    # Schema allows 0..unbounded datum elements.
    datum: List[Datum] = Field(
        default_factory=list,
        description="List of individual measurement data points.",
    )


# --- Top-Level Measurement Container ---


class MLabel(BaseModel):
    """Contains metadata associated specifically with the overall measurement section.

    Corresponds to the 'mlabel' element within 'measurement' in the FluxML
    specification (`fluxml/configuration/measurement/mlabel`).

    Attributes:
        date: Optional date associated with the measurement definition or collection.
        version: Optional version identifier for the measurement definition.
        comment: Optional comments about the overall measurement setup.
        fluxunit: Optional unit used for flux measurements (e.g., 'mmol/gDW/h').
        poolsizeunit: Optional unit used for pool size measurements (e.g., 'umol/gDW').
        timeunit: Optional unit used for time points (e.g., 'min', 'h').
    """  # noqa: E501

    date: Optional[TsType] = None
    version: Optional[str] = None
    comment: Optional[str] = None
    fluxunit: Optional[str] = None
    poolsizeunit: Optional[str] = None
    timeunit: Optional[str] = None


class Measurement(BaseModel):
    """Encapsulates both the definition (model) and the actual data of measurements.

    Corresponds to the 'measurement' element within 'configuration' in the
    FluxML specification (`fluxml/configuration/measurement`).

    Attributes:
        mlabel: Optional metadata pertaining to the entire measurement section.
        model: A `MeasurementModel` object defining the types of measurements included.
        data: A `MeasurementData` object containing the actual measured values and their metadata.
    """  # noqa: E501

    mlabel: Optional[MLabel] = Field(
        None, description="Optional metadata for the measurement section."
    )
    model: MeasurementModel = Field(
        ..., description="Defines the types of measurements included."
    )
    data: MeasurementData = Field(
        ..., description="Contains the actual measurement data points."
    )


# --- Simulation Elements ---


class FluxValue(BaseModel):
    """Represents a flux value (net or exchange), potentially estimated or fixed, with bounds.

    Corresponds to the 'fluxvalue' element within 'variables' in the FluxML
    specification (`fluxml/configuration/simulation/variables/fluxvalue`).

    Attributes:
        flux: Identifier of the flux (references `Reaction.id` or a custom name, e.g., for combined fluxes).
        type: Specifies whether it's a 'net' or 'xch' (exchange) flux.
        value: The actual numerical flux value (can be fixed, an initial guess, or an estimated result).
        lo: Optional lower bound constraint for the flux value during estimation or analysis.
        hi: Optional upper bound constraint for the flux value during estimation or analysis.
        inc: Optional increment step size (usage depends on context, e.g., parameter scans).
        edweight: Weight associated with this flux for experimental design optimization (0.0 to 1.0). Defaults to 1.0.
    """  # noqa: E501

    flux: str = Field(
        ...,
        description="Identifier of the flux"
        " (references Reaction.id or custom name).",
    )
    type: FluxValueTypeEnum = Field(
        ...,
        description="Specifies whether it's a 'net' or 'xch' (exchange) flux.",
    )
    value: float = Field(
        ..., description="The actual flux value (estimated or fixed)."
    )
    lo: Optional[float] = Field(
        None, description="Optional lower bound for the flux value."
    )
    hi: Optional[float] = Field(
        None, description="Optional upper bound for the flux value."
    )
    inc: Optional[float] = Field(
        None,
        description="Optional increment step "
        "(usage depends on context, e.g., exploration).",
    )
    edweight: ExperimentDesignWeightType = Field(
        1.0,
        description="Weight used in experimental design optimization "
        "(0.0 to 1.0).",
    )


class PoolSizeValue(BaseModel):
    """Represents a pool size value, potentially estimated or fixed, with bounds.

    Corresponds to the 'poolsizevalue' element within 'variables' in the FluxML
    specification (`fluxml/configuration/simulation/variables/poolsizevalue`).

    Attributes:
        pool: Identifier of the pool (references `Pool.id`).
        value: The actual numerical pool size value (can be fixed, an initial guess, or an estimated result).
        lo: Optional lower bound constraint for the pool size during estimation or analysis.
        hi: Optional upper bound constraint for the pool size during estimation or analysis.
        inc: Optional increment step size (usage depends on context, e.g., parameter scans).
        edweight: Weight associated with this pool size for experimental design optimization (0.0 to 1.0). Defaults to 1.0.
    """  # noqa: E501

    pool: str = Field(
        ..., description="Identifier of the pool (references Pool.id)."
    )
    value: float = Field(..., description="The actual pool size value.")
    lo: Optional[float] = Field(
        None, description="Optional lower bound for the pool size."
    )
    hi: Optional[float] = Field(
        None, description="Optional upper bound for the pool size."
    )
    inc: Optional[float] = Field(None, description="Optional increment step.")
    edweight: ExperimentDesignWeightType = Field(
        1.0, description="Weight used in experimental design optimization."
    )


class SimulationVariables(BaseModel):
    """Container for simulation variables, including flux and pool size values.

    Corresponds to the 'variables' element within 'simulation' in the FluxML
    specification (`fluxml/configuration/simulation/variables`). These typically
    represent the free parameters to be estimated or the results of an estimation.

    Attributes:
        flux_values: List of `FluxValue` objects representing the net and exchange fluxes.
        pool_size_values: List of `PoolSizeValue` objects representing the metabolite pool sizes (relevant for non-stationary or specific analyses).
    """  # noqa: E501

    flux_values: List[FluxValue] = Field(
        default_factory=list,
        description="List of flux values (net and exchange).",
    )
    pool_size_values: List[PoolSizeValue] = Field(
        default_factory=list, description="List of pool size values."
    )


class Simulation(BaseModel):
    """Contains settings and results related to a specific simulation or estimation run.

    Corresponds to the 'simulation' element within 'configuration' in the FluxML
    specification (`fluxml/configuration/simulation`).

    Attributes:
        type: Type of simulation balancing ('auto', 'explicit', 'full'). Defaults to 'auto'.
        method: Computational method used for labeling simulation ('auto', 'cumomer', 'emu'). Defaults to 'auto'.
        model: Optional `MeasurementModel` defining the outputs to be simulated or compared against data. The structure is reused from measurement definitions, but the semantic meaning here relates to simulation outputs.
        variables: Optional `SimulationVariables` containing the parameters (fluxes, pool sizes) used for or resulting from the simulation/estimation.
    """  # noqa: E501

    type: SimulationTypeEnum = Field(
        "auto", description="Type of simulation balancing."
    )
    method: SimulationMethodEnum = Field(
        "auto", description="Computational method used for labeling simulation."
    )
    # Schema allows reusing the MeasurementModel structure here,
    # semantic meaning might differ (defines simulated outputs).
    model: Optional[MeasurementModel] = Field(
        None, description="Optional model defining simulated outputs."
    )
    variables: Optional[SimulationVariables] = Field(
        None, description="Optional simulation parameters or results."
    )


# --- Configuration Container ---


class Configuration(BaseModel):
    """Represents a specific experimental setup, simulation scenario, or dataset.

    Corresponds to the 'configuration' element in the FluxML specification
    (`fluxml/configuration`). A FluxML document can contain multiple configurations.

    Attributes:
        name: Unique name identifying this configuration (e.g., 'WT_Glucose', 'Mutant_Glycerol_TimeCourse').
        stationary: Boolean indicating if the system is assumed to be at isotopic steady state for this configuration. Defaults to True.
        time: Optional specific time point associated with this configuration (e.g., for a steady-state snapshot or a specific time in non-stationary data).
        comment: Optional comments describing this specific configuration.
        inputs: List containing at least one `Input` object defining the labeled substrate(s) used.
        constraints: Optional `Constraints` specific to this configuration.
        measurement: Optional `Measurement` object containing the measurement model and data for this configuration.
        simulation: Optional `Simulation` object containing simulation settings and/or results for this configuration.
    """  # noqa: E501

    name: str = Field(..., description="Unique name for this configuration.")
    stationary: bool = Field(
        True,
        description="Indicates if the system is assumed to be at isotopic "
        " steady state.",
    )
    time: Optional[float] = Field(
        None,
        description="Optional specific time point associated with this "
        " configuration.",
    )
    comment: Optional[str] = Field(
        None, description="Optional comments about this configuration."
    )
    # Schema requires at least one input element.
    inputs: Annotated[List[Input], Field(min_length=1)] = Field(
        ..., description="List of labeled substrate inputs (minimum of 1)."
    )
    constraints: Optional[Constraints] = Field(
        None, description="Optional constraints specific to this configuration."
    )
    measurement: Optional[Measurement] = Field(
        None, description="Optional measurement model and data."
    )
    simulation: Optional[Simulation] = Field(
        None, description="Optional simulation settings and/or results."
    )


# --- Root Elements ---


class FluxML(BaseModel):
    """Root element for FluxML documents defining metabolic models and MFA results.

    This class represents the top-level 'fluxml' element and encapsulates the
    entire model definition, experimental configurations, measurements, and
    simulation results.

    Attributes:
        info: Optional metadata about the FluxML document itself.
        reaction_network: The definition of the metabolic network structure, including pools and reactions.
        constraints: Optional global constraints that apply across all configurations.
        configurations: List of `Configuration` objects, each representing a specific experimental setup or simulation scenario.
    """  # noqa: E501

    info: Optional[Info] = Field(
        None, description="Optional metadata about the FluxML document."
    )
    reaction_network: ReactionNetwork = Field(
        ..., description="Definition of the metabolic reaction network."
    )
    constraints: Optional[Constraints] = Field(
        None,
        description="Optional constraints applying to all configurations.",
    )
    configurations: List[Configuration] = Field(
        default_factory=list,
        description="List of experimental/simulation configurations.",
    )

    # TODO: Implement model_validator to check key/keyref integrity
    #  if strict validation is needed.
    # Example: Check if all Pool.id referenced in Reaction, Input,
    #  PoolSizeValue etc. exist.


class Mixture(BaseModel):
    """Root element for defining substrate mixtures, often used in experimental design.

    Corresponds to the top-level 'mixture' element in the FluxML specification.
    This structure is typically used to define or evaluate optimal labeling
    strategies by combining different labeled inputs.

    Attributes:
        inputs: List containing at least one `Input` object, defining the components and their labeling within the mixture. Each input might have an 'id' attribute here.
        objvalue: Optional objective function value associated with this mixture (e.g., cost, information content from experimental design).
    """  # noqa: E501

    # Schema requires at least one input.
    inputs: Annotated[List[Input], Field(min_length=1)] = Field(
        ...,
        description="List of input components defining the mixture "
        " (minimum of 1).",
    )
    objvalue: Optional[float] = Field(
        None,
        description="Optional objective function value associated with this "
        " mixture.",
    )
