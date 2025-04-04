from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime
from typing import Optional, List, Any, Type, TypeVar, Dict
from pydantic import ValidationError
from fluxomics_data_model.fluxml_gemini import (
    FluxML,
    Info,
    ReactionNetwork,
    MetabolitePools,
    Pool,
    Reaction,
    Reduct,
    RProduct,
    Variant,
    Constraints,
    NetConstraint,
    XchConstraint,
    Configuration,
    Input,
    Label,
    MathExpression,
    PsizeConstraint,
    Sep,
    TsType,
    AtomType,
    CfgType,
    AtomCfgType,
    LabelCfgType,
    InputTypeEnum,
    Textual,
)

T = TypeVar("T")


def _parse_optional_datetime(text: Optional[str]) -> Optional[TsType]:
    """Safely parses a datetime string."""
    if text is None:
        return None
    try:
        # Adjust format string if necessary based on actual data
        return datetime.strptime(text.strip(), "%Y-%m-%d %H:%M:%S")
    except (ValueError, TypeError):
        return None  # Or raise a more specific error


def _get_attrib(
    element: ET.Element,
    name: str,
    expected_type: Type[T] = str,
    default: Optional[T] = None,
) -> Optional[T]:
    """Safely gets an attribute and converts its type."""
    value = element.get(name)
    if value is None:
        return default
    try:
        return expected_type(value)
    except (ValueError, TypeError):
        print(
            f"Warning: Could not convert attribute '{name}' value '{value}' to {expected_type}. Using default: {default}"
        )
        return default


def _get_bool_attrib(
    element: ET.Element, name: str, default: bool = True
) -> bool:
    """Gets a boolean attribute, handling 'true'/'false'."""
    value = element.get(name)
    if value is None:
        return default
    return value.lower() == "true"


def _get_child_text(element: ET.Element, tag: str) -> Optional[str]:
    """Gets the text content of a direct child element."""
    child = element.find(tag)
    if child is not None and child.text:
        return child.text.strip()
    return None


def _parse_constraints_content(
    element: Optional[ET.Element],
) -> Optional[MathExpression]:
    """Parses the content of a constraint element (net, xch, psize)."""
    if element is None:
        return None
    textual_elem = element.find("textual")
    # mathml_elem = element.find('.//{http://www.w3.org/1998/Math/MathML}math') # Namespace needed for MathML

    if textual_elem is not None and textual_elem.text:
        # For now, store the raw string. Sympy conversion happens later.
        return textual_elem.text.strip()
    # elif mathml_elem is not None:
    #     # Placeholder: Here you'd parse MathML to string or SymPy object
    #     # math_str = ET.tostring(mathml_elem, encoding='unicode')
    #     # return parse_mathml_to_sympy(math_str) # Requires a MathML parser
    #     print(f"Warning: MathML parsing not implemented yet for element {element.tag}")
    #     return None # Or store raw MathML string
    return None


def _parse_label(label_element: ET.Element) -> Label:
    """Parses a <label> element."""
    cfg = _get_attrib(label_element, "cfg", str)
    if cfg is None:
        raise ValueError(
            f"Missing required attribute 'cfg' in label: {ET.tostring(label_element, encoding='unicode')}"
        )

    value_str = label_element.text.strip() if label_element.text else None
    value: Optional[MathExpression] = None
    if value_str:
        try:
            value = float(value_str)
        except ValueError:
            # Keep as string if not a float (might be math expression later)
            value = value_str
            print(
                f"Warning: Label value '{value_str}' is not a float, storing as string."
            )

    # Check for <sep/> children (not in standard.xml example)
    separators = [Sep() for _ in label_element.findall("sep")]

    return Label(
        cfg=cfg,
        purity=_get_attrib(label_element, "purity", str),
        cost=_get_attrib(label_element, "cost", float),
        value=value,
        separators=separators if separators else None,  # Store None if empty
    )


def _parse_input(input_element: ET.Element) -> Input:
    """Parses an <input> element."""
    pool = _get_attrib(input_element, "pool", str)
    if pool is None:
        raise ValueError(
            f"Missing required attribute 'pool' in input: {ET.tostring(input_element, encoding='unicode')}"
        )

    return Input(
        pool=pool,
        input_type=_get_attrib(
            input_element, "type", str, "isotopomer"
        ),  # Default from schema
        profile=_get_attrib(input_element, "profile", str),
        id=_get_attrib(input_element, "id", str),  # Used in <mixture>
        labels=[
            _parse_label(label_elem)
            for label_elem in input_element.findall("label")
        ],
    )


def _parse_variant(variant_element: ET.Element) -> Variant:
    """Parses a <variant> element."""
    cfg = _get_attrib(variant_element, "cfg", str)
    if cfg is None:
        raise ValueError(
            f"Missing required attribute 'cfg' in variant: {ET.tostring(variant_element, encoding='unicode')}"
        )
    return Variant(cfg=cfg, ratio=_get_attrib(variant_element, "ratio", float))


def _parse_reduct(reduct_element: ET.Element) -> Reduct:
    """Parses a <reduct> element."""
    pool_id = _get_attrib(reduct_element, "id", str)
    if pool_id is None:
        raise ValueError(
            f"Missing required attribute 'id' in reduct: {ET.tostring(reduct_element, encoding='unicode')}"
        )

    return Reduct(
        # Map XML 'id' to model 'pool_id'
        pool_id=pool_id,
        cfg=_get_attrib(reduct_element, "cfg", str),
        variants=[
            _parse_variant(var_elem)
            for var_elem in reduct_element.findall("variant")
        ],
    )


def _parse_rproduct(rproduct_element: ET.Element) -> RProduct:
    """Parses an <rproduct> element."""
    pool_id = _get_attrib(rproduct_element, "id", str)
    if pool_id is None:
        raise ValueError(
            f"Missing required attribute 'id' in rproduct: {ET.tostring(rproduct_element, encoding='unicode')}"
        )

    return RProduct(
        # Map XML 'id' to model 'pool_id'
        pool_id=pool_id,
        cfg=_get_attrib(rproduct_element, "cfg", str),
        variants=[
            _parse_variant(var_elem)
            for var_elem in rproduct_element.findall("variant")
        ],
    )


def _parse_reaction(reaction_element: ET.Element) -> Reaction:
    """Parses a <reaction> element."""
    rxn_id = _get_attrib(reaction_element, "id", str)
    if rxn_id is None:
        raise ValueError(
            f"Missing required attribute 'id' in reaction: {ET.tostring(reaction_element, encoding='unicode')}"
        )

    # Handle annotations if needed (not in standard.xml example)
    annotations = []
    # for annot_elem in reaction_element.findall('annotation'):
    #     annotations.append(Annotation(name=annot_elem.get('name'), value=annot_elem.text))

    return Reaction(
        id=rxn_id,
        bidirectional=_get_bool_attrib(
            reaction_element, "bidirectional", True
        ),  # Default from schema
        annotations=annotations,
        reducts=[
            _parse_reduct(reduct_elem)
            for reduct_elem in reaction_element.findall("reduct")
        ],
        rproducts=[
            _parse_rproduct(rproduct_elem)
            for rproduct_elem in reaction_element.findall("rproduct")
        ],
    )


def _parse_pool(pool_element: ET.Element) -> Pool:
    """Parses a <pool> element."""
    pool_id = _get_attrib(pool_element, "id", str)
    if pool_id is None:
        raise ValueError(
            f"Missing required attribute 'id' in pool: {ET.tostring(pool_element, encoding='unicode')}"
        )

    # Handle annotations if needed (not in standard.xml example)
    annotations = []
    # for annot_elem in pool_element.findall('annotation'):
    #     annotations.append(Annotation(name=annot_elem.get('name'), value=annot_elem.text))

    return Pool(
        id=pool_id,
        atoms=_get_attrib(pool_element, "atoms", int, 0),  # Default from schema
        size=_get_attrib(
            pool_element, "size", float, 1.0
        ),  # Default from schema
        cfg=_get_attrib(pool_element, "cfg", str, "0"),  # Default from schema
        annotations=annotations,
    )


def _parse_metabolite_pools(metpools_element: ET.Element) -> MetabolitePools:
    """Parses the <metabolitepools> element."""
    pools = [
        _parse_pool(pool_elem) for pool_elem in metpools_element.findall("pool")
    ]
    return MetabolitePools(pools=pools)


def _parse_reaction_network(network_element: ET.Element) -> ReactionNetwork:
    """Parses the <reactionnetwork> element."""
    metpools_elem = network_element.find("metabolitepools")
    if metpools_elem is None:
        raise ValueError(
            "Missing required <metabolitepools> element in <reactionnetwork>"
        )

    metabolite_pools = _parse_metabolite_pools(metpools_elem)
    reactions = [
        _parse_reaction(rxn_elem)
        for rxn_elem in network_element.findall("reaction")
    ]

    return ReactionNetwork(
        metabolite_pools=metabolite_pools, reactions=reactions
    )


def _parse_constraints(constraints_element: ET.Element) -> Constraints:
    """Parses a <constraints> element (global or local)."""
    # --- Adapt tag names based on standard.xml example ---
    # Schema uses <net>, <xch>, <psize>
    # Example uses <netto>, <exchange>
    net_elem = constraints_element.find("netto")  # Use 'netto' from example
    xch_elem = constraints_element.find(
        "exchange"
    )  # Use 'exchange' from example
    psize_elem = constraints_element.find("psize")
    # ---

    net_content = _parse_constraints_content(net_elem)
    xch_content = _parse_constraints_content(xch_elem)
    psize_content = _parse_constraints_content(psize_elem)

    return Constraints(
        # Map parsed content to schema-defined model fields 'net', 'xch'
        net=NetConstraint(content=net_content) if net_content else None,
        xch=XchConstraint(content=xch_content) if xch_content else None,
        psize=PsizeConstraint(content=psize_content) if psize_content else None,
    )


def _parse_configuration(config_element: ET.Element) -> Configuration:
    """Parses a <configuration> element."""
    name = _get_attrib(config_element, "name", str)
    if name is None:
        raise ValueError(
            f"Missing required attribute 'name' in configuration: {ET.tostring(config_element, encoding='unicode')}"
        )

    # Parse local constraints if present
    local_constraints_elem = config_element.find("constraints")
    local_constraints = (
        _parse_constraints(local_constraints_elem)
        if local_constraints_elem is not None
        else None
    )

    # Parse measurement and simulation (not in standard.xml example)
    measurement = None  # Placeholder
    simulation = None  # Placeholder

    return Configuration(
        name=name,
        stationary=_get_bool_attrib(
            config_element, "stationary", True
        ),  # Default from schema
        time=_get_attrib(config_element, "time", float),
        comment=_get_child_text(config_element, "comment"),
        inputs=[
            _parse_input(input_elem)
            for input_elem in config_element.findall("input")
        ],
        constraints=local_constraints,
        measurement=measurement,  # Add parsing if needed
        simulation=simulation,  # Add parsing if needed
    )


def _parse_info(info_element: ET.Element) -> Info:
    """Parses the <info> element."""
    return Info(
        name=_get_child_text(info_element, "name"),
        version=_get_child_text(info_element, "version"),
        date=_parse_optional_datetime(_get_child_text(info_element, "date")),
        comment=_get_child_text(info_element, "comment"),
        signature=None,  # signature uses base64Binary, needs specific parsing if present
        modeler=_get_child_text(info_element, "modeler"),
        strain=_get_child_text(info_element, "strain"),
    )


# --- Main Parsing Function ---


def parse_fluxml_xml(xml_file_path: str) -> FluxML:
    """
    Parses a FluxML XML file into Pydantic models.

    Args:
        xml_file_path: Path to the FluxML XML file.

    Returns:
        An instance of the FluxML Pydantic model.

    Raises:
        FileNotFoundError: If the XML file does not exist.
        ET.ParseError: If the XML is malformed.
        ValueError: If required elements or attributes are missing or invalid.
        ValidationError: If the parsed data doesn't match Pydantic model constraints.
    """
    try:
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
    except FileNotFoundError:
        raise FileNotFoundError(f"XML file not found: {xml_file_path}")
    except ET.ParseError as e:
        raise ET.ParseError(f"Error parsing XML file {xml_file_path}: {e}")

    if root.tag != "fluxml":
        # Handle potential namespace if present, e.g. {http://www.13cflux.net/fluxml}fluxml
        if not root.tag.endswith("fluxml"):
            raise ValueError(
                f"Root element is not <fluxml>, found <{root.tag}>"
            )

    # Parse Info
    info_elem = root.find("info")
    info = _parse_info(info_elem) if info_elem is not None else None

    # Parse Reaction Network (Required)
    network_elem = root.find("reactionnetwork")
    if network_elem is None:
        raise ValueError(
            "Missing required <reactionnetwork> element in <fluxml>"
        )
    reaction_network = _parse_reaction_network(network_elem)

    # Parse Global Constraints
    global_constraints_elem = root.find("constraints")
    global_constraints = (
        _parse_constraints(global_constraints_elem)
        if global_constraints_elem is not None
        else None
    )

    # Parse Configurations
    configurations = [
        _parse_configuration(conf_elem)
        for conf_elem in root.findall("configuration")
    ]

    # Validate and instantiate the final FluxML model
    try:
        fluxml_model = FluxML(
            info=info,
            reaction_network=reaction_network,
            constraints=global_constraints,
            configurations=configurations,
        )
        return fluxml_model
    except ValidationError as e:
        print(f"Pydantic validation error constructing FluxML model: {e}")
        raise  # Re-raise the validation error


# --- Example Usage ---
if __name__ == "__main__":
    here = Path(__file__).parent
    xml_file = here.parent.parent / "data" / "standard.xml"
    try:
        print(f"Parsing {xml_file}...")
        fluxml_data = parse_fluxml_xml(xml_file)
        print("Parsing successful!")

        # Optionally print some parsed data
        print(
            f"\nModel Name: {fluxml_data.info.name if fluxml_data.info else 'N/A'}"
        )
        print(
            f"Number of Pools: {len(fluxml_data.reaction_network.metabolite_pools.pools)}"
        )
        print(
            f"Number of Reactions: {len(fluxml_data.reaction_network.reactions)}"
        )
        print(f"Number of Configurations: {len(fluxml_data.configurations)}")

        if fluxml_data.configurations:
            config = fluxml_data.configurations[0]
            print(f"\nFirst Configuration Name: {config.name}")
            print(f"  Number of Inputs: {len(config.inputs)}")
            if config.inputs:
                print(f"    First Input Pool: {config.inputs[0].pool}")
                print(
                    f"    Number of Labels in First Input: {len(config.inputs[0].labels)}"
                )
            if config.constraints and config.constraints.net:
                print(
                    f"  Local Net Constraints: {config.constraints.net.content[:100]}..."
                )  # Print start of constraints
            elif fluxml_data.constraints and fluxml_data.constraints.net:
                print(
                    f"  Global Net Constraints: {fluxml_data.constraints.net.content[:100]}..."
                )

        # # You can also dump to JSON (requires Pydantic v2)
        # try:
        #     print("\nDumping model to JSON (first few levels)...")
        #     # Be careful, this can be very large for complex models
        #     print(fluxml_data.model_dump_json(indent=2))
        # except AttributeError:
        #      print("(model_dump_json requires Pydantic v2)")
        # except Exception as json_e:
        #      print(f"Error dumping JSON: {json_e}")

    except (FileNotFoundError, ET.ParseError, ValueError, ValidationError) as e:
        print(f"\nError parsing {xml_file}: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        import traceback

        traceback.print_exc()
