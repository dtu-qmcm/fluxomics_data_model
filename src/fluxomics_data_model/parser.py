"""
FluxML XML parser for reading FluxML files.
"""

import xml.etree.ElementTree as ET
from typing import Optional
from datetime import datetime

from .models import (
    FluxML,
    Info,
    ReactionNetwork,
    Metabolite,
    Metabolites,
    Reaction,
    Reduct,
    RProduct,
    Variant,
    Experiments,
    Tracers,
    Label,
    Measurement,
    MeasurementData,
    MeasurementModel,
    Datum,
    LabelingMeasurement,
    Group,
    FluxMeasurement,
    NetFlux,
    XchFlux,
    MetaboliteSizeMeasurement,
    MetaboliteSize,
    Simulation,
    Constraints,
    NetConstraints,
    XchConstraints,
    MetaboliteSizeConstraints,
    Annotation,
    TextualOrMath,
    ErrorModel,
)


class FluxMLParser:
    """Parser for FluxML XML files."""

    def __init__(self):
        self.namespaces = {
            "fluxml": "http://www.13cflux.net/fluxml",
            "mml": "http://www.w3.org/1998/Math/MathML",
        }

    def parse_file(self, file_path: str) -> FluxML:
        """Parse a FluxML file and return a FluxML object."""
        try:
            # Try to parse with different encodings
            encodings = ["utf-8", "latin-1", "iso-8859-1", "cp1252"]

            for encoding in encodings:
                try:
                    with open(file_path, "r", encoding=encoding) as f:
                        content = f.read()

                    # Parse the XML content
                    root = ET.fromstring(content)
                    return self._parse_fluxml(root)

                except UnicodeDecodeError:
                    continue
                except ET.ParseError:
                    # If XML parsing fails, try next encoding
                    continue

            # If all encodings failed, try binary mode with ET.parse
            try:
                tree = ET.parse(file_path)
                root = tree.getroot()
                return self._parse_fluxml(root)
            except ET.ParseError as e:
                raise ValueError(f"XML parsing error in {file_path}: {e}")

        except Exception as e:
            raise ValueError(f"Error parsing {file_path}: {e}")

    def _parse_fluxml(self, root: ET.Element) -> FluxML:
        """Parse the root fluxml element."""
        # Handle namespace prefix
        if root.tag.startswith("{"):
            # Namespaced element
            ns_prefix = root.tag.split("}")[0] + "}"
        else:
            # No namespace
            ns_prefix = ""

        # Parse components
        info = self._parse_info(root.find(f"{ns_prefix}info"))
        reactionnetwork = self._parse_reactionnetwork(
            root.find(f"{ns_prefix}reactionnetwork")
        )
        constraints = self._parse_constraints(
            root.find(f"{ns_prefix}constraints")
        )

        # Parse experiments
        experiments = []
        for exp_elem in root.findall(f"{ns_prefix}configuration"):
            experiments.append(self._parse_experiments(exp_elem))

        return FluxML(
            info=info,
            reactionnetwork=reactionnetwork,
            constraints=constraints,
            experiments=experiments,
        )

    def _parse_info(self, info_elem: Optional[ET.Element]) -> Optional[Info]:
        """Parse info element."""
        if info_elem is None:
            return None

        # Get namespace prefix
        ns_prefix = self._get_namespace_prefix(info_elem)

        name = self._get_text(info_elem.find(f"{ns_prefix}name"))
        version = self._get_text(info_elem.find(f"{ns_prefix}version"))
        date_str = self._get_text(info_elem.find(f"{ns_prefix}date"))
        comment = self._get_text(info_elem.find(f"{ns_prefix}comment"))
        modeler = self._get_text(info_elem.find(f"{ns_prefix}modeler"))
        strain = self._get_text(info_elem.find(f"{ns_prefix}strain"))

        # Parse date
        date = None
        if date_str:
            try:
                date = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                pass  # Invalid date format

        return Info(
            name=name,
            version=version,
            date=date,
            comment=comment,
            modeler=modeler,
            strain=strain,
        )

    def _parse_reactionnetwork(self, rn_elem: ET.Element) -> ReactionNetwork:
        """Parse reactionnetwork element."""
        if rn_elem is None:
            raise ValueError("reactionnetwork element is required")

        ns_prefix = self._get_namespace_prefix(rn_elem)

        # Parse metabolites
        pools_elem = rn_elem.find(f"{ns_prefix}metabolitepools")
        if pools_elem is None:
            raise ValueError("metabolitepools element is required")

        metabolites = []
        for pool_elem in pools_elem.findall(f"{ns_prefix}pool"):
            metabolites.append(self._parse_metabolite(pool_elem))

        metabolites_collection = Metabolites(metabolites=metabolites)

        # Parse reactions
        reactions = []
        for reaction_elem in rn_elem.findall(f"{ns_prefix}reaction"):
            reactions.append(self._parse_reaction(reaction_elem))

        return ReactionNetwork(
            metabolites=metabolites_collection, reactions=reactions
        )

    def _parse_metabolite(self, pool_elem: ET.Element) -> Metabolite:
        """Parse metabolite (pool) element."""
        metabolite_id = pool_elem.get("id")
        if not metabolite_id:
            raise ValueError("Metabolite (pool) must have an id attribute")

        atoms = int(pool_elem.get("atoms", "0"))
        size = float(pool_elem.get("size", "1.0"))
        cfg = pool_elem.get("cfg", "0")

        # Parse annotations
        annotations = []
        ns_prefix = self._get_namespace_prefix(pool_elem)
        for ann_elem in pool_elem.findall(f"{ns_prefix}annotation"):
            annotations.append(self._parse_annotation(ann_elem))

        return Metabolite(
            id=metabolite_id,
            atoms=atoms,
            size=size,
            cfg=cfg,
            annotations=annotations,
        )

    def _parse_reaction(self, reaction_elem: ET.Element) -> Reaction:
        """Parse reaction element."""
        reaction_id = reaction_elem.get("id")
        if not reaction_id:
            raise ValueError("Reaction must have an id attribute")

        bidirectional = (
            reaction_elem.get("bidirectional", "true").lower() == "true"
        )

        ns_prefix = self._get_namespace_prefix(reaction_elem)

        # Parse annotations
        annotations = []
        for ann_elem in reaction_elem.findall(f"{ns_prefix}annotation"):
            annotations.append(self._parse_annotation(ann_elem))

        # Parse reducts
        reducts = []
        for reduct_elem in reaction_elem.findall(f"{ns_prefix}reduct"):
            reducts.append(self._parse_reduct(reduct_elem))

        # Parse rproducts
        rproducts = []
        for rproduct_elem in reaction_elem.findall(f"{ns_prefix}rproduct"):
            rproducts.append(self._parse_rproduct(rproduct_elem))

        return Reaction(
            id=reaction_id,
            bidirectional=bidirectional,
            annotations=annotations,
            reducts=reducts,
            rproducts=rproducts,
        )

    def _parse_reduct(self, reduct_elem: ET.Element) -> Reduct:
        """Parse reduct element."""
        reduct_id = reduct_elem.get("id")
        if not reduct_id:
            raise ValueError("Reduct must have an id attribute")

        cfg = reduct_elem.get("cfg")

        # Parse variants
        variants = []
        ns_prefix = self._get_namespace_prefix(reduct_elem)
        for variant_elem in reduct_elem.findall(f"{ns_prefix}variant"):
            variants.append(self._parse_variant(variant_elem))

        return Reduct(id=reduct_id, cfg=cfg, variants=variants)

    def _parse_rproduct(self, rproduct_elem: ET.Element) -> RProduct:
        """Parse rproduct element."""
        rproduct_id = rproduct_elem.get("id")
        if not rproduct_id:
            raise ValueError("RProduct must have an id attribute")

        cfg = rproduct_elem.get("cfg")

        # Parse variants
        variants = []
        ns_prefix = self._get_namespace_prefix(rproduct_elem)
        for variant_elem in rproduct_elem.findall(f"{ns_prefix}variant"):
            variants.append(self._parse_variant(variant_elem))

        return RProduct(id=rproduct_id, cfg=cfg, variants=variants)

    def _parse_variant(self, variant_elem: ET.Element) -> Variant:
        """Parse variant element."""
        cfg = variant_elem.get("cfg")
        if not cfg:
            raise ValueError("Variant must have a cfg attribute")

        ratio_str = variant_elem.get("ratio")
        ratio = float(ratio_str) if ratio_str else None

        return Variant(cfg=cfg, ratio=ratio)

    def _parse_annotation(self, ann_elem: ET.Element) -> Annotation:
        """Parse annotation element."""
        name = ann_elem.get("name")
        if not name:
            raise ValueError("Annotation must have a name attribute")

        content = ann_elem.text
        return Annotation(name=name, content=content)

    def _parse_constraints(
        self, constraints_elem: Optional[ET.Element]
    ) -> Optional[Constraints]:
        """Parse constraints element."""
        if constraints_elem is None:
            return None

        ns_prefix = self._get_namespace_prefix(constraints_elem)

        # Parse net constraints
        net_elem = constraints_elem.find(f"{ns_prefix}net")
        net = (
            self._parse_net_constraints(net_elem)
            if net_elem is not None
            else None
        )

        # Parse xch constraints
        xch_elem = constraints_elem.find(f"{ns_prefix}xch")
        xch = (
            self._parse_xch_constraints(xch_elem)
            if xch_elem is not None
            else None
        )

        # Parse psize constraints
        psize_elem = constraints_elem.find(f"{ns_prefix}psize")
        psize = (
            self._parse_metabolitesize_constraints(psize_elem)
            if psize_elem is not None
            else None
        )

        return Constraints(net=net, xch=xch, metabolitesize=psize)

    def _parse_net_constraints(self, net_elem: ET.Element) -> NetConstraints:
        """Parse net constraints element."""
        expression = self._parse_textual_or_math(net_elem)
        return NetConstraints(expression=expression)

    def _parse_xch_constraints(self, xch_elem: ET.Element) -> XchConstraints:
        """Parse xch constraints element."""
        expression = self._parse_textual_or_math(xch_elem)
        return XchConstraints(expression=expression)

    def _parse_metabolitesize_constraints(
        self, psize_elem: ET.Element
    ) -> MetaboliteSizeConstraints:
        """Parse psize constraints element."""
        expression = self._parse_textual_or_math(psize_elem)
        return MetaboliteSizeConstraints(expression=expression)

    def _parse_textual_or_math(self, elem: ET.Element) -> TextualOrMath:
        """Parse textual or MathML content."""
        ns_prefix = self._get_namespace_prefix(elem)

        textual_elem = elem.find(f"{ns_prefix}textual")
        textual = textual_elem.text if textual_elem is not None else None

        # For MathML, we'd need to serialize the subtree
        math_elem = elem.find(f"{ns_prefix}math")
        mathml = None
        if math_elem is not None:
            mathml = ET.tostring(math_elem, encoding="unicode")

        # If neither textual nor mathml, use element text
        if textual is None and mathml is None:
            textual = elem.text

        return TextualOrMath(textual=textual, mathml=mathml)

    def _parse_experiments(self, config_elem: ET.Element) -> Experiments:
        """Parse configuration element into an Experiments object."""
        name = config_elem.get("name")
        if not name:
            raise ValueError(
                "Experiment (configuration) must have a name attribute"
            )

        stationary = config_elem.get("stationary", "true").lower() == "true"
        time_str = config_elem.get("time")
        time = float(time_str) if time_str else None

        ns_prefix = self._get_namespace_prefix(config_elem)

        # Parse comment
        comment_elem = config_elem.find(f"{ns_prefix}comment")
        comment = comment_elem.text if comment_elem is not None else None

        # Parse tracers (inputs)
        tracers = []
        for input_elem in config_elem.findall(f"{ns_prefix}input"):
            tracers.append(self._parse_tracer(input_elem))

        # Parse constraints (optional)
        constraints_elem = config_elem.find(f"{ns_prefix}constraints")
        constraints = (
            self._parse_constraints(constraints_elem)
            if constraints_elem is not None
            else None
        )

        # Parse measurement (optional)
        measurement_elem = config_elem.find(f"{ns_prefix}measurement")
        measurement = (
            self._parse_measurement(measurement_elem)
            if measurement_elem is not None
            else None
        )

        # Parse simulation (optional)
        simulation_elem = config_elem.find(f"{ns_prefix}simulation")
        simulation = (
            self._parse_simulation(simulation_elem)
            if simulation_elem is not None
            else None
        )

        return Experiments(
            name=name,
            stationary=stationary,
            time=time,
            comment=comment,
            tracers=tracers,
            constraints=constraints,
            measurement=measurement,
            simulation=simulation,
        )

    def _parse_tracer(self, input_elem: ET.Element) -> Tracers:
        """Parse input element into a Tracers object."""
        metabolite = input_elem.get("pool")
        if not metabolite:
            raise ValueError("Tracer (input) must have a pool attribute")

        input_id = input_elem.get("id")
        input_type = input_elem.get("type", "isotopomer")
        profile = input_elem.get("profile")

        # Parse labels
        labels = []
        ns_prefix = self._get_namespace_prefix(input_elem)
        for label_elem in input_elem.findall(f"{ns_prefix}label"):
            labels.append(self._parse_label(label_elem))

        return Tracers(
            id=input_id,
            metabolite=metabolite,
            type=input_type,
            profile=profile,
            labels=labels,
        )

    def _parse_label(self, label_elem: ET.Element) -> Label:
        """Parse label element."""
        cfg = label_elem.get("cfg")
        if not cfg:
            raise ValueError("Label must have a cfg attribute")

        purity = label_elem.get("purity")
        cost_str = label_elem.get("cost")
        cost = float(cost_str) if cost_str else None

        content = label_elem.text

        return Label(cfg=cfg, purity=purity, cost=cost, content=content)

    def _parse_measurement(self, measurement_elem: ET.Element) -> Measurement:
        """Parse measurement element according to FluxML schema."""
        ns_prefix = self._get_namespace_prefix(measurement_elem)

        # Parse model section (required)
        model_elem = measurement_elem.find(f"{ns_prefix}model")
        if model_elem is None:
            raise ValueError("measurement must have a model element")
        model = self._parse_measurement_model(model_elem)

        # Parse data section (required)
        data_elem = measurement_elem.find(f"{ns_prefix}data")
        if data_elem is None:
            raise ValueError("measurement must have a data element")
        data = self._parse_measurement_data(data_elem)

        return Measurement(model=model, data=data)

    def _parse_measurement_model(
        self, model_elem: ET.Element
    ) -> MeasurementModel:
        """Parse measurement model element."""
        ns_prefix = self._get_namespace_prefix(model_elem)

        # Parse labelingmeasurement (optional)
        labeling_elem = model_elem.find(f"{ns_prefix}labelingmeasurement")
        labeling_measurement = None
        if labeling_elem is not None:
            labeling_measurement = self._parse_labeling_measurement(
                labeling_elem
            )

        # Parse fluxmeasurement (optional)
        flux_elem = model_elem.find(f"{ns_prefix}fluxmeasurement")
        flux_measurement = None
        if flux_elem is not None:
            flux_measurement = self._parse_flux_measurement(flux_elem)

        # Parse poolsizemeasurement (optional)
        poolsize_elem = model_elem.find(f"{ns_prefix}poolsizemeasurement")
        metabolitesize_measurement = None
        if poolsize_elem is not None:
            metabolitesize_measurement = self._parse_metabolitesize_measurement(
                poolsize_elem
            )

        return MeasurementModel(
            labeling_measurement=labeling_measurement,
            flux_measurement=flux_measurement,
            metabolitesize_measurement=metabolitesize_measurement,
        )

    def _parse_labeling_measurement(
        self, labeling_elem: ET.Element
    ) -> LabelingMeasurement:
        """Parse labelingmeasurement element."""
        ns_prefix = self._get_namespace_prefix(labeling_elem)

        groups = []
        for group_elem in labeling_elem.findall(f"{ns_prefix}group"):
            groups.append(self._parse_group(group_elem))

        return LabelingMeasurement(groups=groups)

    def _parse_group(self, group_elem: ET.Element) -> Group:
        """Parse group element."""
        group_id = group_elem.get("id")
        if not group_id:
            raise ValueError("Group must have an id attribute")

        times = group_elem.get("times")
        scale = group_elem.get("scale", "auto")

        ns_prefix = self._get_namespace_prefix(group_elem)

        # Parse errormodel (optional)
        errormodel_elem = group_elem.find(f"{ns_prefix}errormodel")
        errormodel = None
        if errormodel_elem is not None:
            errormodel = ErrorModel(
                expression=self._parse_textual_or_math(errormodel_elem)
            )

        # Parse expression (textual or math)
        expression = self._parse_textual_or_math(group_elem)

        return Group(
            id=group_id,
            times=times,
            scale=scale,
            errormodel=errormodel,
            expression=expression,
        )

    def _parse_flux_measurement(self, flux_elem: ET.Element) -> FluxMeasurement:
        """Parse fluxmeasurement element."""
        ns_prefix = self._get_namespace_prefix(flux_elem)

        net_fluxes = []
        for netflux_elem in flux_elem.findall(f"{ns_prefix}netflux"):
            net_fluxes.append(self._parse_netflux(netflux_elem))

        xch_fluxes = []
        for xchflux_elem in flux_elem.findall(f"{ns_prefix}xchflux"):
            xch_fluxes.append(self._parse_xchflux(xchflux_elem))

        return FluxMeasurement(net_fluxes=net_fluxes, xch_fluxes=xch_fluxes)

    def _parse_netflux(self, netflux_elem: ET.Element) -> NetFlux:
        """Parse netflux element."""
        netflux_id = netflux_elem.get("id")
        if not netflux_id:
            raise ValueError("NetFlux must have an id attribute")

        ns_prefix = self._get_namespace_prefix(netflux_elem)

        # Parse errormodel (optional)
        errormodel_elem = netflux_elem.find(f"{ns_prefix}errormodel")
        errormodel = None
        if errormodel_elem is not None:
            errormodel = ErrorModel(
                expression=self._parse_textual_or_math(errormodel_elem)
            )

        # Parse expression
        expression = self._parse_textual_or_math(netflux_elem)

        return NetFlux(
            id=netflux_id, errormodel=errormodel, expression=expression
        )

    def _parse_xchflux(self, xchflux_elem: ET.Element) -> XchFlux:
        """Parse xchflux element."""
        xchflux_id = xchflux_elem.get("id")
        if not xchflux_id:
            raise ValueError("XchFlux must have an id attribute")

        ns_prefix = self._get_namespace_prefix(xchflux_elem)

        # Parse errormodel (optional)
        errormodel_elem = xchflux_elem.find(f"{ns_prefix}errormodel")
        errormodel = None
        if errormodel_elem is not None:
            errormodel = ErrorModel(
                expression=self._parse_textual_or_math(errormodel_elem)
            )

        # Parse expression
        expression = self._parse_textual_or_math(xchflux_elem)

        return XchFlux(
            id=xchflux_id, errormodel=errormodel, expression=expression
        )

    def _parse_metabolitesize_measurement(
        self, poolsize_elem: ET.Element
    ) -> MetaboliteSizeMeasurement:
        """Parse poolsizemeasurement element."""
        ns_prefix = self._get_namespace_prefix(poolsize_elem)

        metabolite_sizes = []
        for poolsize_item_elem in poolsize_elem.findall(f"{ns_prefix}poolsize"):
            metabolite_sizes.append(
                self._parse_metabolitesize(poolsize_item_elem)
            )

        return MetaboliteSizeMeasurement(metabolite_sizes=metabolite_sizes)

    def _parse_metabolitesize(
        self, poolsize_elem: ET.Element
    ) -> MetaboliteSize:
        """Parse poolsize element."""
        metabolitesize_id = poolsize_elem.get("id")
        if not metabolitesize_id:
            raise ValueError(
                "MetaboliteSize (poolsize) must have an id attribute"
            )

        ns_prefix = self._get_namespace_prefix(poolsize_elem)

        # Parse errormodel (optional)
        errormodel_elem = poolsize_elem.find(f"{ns_prefix}errormodel")
        errormodel = None
        if errormodel_elem is not None:
            errormodel = ErrorModel(
                expression=self._parse_textual_or_math(errormodel_elem)
            )

        # Parse expression
        expression = self._parse_textual_or_math(poolsize_elem)

        return MetaboliteSize(
            id=metabolitesize_id, errormodel=errormodel, expression=expression
        )

    def _parse_measurement_data(self, data_elem: ET.Element) -> MeasurementData:
        """Parse measurement data element."""
        ns_prefix = self._get_namespace_prefix(data_elem)

        data = []
        for datum_elem in data_elem.findall(f"{ns_prefix}datum"):
            data.append(self._parse_datum(datum_elem))

        return MeasurementData(data=data)

    def _parse_datum(self, datum_elem: ET.Element) -> Datum:
        """Parse datum element."""
        datum_id = datum_elem.get("id")
        if not datum_id:
            raise ValueError("Datum must have an id attribute")

        stddev_str = datum_elem.get("stddev")
        if not stddev_str:
            raise ValueError("Datum must have a stddev attribute")
        stddev = float(stddev_str.strip())

        # Parse optional attributes
        row_str = datum_elem.get("row")
        row = int(row_str) if row_str else None

        time_str = datum_elem.get("time")
        time = float(time_str) if time_str else None

        weight = datum_elem.get("weight")

        pos_str = datum_elem.get("pos")
        pos = int(pos_str) if pos_str else None

        datum_type = datum_elem.get("type")

        # Parse value from element text
        value_text = datum_elem.text
        if not value_text:
            raise ValueError("Datum must have a value in its text content")
        value = float(value_text.strip())

        return Datum(
            id=datum_id,
            value=value,
            stddev=stddev,
            row=row,
            time=time,
            weight=weight,
            pos=pos,
            type=datum_type,
        )

    def _parse_simulation(self, simulation_elem: ET.Element) -> Simulation:
        """Parse simulation element (simplified)."""
        sim_type = simulation_elem.get("type", "auto")
        method = simulation_elem.get("method", "auto")

        return Simulation(type=sim_type, method=method)

    def _get_namespace_prefix(self, elem: ET.Element) -> str:
        """Get namespace prefix for element."""
        if elem.tag.startswith("{"):
            return elem.tag.split("}")[0] + "}"
        return ""

    def _get_text(self, elem: Optional[ET.Element]) -> Optional[str]:
        """Get text content from element."""
        if elem is None:
            return None
        return elem.text.strip() if elem.text else None


def parse_fluxml_file(file_path: str) -> FluxML:
    """Parse a FluxML file and return a FluxML object."""
    parser = FluxMLParser()
    return parser.parse_file(file_path)
