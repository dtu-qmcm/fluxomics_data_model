"""Tests for MTF (influx_si) format parser."""

import pytest
from pathlib import Path

from fluxomics_data_converter.io import MTFParser, parse_mtf
from fluxomics_data_converter.core.core import FluxomicsData


# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "MTF_tests"
TOY_BASE = TEST_DATA_DIR / "toy_test" / "toy"


class TestMTFParserBasic:
    """Basic MTF parser functionality tests."""

    def test_parse_toy_network(self):
        """Test parsing the toy network file."""
        model = parse_mtf(TOY_BASE)

        assert isinstance(model, FluxomicsData)
        assert model.metadata.name == "toy"

        # V1-V19 with V4 split into two variants (V4_1, V4_2) -> 20 reactions.
        assert len(model.model.reactions) == 20
        assert len(model.model.metabolites) == 20
        # V17-V19 carry no atom transitions, so they get no atom mapping.
        assert len(model.model.atom_mappings) == 17

    def test_parse_toy_constraints(self):
        """Test parsing toy constraints."""
        model = parse_mtf(TOY_BASE)

        assert model.constraints is not None
        assert model.constraints.net is not None
        # Two biomass-coupling equalities, V1 fixed, and V4_1-V4_2 symmetry.
        # The commented '# NET V4 >= 0.2' line must be skipped.
        assert len(model.constraints.net.formulas) == 4

        # The toy model has no exchange (XCH) constraints.
        assert model.constraints.xch is None

    def test_parse_toy_tracers(self):
        """Test parsing toy tracer specifications."""
        model = parse_mtf(TOY_BASE)

        assert len(model.experiments) == 1
        experiment = model.experiments[0]

        assert len(experiment.tracers) == 3

        # Metabolite A is labeled with two isotopomers: 100 (0.8) and 111 (0.2).
        a_tracer = None
        for tracer in experiment.tracers:
            if tracer.metabolite == "A":
                a_tracer = tracer
                break

        assert a_tracer is not None
        assert len(a_tracer.labels) == 2

    def test_parse_toy_measurements(self):
        """Test parsing toy measurements."""
        model = parse_mtf(TOY_BASE)

        experiment = model.experiments[0]
        assert experiment.measurement is not None

        # Labeling (MS) measurements from the .miso file.
        assert experiment.measurement.model.labeling_measurement is not None
        groups = experiment.measurement.model.labeling_measurement.groups
        assert len(groups) == 4

        # Flux measurements from the .mflux file (BIOMASS_YIELD, V11, V8).
        assert experiment.measurement.model.flux_measurement is not None
        net_fluxes = experiment.measurement.model.flux_measurement.net_fluxes
        assert len(net_fluxes) == 3


class TestMTFParserReactions:
    """Tests for reaction parsing from .netw files."""

    def test_reaction_reversibility(self):
        """Test that reversible reactions are correctly identified."""
        model = parse_mtf(TOY_BASE)

        # V3 is reversible (<->).
        v3 = model.model.reactions.get_by_id("V3")
        assert v3.reversibility is True

        # V1 is irreversible (->).
        v1 = model.model.reactions.get_by_id("V1")
        assert v1.reversibility is False

    def test_reaction_reactants_products(self):
        """Test that reactants and products are correctly parsed."""
        model = parse_mtf(TOY_BASE)

        # V1: A -> B
        v1 = model.model.reactions.get_by_id("V1")
        assert v1.reactants == ["A"]
        assert v1.products == ["B"]

        # V15: I <-> L + F
        v15 = model.model.reactions.get_by_id("V15")
        assert v15.reactants == ["I"]
        assert v15.products == ["L", "F"]

    def test_atom_mapping_created(self):
        """Test that atom transitions are created for reactions."""
        model = parse_mtf(TOY_BASE)

        # Check atom transition exists.
        assert "V1" in model.model.atom_mappings
        v1_mapping = model.model.atom_mappings["V1"]

        assert v1_mapping.reaction_id == "V1"
        assert v1_mapping.reactants == ["A"]
        assert v1_mapping.products == ["B"]


class TestMTFParserConstraints:
    """Tests for constraint parsing from .cnstr files."""

    def test_equality_constraints(self):
        """Test parsing equality constraints."""
        model = parse_mtf(TOY_BASE)

        # Find the V4_1 - V4_2 == 0 symmetry constraint.
        net_formulas = model.constraints.net.formulas
        v4_constraint = None
        for formula in net_formulas:
            if "V4_1 - V4_2" in formula.expression:
                v4_constraint = formula
                break

        assert v4_constraint is not None
        assert "= 0" in v4_constraint.expression

    def test_commented_constraints_skipped(self):
        """Commented (#) constraint lines must not be parsed."""
        model = parse_mtf(TOY_BASE)

        # The '# NET V4 >= 0.2' line is commented out, so no inequality
        # constraint should have been parsed.
        net_formulas = model.constraints.net.formulas
        assert all(">=" not in f.expression for f in net_formulas)
        assert all("<=" not in f.expression for f in net_formulas)


class TestMTFParserConvenienceFunctions:
    """Tests for convenience functions and edge cases."""

    def test_parse_with_extension(self):
        """Test that parsing works with full file path including extension."""
        model = parse_mtf(TEST_DATA_DIR / "toy_test" / "toy.netw")

        assert isinstance(model, FluxomicsData)
        assert len(model.model.reactions) == 20

    def test_parser_instance(self):
        """Test using MTFParser class directly."""
        parser = MTFParser()
        model = parser.parse(TOY_BASE)

        assert isinstance(model, FluxomicsData)

    def test_missing_network_file(self):
        """Test error when network file is missing."""
        with pytest.raises(FileNotFoundError):
            parse_mtf(TEST_DATA_DIR / "nonexistent" / "model")

    def test_model_repr(self):
        """Test model string representation."""
        model = parse_mtf(TOY_BASE)

        repr_str = repr(model)
        assert "Fluxomics Data Converter Summary" in repr_str
        assert "toy" in repr_str
        assert "Reactions" in repr_str
