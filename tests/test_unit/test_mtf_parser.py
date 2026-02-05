"""
Tests for MTF (influx_si) format parser.
"""

import pytest
from pathlib import Path

from fluxomics_data_model.io import MTFParser, parse_mtf
from fluxomics_data_model.core.core import FluxomicsDataModel


# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "MTF_tests"


class TestMTFParserBasic:
    """Basic MTF parser functionality tests."""

    def test_parse_ecoli_network(self):
        """Test parsing E.coli network file."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        assert isinstance(model, FluxomicsDataModel)
        assert model.info.name == "e_coli"

        # Check reactions were parsed
        assert len(model.model.reactions) == 86
        assert len(model.model.metabolites) == 81
        assert len(model.model.atom_mappings) == 86

    def test_parse_ecoli_constraints(self):
        """Test parsing E.coli constraints."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        assert model.constraints is not None
        assert model.constraints.net is not None
        assert len(model.constraints.net.formulas) == 11

        # Check XCH constraints
        assert model.constraints.xch is not None
        assert len(model.constraints.xch.formulas) == 1

    def test_parse_ecoli_tracers(self):
        """Test parsing E.coli tracer specifications."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        assert len(model.experiments) == 1
        experiment = model.experiments[0]

        # Check tracers
        assert len(experiment.tracers) > 0

        # Find Gluc_U tracer (uniformly labeled glucose)
        gluc_u_tracer = None
        for tracer in experiment.tracers:
            if tracer.metabolite == "Gluc_U":
                gluc_u_tracer = tracer
                break

        assert gluc_u_tracer is not None
        assert len(gluc_u_tracer.labels) == 2  # 111111 and 000000

    def test_parse_ecoli_measurements(self):
        """Test parsing E.coli measurements."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        experiment = model.experiments[0]
        assert experiment.measurement is not None

        # Check labeling measurements (MS data)
        assert experiment.measurement.model.labeling_measurement is not None
        groups = experiment.measurement.model.labeling_measurement.groups
        assert len(groups) > 0

        # Check flux measurements
        assert experiment.measurement.model.flux_measurement is not None
        net_fluxes = experiment.measurement.model.flux_measurement.net_fluxes
        assert len(net_fluxes) == 1  # out_Ac flux measurement


class TestMTFParserReactions:
    """Tests for reaction parsing from .netw files."""

    def test_reaction_reversibility(self):
        """Test that reversible reactions are correctly identified."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        # pgi should be reversible (<->)
        pgi = model.model.reactions.get_by_id("pgi")
        assert pgi.reversibility is True

        # pfk should be irreversible (->)
        pfk = model.model.reactions.get_by_id("pfk")
        assert pfk.reversibility is False

    def test_reaction_reactants_products(self):
        """Test that reactants and products are correctly parsed."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        # pgi: Glc6P <-> Fru6P
        pgi = model.model.reactions.get_by_id("pgi")
        assert pgi.reactants == ["Glc6P"]
        assert pgi.products == ["Fru6P"]

        # ald: FruBP <-> GA3P + GA3P
        ald = model.model.reactions.get_by_id("ald")
        assert ald.reactants == ["FruBP"]
        assert ald.products == ["GA3P", "GA3P"]

    def test_atom_mapping_created(self):
        """Test that atom mappings are created for reactions."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        # Check atom mapping exists
        assert "pgi" in model.model.atom_mappings
        pgi_mapping = model.model.atom_mappings["pgi"]

        assert pgi_mapping.reaction_id == "pgi"
        assert pgi_mapping.reactants == ["Glc6P"]
        assert pgi_mapping.products == ["Fru6P"]


class TestMTFParserConstraints:
    """Tests for constraint parsing from .cnstr files."""

    def test_equality_constraints(self):
        """Test parsing equality constraints."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        # Find fum_a-fum_b == 0 constraint
        net_formulas = model.constraints.net.formulas
        fum_constraint = None
        for formula in net_formulas:
            if "fum_a-fum_b" in formula.expression:
                fum_constraint = formula
                break

        assert fum_constraint is not None
        assert "= 0" in fum_constraint.expression

    def test_inequality_constraints(self):
        """Test parsing inequality constraints."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        # Find pyk >= 1 constraint
        net_formulas = model.constraints.net.formulas
        pyk_constraint = None
        for formula in net_formulas:
            if "pyk" in formula.expression and ">=" in formula.expression:
                pyk_constraint = formula
                break

        assert pyk_constraint is not None
        assert ">= 1" in pyk_constraint.expression


class TestMTFParserConvenienceFunctions:
    """Tests for convenience functions and edge cases."""

    def test_parse_with_extension(self):
        """Test that parsing works with full file path including extension."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli.netw")

        assert isinstance(model, FluxomicsDataModel)
        assert len(model.model.reactions) == 86

    def test_parser_instance(self):
        """Test using MTFParser class directly."""
        parser = MTFParser()
        model = parser.parse(TEST_DATA_DIR / "Ecoli" / "e_coli")

        assert isinstance(model, FluxomicsDataModel)

    def test_missing_network_file(self):
        """Test error when network file is missing."""
        with pytest.raises(FileNotFoundError):
            parse_mtf(TEST_DATA_DIR / "nonexistent" / "model")

    def test_model_repr(self):
        """Test model string representation."""
        model = parse_mtf(TEST_DATA_DIR / "Ecoli" / "e_coli")

        repr_str = repr(model)
        assert "Fluxomics Data Model Summary" in repr_str
        assert "e_coli" in repr_str
        assert "Reactions" in repr_str
