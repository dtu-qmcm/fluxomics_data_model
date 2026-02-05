"""
Tests for Freeflux tabular format parser.
"""

import pytest
from pathlib import Path

from fluxomics_data_model.io import FreefluxParser, parse_freeflux
from fluxomics_data_model.core.core import FluxomicsDataModel


# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "freeflux_test"


class TestFreefluxParserToy:
    """Tests for parsing the toy example (TSV files)."""

    def test_parse_toy_network(self):
        """Test parsing toy network from TSV files."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        assert isinstance(model, FluxomicsDataModel)
        assert model.info.name == "toy"

        # Check reactions were parsed
        assert len(model.model.reactions) > 0
        assert len(model.model.metabolites) > 0

    def test_parse_toy_reactions(self):
        """Test that reactions are correctly parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Check that we have reactions with atom mappings
        assert len(model.model.atom_mappings) > 0

    def test_parse_toy_fluxes(self):
        """Test parsing flux values from toy example."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Should have an experiment with simulation data
        assert len(model.experiments) > 0
        experiment = model.experiments[0]

        # Check simulation has flux values
        if experiment.simulation:
            flux_values = experiment.simulation.variables.flux_values
            assert len(flux_values) > 0

    def test_parse_toy_concentrations(self):
        """Test parsing metabolite concentrations."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        experiment = model.experiments[0]

        # Check simulation has metabolite size values
        if experiment.simulation:
            metab_values = experiment.simulation.variables.metabolitesize_values
            assert len(metab_values) > 0

    def test_parse_toy_measured_mdvs(self):
        """Test parsing steady-state MDV measurements."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        experiment = model.experiments[0]
        assert experiment.measurement is not None

        # Check labeling measurements exist
        labeling = experiment.measurement.model.labeling_measurement
        assert labeling is not None
        assert len(labeling.groups) > 0

        # Check measurement data
        data = experiment.measurement.data.data
        assert len(data) > 0

    def test_parse_toy_measured_fluxes(self):
        """Test parsing measured flux values."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        experiment = model.experiments[0]
        assert experiment.measurement is not None

        # Check flux measurements
        flux_meas = experiment.measurement.model.flux_measurement
        assert flux_meas is not None
        assert len(flux_meas.net_fluxes) > 0

    def test_parse_toy_inst_mdvs(self):
        """Test parsing time-course MDV measurements."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        experiment = model.experiments[0]

        # Toy example has inst_MDVs, so should be non-stationary
        assert experiment.stationary is False


class TestFreefluxParserEcoli:
    """Tests for parsing E. coli example (XLSX files)."""

    def test_parse_ecoli_synthetic(self):
        """Test parsing E. coli synthetic data."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "synthetic_data")

        assert isinstance(model, FluxomicsDataModel)
        assert model.info.name == "synthetic_data"

        # Check reactions were parsed
        assert len(model.model.reactions) > 0
        assert len(model.model.metabolites) > 0

    def test_parse_ecoli_experimental(self):
        """Test parsing E. coli experimental data."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "experimental_data")

        assert isinstance(model, FluxomicsDataModel)

        # Check reactions were parsed
        assert len(model.model.reactions) > 0

    def test_ecoli_mdv_fragments(self):
        """Test that MDV fragments are correctly parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "synthetic_data")

        experiment = model.experiments[0]
        labeling = experiment.measurement.model.labeling_measurement

        # Check some known fragments exist
        group_ids = [g.id for g in labeling.groups]

        # E. coli should have amino acid fragments
        assert any("Ala" in gid for gid in group_ids)
        assert any("Glu" in gid for gid in group_ids)


class TestFreefluxParserSynechocystis:
    """Tests for parsing Synechocystis example (XLSX with inst_MDVs)."""

    def test_parse_synechocystis_synthetic(self):
        """Test parsing Synechocystis synthetic data."""
        model = parse_freeflux(TEST_DATA_DIR / "synechocystis" / "synthetic_data")

        assert isinstance(model, FluxomicsDataModel)

        # Check reactions were parsed
        assert len(model.model.reactions) > 0

    def test_synechocystis_inst_mdvs(self):
        """Test that time-course MDVs are parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "synechocystis" / "synthetic_data")

        experiment = model.experiments[0]

        # Synechocystis has inst_MDVs
        assert experiment.stationary is False

        # Check labeling measurements have time info
        labeling = experiment.measurement.model.labeling_measurement
        assert labeling is not None

        # Check data includes time points
        data = experiment.measurement.data.data
        time_values = [d.time for d in data if d.time is not None]
        assert len(time_values) > 0


class TestFreefluxParserAtomMappings:
    """Tests for atom mapping parsing."""

    def test_atom_mapping_created(self):
        """Test that atom mappings are created for reactions."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Check that atom mappings exist
        assert len(model.model.atom_mappings) > 0

    def test_atom_mapping_structure(self):
        """Test atom mapping structure."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Get first atom mapping
        first_mapping = list(model.model.atom_mappings.values())[0]

        assert first_mapping.reaction_id is not None
        assert len(first_mapping.reactants) > 0
        assert len(first_mapping.products) > 0
        assert len(first_mapping.maps) > 0


class TestFreefluxParserEdgeCases:
    """Tests for edge cases and error handling."""

    def test_missing_reactions_file(self):
        """Test error when reactions file is missing."""
        with pytest.raises(FileNotFoundError):
            parse_freeflux(TEST_DATA_DIR / "nonexistent")

    def test_parser_instance(self):
        """Test using FreefluxParser class directly."""
        parser = FreefluxParser()
        model = parser.parse(TEST_DATA_DIR / "toy")

        assert isinstance(model, FluxomicsDataModel)

    def test_model_repr(self):
        """Test model string representation."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        repr_str = repr(model)
        assert "Fluxomics Data Model Summary" in repr_str
        assert "Reactions" in repr_str


class TestFreefluxParserMDVParsing:
    """Tests for MDV value parsing."""

    def test_mdv_values_parsed(self):
        """Test that MDV values are correctly parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "synthetic_data")

        experiment = model.experiments[0]
        data = experiment.measurement.data.data

        # Check that we have MDV components (M0, M1, M2, etc.)
        mdv_data = [d for d in data if d.pos is not None]
        assert len(mdv_data) > 0

        # Check values are in valid range [0, 1] for MDVs
        for datum in mdv_data:
            if datum.value is not None:
                assert 0 <= datum.value <= 1, f"MDV value {datum.value} out of range"

    def test_mdv_stddev_parsed(self):
        """Test that MDV standard deviations are parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "synthetic_data")

        experiment = model.experiments[0]
        data = experiment.measurement.data.data

        # Check that stddev values exist
        for datum in data:
            if datum.pos is not None:
                assert datum.stddev is not None
                assert datum.stddev > 0


class TestFreefluxFormatConsistency:
    """Tests to ensure Freeflux format produces consistent results with other formats."""

    def test_reaction_has_required_fields(self):
        """Test that reactions have all required fields."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        for reaction in model.model.reactions:
            assert reaction.id is not None
            assert reaction.reactants is not None
            assert reaction.products is not None
            assert reaction.reversibility is not None

    def test_metabolite_has_required_fields(self):
        """Test that metabolites have all required fields."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        for metabolite in model.model.metabolites:
            assert metabolite.id is not None

    def test_flux_value_structure(self):
        """Test that flux values have proper structure."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        experiment = model.experiments[0]
        if experiment.simulation:
            for flux_val in experiment.simulation.variables.flux_values:
                assert flux_val.flux is not None
                assert flux_val.value is not None
