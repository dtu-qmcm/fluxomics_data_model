"""Tests for Freeflux tabular format parser."""

import pytest
from pathlib import Path

from fluxomics_data_converter.io import FreefluxParser, parse_freeflux
from fluxomics_data_converter.core.core import FluxomicsData


# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "freeflux_tests"


class TestFreefluxParserToy:
    """Tests for parsing the toy example (TSV files)."""

    def test_parse_toy_network(self):
        """Test parsing toy network from TSV files."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        assert isinstance(model, FluxomicsData)
        assert model.metadata.name == "toy"

        # Check reactions were parsed
        assert len(model.model.reactions) > 0
        assert len(model.model.metabolites) > 0

    def test_parse_toy_reactions(self):
        """Test that reactions are correctly parsed."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Check that we have reactions with atom transitions
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

        assert isinstance(model, FluxomicsData)
        assert model.metadata.name == "synthetic_data"

        # Check reactions were parsed
        assert len(model.model.reactions) > 0
        assert len(model.model.metabolites) > 0

    def test_parse_ecoli_experimental(self):
        """Test parsing E. coli experimental data."""
        model = parse_freeflux(TEST_DATA_DIR / "ecoli" / "experimental_data")

        assert isinstance(model, FluxomicsData)

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
        model = parse_freeflux(
            TEST_DATA_DIR / "synechocystis" / "synthetic_data"
        )

        assert isinstance(model, FluxomicsData)

        # Check reactions were parsed
        assert len(model.model.reactions) > 0

    def test_synechocystis_inst_mdvs(self):
        """Test that time-course MDVs are parsed."""
        model = parse_freeflux(
            TEST_DATA_DIR / "synechocystis" / "synthetic_data"
        )

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


class TestFreefluxParserAtomTransitions:
    """Tests for atom transition parsing."""

    def test_atom_mapping_created(self):
        """Test that atom transitions are created for reactions."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Check that atom transitions exist
        assert len(model.model.atom_mappings) > 0

    def test_atom_mapping_structure(self):
        """Test atom transition structure."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        # Get first atom transition
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

        assert isinstance(model, FluxomicsData)

    def test_model_repr(self):
        """Test model string representation."""
        model = parse_freeflux(TEST_DATA_DIR / "toy")

        repr_str = repr(model)
        assert "Fluxomics Data Converter Summary" in repr_str
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
                assert (
                    0 <= datum.value <= 1
                ), f"MDV value {datum.value} out of range"

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_tsv(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n")


def _minimal_reactions_tsv(directory: Path) -> None:
    """Write the smallest valid reactions.tsv so parse() can proceed."""
    _write_tsv(
        directory / "reactions.tsv",
        [
            "#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility",
            "R1\tA (abc)\tB (abc)\t0",
        ],
    )


# ---------------------------------------------------------------------------
# _parse_compounds — space-separated format (the bug we fixed)
# ---------------------------------------------------------------------------


class TestParseCompounds:
    """Unit-tests for _parse_compounds; exercises the \\s* regex fix."""

    def setup_method(self):
        self.parser = FreefluxParser()

    def test_compact_format(self):
        """Name(atoms) — original format still works."""
        ids, atoms = self.parser._parse_compounds("A(abc)")
        assert ids == ["A"]
        assert atoms == [("A", "abc")]

    def test_space_before_parens(self):
        """Name (atoms) — space between name and parentheses."""
        ids, atoms = self.parser._parse_compounds("A (abc)")
        assert ids == ["A"]
        assert atoms == [("A", "abc")]

    def test_multiple_compounds_with_spaces(self):
        """B (abc) + C (de) — multiple space-separated compounds."""
        ids, atoms = self.parser._parse_compounds("B (abc) + C (de)")
        assert ids == ["B", "C"]
        assert dict(atoms) == {"B": "abc", "C": "de"}

    def test_repeated_compound(self):
        """Same metabolite appearing twice (e.g. E (a) + E (e))."""
        ids, atoms = self.parser._parse_compounds("E (a) + E (e)")
        assert ids == ["E", "E"]
        assert len(atoms) == 2

    def test_symmetric_atoms(self):
        """SUCC(abcd,dcba) — comma variants preserved."""
        ids, atoms = self.parser._parse_compounds("SUCC(abcd,dcba)")
        assert ids == ["SUCC"]
        assert atoms[0][1] == "abcd,dcba"

    def test_no_atoms(self):
        """Compound with no atom annotation still yields a compound ID."""
        ids, atoms = self.parser._parse_compounds("CO2")
        assert ids == ["CO2"]
        assert atoms == []


# ---------------------------------------------------------------------------
# _find_file — list patterns and alias resolution
# ---------------------------------------------------------------------------


class TestFindFile:
    """Unit-tests for _find_file; exercises list-of-stems behaviour."""

    def setup_method(self):
        self.parser = FreefluxParser()

    def test_finds_primary_name(self, tmp_path):
        (tmp_path / "reactions.tsv").touch()
        result = self.parser._find_file(tmp_path, "reactions")
        assert result == tmp_path / "reactions.tsv"

    def test_finds_singular_alias(self, tmp_path):
        """'reaction.tsv' should be found when looking for 'reactions'."""
        (tmp_path / "reaction.tsv").touch()
        result = self.parser._find_file(tmp_path, "reactions")
        assert result == tmp_path / "reaction.tsv"

    def test_finds_measured_mid_alias(self, tmp_path):
        """measured_MID.tsv must resolve for the 'measured_MDVs' key."""
        (tmp_path / "measured_MID.tsv").touch()
        result = self.parser._find_file(tmp_path, "measured_MDVs")
        assert result == tmp_path / "measured_MID.tsv"

    def test_finds_measured_mdv_alias(self, tmp_path):
        """measured_MDV.tsv (no s) must also resolve."""
        (tmp_path / "measured_MDV.tsv").touch()
        result = self.parser._find_file(tmp_path, "measured_MDVs")
        assert result == tmp_path / "measured_MDV.tsv"

    def test_finds_inst_mid_alias(self, tmp_path):
        """measured_inst_MID.tsv must resolve for the inst key."""
        (tmp_path / "measured_inst_MID.tsv").touch()
        result = self.parser._find_file(tmp_path, "measured_inst_MDVs")
        assert result == tmp_path / "measured_inst_MID.tsv"

    def test_returns_none_when_absent(self, tmp_path):
        result = self.parser._find_file(tmp_path, "measured_MDVs")
        assert result is None

    def test_primary_name_preferred_over_alias(self, tmp_path):
        """When both exist, the primary (first) name takes precedence."""
        (tmp_path / "measured_MDVs.tsv").touch()
        (tmp_path / "measured_MID.tsv").touch()
        result = self.parser._find_file(tmp_path, "measured_MDVs")
        assert result == tmp_path / "measured_MDVs.tsv"

    def test_csv_extension_found(self, tmp_path):
        """Parser resolves .csv as well as .tsv."""
        (tmp_path / "fluxes.csv").touch()
        result = self.parser._find_file(tmp_path, "fluxes")
        assert result == tmp_path / "fluxes.csv"


# ---------------------------------------------------------------------------
# _parse_label_input — new function
# ---------------------------------------------------------------------------


class TestParseLabelInput:
    """Tests for the new _parse_label_input method."""

    def setup_method(self):
        self.parser = FreefluxParser()

    def test_returns_empty_when_no_file(self, tmp_path):
        result = self.parser._parse_label_input(tmp_path)
        assert result == []

    def test_basic_single_pattern(self, tmp_path):
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t111\t1.0\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert len(tracers) == 1
        tracer = tracers[0]
        assert tracer.metabolite == "A"
        assert len(tracer.labels) == 1
        assert tracer.labels[0].labeled_pattern == "111"
        assert tracer.labels[0].fraction == pytest.approx(1.0)
        assert tracer.labels[0].purity == pytest.approx(1.0)

    def test_quoted_pattern_stripped(self, tmp_path):
        """Patterns stored with surrounding quotes (e.g. '010') are cleaned."""
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t'010'\t1.0\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert tracers[0].labels[0].labeled_pattern == "010"

    def test_comma_separated_patterns_in_one_row(self, tmp_path):
        """Multiple patterns in a single cell become multiple LabelCompositions."""
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "AcCoA\t01,11\t0.25,0.25\t1.0,1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert len(tracers) == 1
        assert len(tracers[0].labels) == 2
        assert tracers[0].labels[0].labeled_pattern == "01"
        assert tracers[0].labels[1].labeled_pattern == "11"
        assert tracers[0].labels[0].fraction == pytest.approx(0.25)
        assert tracers[0].labels[1].fraction == pytest.approx(0.25)

    def test_multiple_rows_same_metabolite_merged(self, tmp_path):
        """Two rows for the same metabolite are merged into one Tracers entry."""
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "Glc\t100000\t0.5\t1.0",
                "Glc\t111111\t0.5\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert len(tracers) == 1
        assert tracers[0].metabolite == "Glc"
        assert len(tracers[0].labels) == 2
        patterns = [lb.labeled_pattern for lb in tracers[0].labels]
        assert "100000" in patterns
        assert "111111" in patterns

    def test_multiple_substrates(self, tmp_path):
        """Different metabolites produce separate Tracers objects."""
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t111\t1.0\t1.0",
                "B\t000\t1.0\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert len(tracers) == 2
        metabolites = {t.metabolite for t in tracers}
        assert metabolites == {"A", "B"}

    def test_tracer_type_is_isotopomer(self, tmp_path):
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t1\t1.0\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert tracers[0].type == "isotopomer"

    def test_label_input_singular_filename(self, tmp_path):
        """label_inputs.tsv (plural) is also discovered."""
        _write_tsv(
            tmp_path / "label_inputs.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t111\t1.0\t1.0",
            ],
        )
        tracers = self.parser._parse_label_input(tmp_path)
        assert len(tracers) == 1


# ---------------------------------------------------------------------------
# Integration: label_input wired into parse()
# ---------------------------------------------------------------------------


class TestLabelInputIntegration:
    """Verify that label_input feeds through parse() into LabelingExperiments."""

    def test_tracers_present_in_experiment(self, tmp_path):
        _minimal_reactions_tsv(tmp_path)
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t111\t1.0\t1.0",
            ],
        )
        model = parse_freeflux(tmp_path)
        assert len(model.experiments) == 1
        exp = model.experiments[0]
        assert len(exp.tracers) == 1
        assert exp.tracers[0].metabolite == "A"

    def test_tracers_alone_create_experiment(self, tmp_path):
        """An experiment is created even when there are only tracers (no measurements)."""
        _minimal_reactions_tsv(tmp_path)
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t1\t1.0\t1.0",
            ],
        )
        model = parse_freeflux(tmp_path)
        assert len(model.experiments) == 1
        assert model.experiments[0].measurement is None
        assert model.experiments[0].simulation is None

    def test_tracers_combined_with_measurements(self, tmp_path):
        _minimal_reactions_tsv(tmp_path)
        _write_tsv(
            tmp_path / "label_input.tsv",
            [
                "#metabolite_id\tlabeling_pattern\tpercentage\tpurity",
                "A\t111\t1.0\t1.0",
            ],
        )
        _write_tsv(
            tmp_path / "measured_MID.tsv",
            [
                "#fragment_ID\tmean\tsd",
                "B_123\t0.1,0.9,0.0\t0.01,0.01,0.01",
            ],
        )
        model = parse_freeflux(tmp_path)
        exp = model.experiments[0]
        assert len(exp.tracers) == 1
        assert exp.measurement is not None
        labeling = exp.measurement.model.labeling_measurement
        assert labeling is not None
        assert len(labeling.groups) == 1


# ---------------------------------------------------------------------------
# _parse_flux_bounds — new function
# ---------------------------------------------------------------------------


class TestParseFluxBounds:
    """Tests for _parse_flux_bounds / constraints file parsing."""

    def setup_method(self):
        self.parser = FreefluxParser()
        # Pre-populate reactions so 'all' expansion works
        from fluxomics_data_converter.model.reaction import Reaction

        self.parser._reactions = {
            "R1": Reaction(
                id="R1", reactants=["A"], products=["B"], reversibility=False
            ),
            "R2": Reaction(
                id="R2", reactants=["B"], products=["C"], reversibility=False
            ),
        }

    def test_returns_none_when_no_file(self, tmp_path):
        result = self.parser._parse_flux_bounds(tmp_path)
        assert result is None

    def test_basic_bounds(self, tmp_path):
        _write_tsv(
            tmp_path / "flux_bounds.tsv",
            ["#reaction_id\tlo\thi", "R1\t0\t100"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        assert constraints is not None
        exprs = [f.expression for f in constraints.net.formulas]
        assert "R1 >= 0.0" in exprs
        assert "R1 <= 100.0" in exprs

    def test_all_expands_to_all_reactions(self, tmp_path):
        _write_tsv(
            tmp_path / "constraints.tsv",
            ["#reaction_id\tlo\thi", "all\t0\t200"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        exprs = [f.expression for f in constraints.net.formulas]
        assert "R1 >= 0.0" in exprs
        assert "R1 <= 200.0" in exprs
        assert "R2 >= 0.0" in exprs
        assert "R2 <= 200.0" in exprs
        assert len(exprs) == 4  # 2 bounds × 2 reactions

    def test_constraints_filename_alias(self, tmp_path):
        """File named constraints.tsv is found for the flux_bounds key."""
        _write_tsv(
            tmp_path / "constraints.tsv",
            ["#reaction_id\tlo\thi", "R1\t5\t50"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        assert constraints is not None

    def test_one_sided_lower_bound(self, tmp_path):
        """Missing hi column produces only a >= constraint."""
        _write_tsv(
            tmp_path / "flux_bounds.tsv",
            ["#reaction_id\tlo", "R1\t10"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        exprs = [f.expression for f in constraints.net.formulas]
        assert "R1 >= 10.0" in exprs
        assert not any("<=" in e for e in exprs)

    def test_one_sided_upper_bound(self, tmp_path):
        """Missing lo column produces only a <= constraint."""
        _write_tsv(
            tmp_path / "flux_bounds.tsv",
            ["#reaction_id\thi", "R1\t50"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        exprs = [f.expression for f in constraints.net.formulas]
        assert "R1 <= 50.0" in exprs
        assert not any(">=" in e for e in exprs)

    def test_multiple_reactions(self, tmp_path):
        _write_tsv(
            tmp_path / "flux_bounds.tsv",
            ["#reaction_id\tlo\thi", "R1\t0\t100", "R2\t10\t80"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        exprs = {f.expression for f in constraints.net.formulas}
        assert "R1 >= 0.0" in exprs
        assert "R1 <= 100.0" in exprs
        assert "R2 >= 10.0" in exprs
        assert "R2 <= 80.0" in exprs

    def test_returns_net_constraints(self, tmp_path):
        """Flux bounds go into NetConstraints; xch and metabolitesize are None."""
        _write_tsv(
            tmp_path / "flux_bounds.tsv",
            ["#reaction_id\tlo\thi", "R1\t0\t100"],
        )
        constraints = self.parser._parse_flux_bounds(tmp_path)
        assert constraints.net is not None
        assert constraints.xch is None
        assert constraints.metabolitesize is None


# ---------------------------------------------------------------------------
# Integration: flux_bounds wired into parse() → FluxomicsData.constraints
# ---------------------------------------------------------------------------


class TestFluxBoundsIntegration:
    """Verify that flux_bounds feeds through parse() into FluxomicsData.constraints."""

    def test_constraints_present_on_model(self, tmp_path):
        _minimal_reactions_tsv(tmp_path)
        _write_tsv(
            tmp_path / "constraints.tsv",
            ["#reaction_id\tlo\thi", "all\t0\t100"],
        )
        model = parse_freeflux(tmp_path)
        assert model.constraints is not None
        assert model.constraints.net is not None
        exprs = [f.expression for f in model.constraints.net.formulas]
        assert any(">= 0" in e for e in exprs)
        assert any("<= 100" in e for e in exprs)

    def test_all_reactions_get_bounds(self, tmp_path):
        """'all' expands to every reaction parsed from reactions.tsv."""
        _write_tsv(
            tmp_path / "reactions.tsv",
            [
                "#reaction_ID\treactant_IDs(atom)\tproduct_IDs(atom)\treversibility",
                "R1\tA (abc)\tB (abc)\t0",
                "R2\tB (abc)\tC (abc)\t0",
            ],
        )
        _write_tsv(
            tmp_path / "constraints.tsv",
            ["#reaction_id\tlo\thi", "all\t0\t50"],
        )
        model = parse_freeflux(tmp_path)
        exprs = {f.expression for f in model.constraints.net.formulas}
        assert "R1 >= 0.0" in exprs
        assert "R1 <= 50.0" in exprs
        assert "R2 >= 0.0" in exprs
        assert "R2 <= 50.0" in exprs

    def test_no_constraints_file_gives_none(self, tmp_path):
        _minimal_reactions_tsv(tmp_path)
        model = parse_freeflux(tmp_path)
        assert model.constraints is None
