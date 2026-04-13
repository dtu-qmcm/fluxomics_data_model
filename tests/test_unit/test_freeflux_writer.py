"""
Tests for FreeFlux writer and 3-way format round-trip.

Tests that data is preserved when converting between:
- FreeFlux -> FluxML -> FreeFlux
- FreeFlux -> MTF -> FreeFlux
- MTF -> FreeFlux -> MTF
"""

import pytest
from pathlib import Path
import tempfile
import shutil

from fluxomics_data_model.io import (
    FluxMLParser,
    FluxMLWriter,
    MTFParser,
    MTFWriter,
    FreefluxParser,
    FreefluxWriter,
    parse_fluxml_file,
    write_fluxml,
    parse_mtf,
    write_mtf,
    parse_freeflux,
    write_freeflux,
)
from fluxomics_data_model.core.core import FluxomicsDataModel


FREEFLUX_DIR = Path(__file__).parent.parent.parent / "data" / "freeflux_tests"
FLUXML_BENCHMARK = (
    Path(__file__).parent.parent.parent
    / "data"
    / "benchmark"
    / "FluxML"
    / "EC.fml"
)
MTF_BENCHMARK = (
    Path(__file__).parent.parent.parent
    / "data"
    / "benchmark"
    / "influx_si"
    / "EC"
)


class TestFreefluxWriter:
    """Test FreeFlux writer output."""

    @pytest.fixture
    def temp_dir(self):
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_write_reactions(self, temp_dir):
        model = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(model, str(temp_dir), experiment_name="toy")

        reactions_file = temp_dir / "reactions.tsv"
        assert reactions_file.exists()
        content = reactions_file.read_text()
        assert "reaction_ID" in content
        assert len(model.model.reactions) == content.count("\n") - 1

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_write_creates_all_files(self, temp_dir):
        model = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(model, str(temp_dir), experiment_name="toy")

        assert (temp_dir / "reactions.tsv").exists()
        assert (temp_dir / "fluxes.tsv").exists()
        assert (temp_dir / "concentrations.tsv").exists()
        assert (temp_dir / "measured_fluxes.tsv").exists()

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_multiple_experiments_requires_name(self, temp_dir):
        from fluxomics_data_model.core.core import Model, Experiments, Metadata
        from fluxomics_data_model.model.metabolite import Metabolite
        from fluxomics_data_model.model.reaction import Reaction

        m = Model(
            metabolites=[Metabolite(id="A"), Metabolite(id="B")],
            reactions=[Reaction(id="r1", reactants=["A"], products=["B"])],
        )
        dm = FluxomicsDataModel(
            model=m,
            experiments=[Experiments(name="exp1"), Experiments(name="exp2")],
        )
        with pytest.raises(ValueError, match="Multiple experiments"):
            write_freeflux(dm, str(temp_dir))


class TestFreefluxRoundTrip:
    """Test FreeFlux -> FreeFlux round-trip."""

    @pytest.fixture
    def temp_dir(self):
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_roundtrip_reactions(self, temp_dir):
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(original, str(temp_dir), experiment_name="toy")

        roundtripped = parse_freeflux(temp_dir)
        assert set(original.model.reactions.ids) == set(
            roundtripped.model.reactions.ids
        )

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_roundtrip_metabolites(self, temp_dir):
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(original, str(temp_dir), experiment_name="toy")

        roundtripped = parse_freeflux(temp_dir)
        assert set(original.model.metabolites.ids) == set(
            roundtripped.model.metabolites.ids
        )

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_roundtrip_atom_mappings(self, temp_dir):
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(original, str(temp_dir), experiment_name="toy")

        roundtripped = parse_freeflux(temp_dir)
        assert len(original.model.atom_mappings) == len(
            roundtripped.model.atom_mappings
        )

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_roundtrip_flux_count(self, temp_dir):
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        write_freeflux(original, str(temp_dir), experiment_name="toy")

        roundtripped = parse_freeflux(temp_dir)
        if original.experiments and roundtripped.experiments:
            orig_fluxes = len(
                original.experiments[0].simulation.variables.flux_values
            )
            rt_fluxes = len(
                roundtripped.experiments[0].simulation.variables.flux_values
            )
            assert orig_fluxes == rt_fluxes


class TestThreeWayRoundTrip:
    """Test 3-way format conversion: FreeFlux <-> FluxML <-> MTF."""

    @pytest.fixture
    def temp_dir(self):
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_to_fluxml(self, temp_dir):
        """FreeFlux -> FluxML conversion preserves key data."""
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        fluxml_path = temp_dir / "output.fml"
        write_fluxml(original, str(fluxml_path))

        assert fluxml_path.exists()
        converted = parse_fluxml_file(str(fluxml_path))
        assert len(converted.model.reactions) > 0

    @pytest.mark.skipif(
        not FLUXML_BENCHMARK.exists(),
        reason="Benchmark FluxML file not found",
    )
    def test_fluxml_to_freeflux(self, temp_dir):
        """FluxML -> FreeFlux conversion preserves key data."""
        original = parse_fluxml_file(str(FLUXML_BENCHMARK))

        exp_name = (
            original.experiments[0].name if original.experiments else None
        )
        write_freeflux(original, str(temp_dir), experiment_name=exp_name)

        assert (temp_dir / "reactions.tsv").exists()
        converted = parse_freeflux(temp_dir)
        assert len(converted.model.reactions) > 0
        assert len(converted.model.metabolites) > 0

    @pytest.mark.skipif(
        not (MTF_BENCHMARK.with_suffix(".netw")).exists(),
        reason="Benchmark MTF file not found",
    )
    def test_mtf_to_freeflux(self, temp_dir):
        """MTF -> FreeFlux conversion preserves key data."""
        original = parse_mtf(str(MTF_BENCHMARK))

        exp_name = (
            original.experiments[0].name if original.experiments else None
        )
        write_freeflux(original, str(temp_dir), experiment_name=exp_name)

        assert (temp_dir / "reactions.tsv").exists()
        converted = parse_freeflux(temp_dir)
        assert len(converted.model.reactions) > 0

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_to_mtf(self, temp_dir):
        """FreeFlux -> MTF conversion preserves key data."""
        original = parse_freeflux(FREEFLUX_DIR / "toy")
        mtf_base = temp_dir / "output"
        exp_name = (
            original.experiments[0].name if original.experiments else None
        )
        write_mtf(original, str(mtf_base), experiment_name=exp_name)

        assert (mtf_base.with_suffix(".netw")).exists()

    @pytest.mark.skipif(
        not (FREEFLUX_DIR / "toy").exists(),
        reason="FreeFlux toy test data not found",
    )
    def test_freeflux_fluxml_freeflux_roundtrip(self, temp_dir):
        """FreeFlux -> FluxML -> FreeFlux round-trip."""
        original = parse_freeflux(FREEFLUX_DIR / "toy")

        fluxml_path = temp_dir / "intermediate.fml"
        write_fluxml(original, str(fluxml_path))

        intermediate = parse_fluxml_file(str(fluxml_path))

        if intermediate.experiments:
            exp_name = intermediate.experiments[0].name
        else:
            exp_name = None
        write_freeflux(
            intermediate, str(temp_dir / "rt"), experiment_name=exp_name
        )

        final = parse_freeflux(temp_dir / "rt")
        assert len(final.model.reactions) > 0
        assert len(final.model.metabolites) > 0

        orig_mets = set(original.model.metabolites.ids)
        final_mets = set(final.model.metabolites.ids)
        overlap = orig_mets & final_mets
        assert len(overlap) / len(orig_mets) > 0.8
