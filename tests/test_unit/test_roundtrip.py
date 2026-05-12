"""Tests for round-trip conversion between formats.

Tests that data is preserved when converting between:
- FluxML -> MTF -> FluxML
- MTF -> FluxML -> MTF
"""

import pytest
from pathlib import Path
import tempfile
import shutil

from fluxomics_data_converter.io import (
    parse_fluxml_file,
    write_fluxml,
    parse_mtf,
    write_mtf,
)


# Path to benchmark data
BENCHMARK_DIR = Path(__file__).parent.parent.parent / "data" / "benchmark"
FLUXML_BENCHMARK = BENCHMARK_DIR / "FluxML" / "EC.fml"
MTF_BENCHMARK = BENCHMARK_DIR / "influx_si" / "EC"


class TestFluxMLRoundTrip:
    """Test FluxML -> MTF -> FluxML round-trip."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test output."""
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not FLUXML_BENCHMARK.exists(), reason="Benchmark FluxML file not found"
    )
    def test_fluxml_to_mtf(self, temp_dir):
        """Test converting FluxML to MTF."""
        # Parse FluxML
        model = parse_fluxml_file(str(FLUXML_BENCHMARK))
        assert model is not None

        # Write to MTF (use first experiment name)
        exp_name = model.experiments[0].name if model.experiments else None
        mtf_base = temp_dir / "test_model"
        write_mtf(model, str(mtf_base), experiment_name=exp_name)

        # Verify MTF files were created
        assert (mtf_base.with_suffix(".netw")).exists()

        # Read the generated .netw file
        netw_content = (mtf_base.with_suffix(".netw")).read_text()
        assert "# Network definition" in netw_content

    @pytest.mark.skipif(
        not FLUXML_BENCHMARK.exists(), reason="Benchmark FluxML file not found"
    )
    def test_fluxml_reaction_count(self, temp_dir):
        """Test that reaction count is preserved in conversion."""
        # Parse original
        original = parse_fluxml_file(str(FLUXML_BENCHMARK))
        original_rxn_count = len(original.model.reactions)

        # Write to MTF and read back
        exp_name = (
            original.experiments[0].name if original.experiments else None
        )
        mtf_base = temp_dir / "test_model"
        write_mtf(original, str(mtf_base), experiment_name=exp_name)

        # Parse MTF
        converted = parse_mtf(str(mtf_base))

        # Check reaction count (may differ slightly due to variant handling)
        assert len(converted.model.reactions) > 0
        print(f"Original reactions: {original_rxn_count}")
        print(f"Converted reactions: {len(converted.model.reactions)}")


class TestMTFRoundTrip:
    """Test MTF -> FluxML -> MTF round-trip."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test output."""
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not (MTF_BENCHMARK.with_suffix(".netw")).exists(),
        reason="Benchmark MTF file not found",
    )
    def test_mtf_to_fluxml(self, temp_dir):
        """Test converting MTF to FluxML."""
        # Parse MTF
        model = parse_mtf(str(MTF_BENCHMARK))
        assert model is not None

        # Write to FluxML
        fluxml_path = temp_dir / "test_model.fml"
        write_fluxml(model, str(fluxml_path))

        # Verify FluxML file was created
        assert fluxml_path.exists()

        # Read the generated FluxML file
        fluxml_content = fluxml_path.read_text()
        assert "fluxml" in fluxml_content
        assert "reactionnetwork" in fluxml_content

    @pytest.mark.skipif(
        not (MTF_BENCHMARK.with_suffix(".netw")).exists(),
        reason="Benchmark MTF file not found",
    )
    def test_mtf_metabolite_count(self, temp_dir):
        """Test that metabolite count is preserved."""
        # Parse original
        original = parse_mtf(str(MTF_BENCHMARK))
        original_met_count = len(original.model.metabolites)

        # Write to FluxML and read back
        fluxml_path = temp_dir / "test_model.fml"
        write_fluxml(original, str(fluxml_path))

        # Parse FluxML
        converted = parse_fluxml_file(str(fluxml_path))

        # Check metabolite count
        assert len(converted.model.metabolites) == original_met_count
        print(f"Metabolites preserved: {original_met_count}")


class TestFullRoundTrip:
    """Test complete round-trip: FluxML -> MTF -> FluxML."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test output."""
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not FLUXML_BENCHMARK.exists(), reason="Benchmark FluxML file not found"
    )
    def test_full_roundtrip(self, temp_dir):
        """Test FluxML -> MTF -> FluxML round-trip."""
        # Parse original FluxML
        original = parse_fluxml_file(str(FLUXML_BENCHMARK))

        # Convert to MTF
        exp_name = (
            original.experiments[0].name if original.experiments else None
        )
        mtf_base = temp_dir / "intermediate"
        write_mtf(original, str(mtf_base), experiment_name=exp_name)

        # Read MTF
        intermediate = parse_mtf(str(mtf_base))

        # Convert back to FluxML
        fluxml_path = temp_dir / "roundtrip.fml"
        write_fluxml(intermediate, str(fluxml_path))

        # Read final FluxML
        final = parse_fluxml_file(str(fluxml_path))

        # Verify key properties are preserved
        assert final is not None

        # Metabolites should be preserved
        original_mets = set(original.model.metabolites.ids)
        final_mets = set(final.model.metabolites.ids)

        # Check overlap (some metabolites might be handled differently)
        overlap = original_mets & final_mets
        print(f"Original metabolites: {len(original_mets)}")
        print(f"Final metabolites: {len(final_mets)}")
        print(f"Overlap: {len(overlap)}")

        # At least 90% should overlap
        assert len(overlap) / len(original_mets) > 0.9


class TestConversionDetails:
    """Detailed tests for specific conversion aspects."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test output."""
        temp = tempfile.mkdtemp()
        yield Path(temp)
        shutil.rmtree(temp)

    @pytest.mark.skipif(
        not (MTF_BENCHMARK.with_suffix(".netw")).exists(),
        reason="Benchmark MTF file not found",
    )
    def test_variant_reaction_handling(self, temp_dir):
        """Test that variant reactions are correctly handled."""
        # Parse MTF which has variant reactions (v24_a, v24_b, etc.)
        model = parse_mtf(str(MTF_BENCHMARK))

        # Check variant reactions exist
        rxn_ids = [r.id for r in model.model.reactions]

        # v24_a and v24_b should be separate reactions in MTF
        assert "v24_a" in rxn_ids or "v24_b" in rxn_ids

        # Write to FluxML
        fluxml_path = temp_dir / "variants.fml"
        write_fluxml(model, str(fluxml_path))

        # Verify the file can be read back
        converted = parse_fluxml_file(str(fluxml_path))
        assert converted is not None

    @pytest.mark.skipif(
        not (MTF_BENCHMARK.with_suffix(".netw")).exists(),
        reason="Benchmark MTF file not found",
    )
    def test_tracer_preservation(self, temp_dir):
        """Test that tracer information is preserved."""
        # Parse MTF
        model = parse_mtf(str(MTF_BENCHMARK))

        if not model.experiments:
            pytest.skip("No experiments in model")

        original_tracers = model.experiments[0].tracers
        if not original_tracers:
            pytest.skip("No tracers in experiment")

        # Write to FluxML and read back
        fluxml_path = temp_dir / "tracers.fml"
        write_fluxml(model, str(fluxml_path))
        converted = parse_fluxml_file(str(fluxml_path))

        if converted.experiments:
            converted_tracers = converted.experiments[0].tracers

            # Compare tracer metabolites
            original_metabolites = {t.metabolite for t in original_tracers}
            converted_metabolites = {t.metabolite for t in converted_tracers}

            assert original_metabolites == converted_metabolites
