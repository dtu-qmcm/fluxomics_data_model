import pytest
import jax.numpy as jnp

from fluxomics_data_model.experiment.tracer import Tracers, LabelComposition


class TestLabelComposition:
    def test_creation(self):
        lc = LabelComposition(labeled_pattern="111111")
        assert lc.labeled_pattern == "111111"

    def test_with_purity(self):
        lc = LabelComposition(labeled_pattern="111111", purity=0.99)
        assert lc.purity == 0.99

    def test_with_fraction(self):
        lc = LabelComposition(labeled_pattern="111111", fraction=0.5)
        assert lc.fraction == 0.5

    def test_invalid_pattern(self):
        with pytest.raises(Exception):
            LabelComposition(labeled_pattern="abc")

    def test_purity_out_of_range(self):
        with pytest.raises(Exception):
            LabelComposition(labeled_pattern="111111", purity=1.5)

    def test_zero_pattern(self):
        lc = LabelComposition(labeled_pattern="000000")
        assert lc.labeled_pattern == "000000"

    def test_mixed_pattern(self):
        lc = LabelComposition(labeled_pattern="110000")
        assert lc.labeled_pattern == "110000"


class TestTracers:
    def test_creation(self):
        t = Tracers(metabolite="glc")
        assert t.metabolite == "glc"
        assert t.labels == []

    def test_with_labels(self):
        t = Tracers(
            metabolite="glc",
            labels=[
                LabelComposition(labeled_pattern="111111", fraction=0.8),
                LabelComposition(labeled_pattern="000000", fraction=0.2),
            ],
        )
        assert len(t.labels) == 2

    def test_fraction_sum_exceeds_one_raises(self):
        with pytest.raises(ValueError, match="exceeds 1.0"):
            Tracers(
                metabolite="glc",
                labels=[
                    LabelComposition(labeled_pattern="111111", fraction=0.6),
                    LabelComposition(labeled_pattern="000000", fraction=0.5),
                ],
            )

    def test_fraction_sum_less_than_one_warns(self):
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Tracers(
                metabolite="glc",
                labels=[
                    LabelComposition(labeled_pattern="111111", fraction=0.3),
                ],
            )
            assert len(w) == 1
            assert (
                "sum to" in str(w[0].message).lower()
                or "expected" in str(w[0].message).lower()
            )

    def test_composition_vector_default(self):
        t = Tracers(metabolite="glc")
        vec = t.composition_vector
        assert vec.shape == (1,)

    def test_with_composition(self):
        t = Tracers(metabolite="glc")
        comp = jnp.array([0.2, 0.8])
        t2 = t.with_composition(comp)
        assert jnp.allclose(t2.composition_vector, comp)

    def test_with_time_profile(self):
        t = Tracers(metabolite="glc")
        times = jnp.array([0.0, 1.0, 2.0])
        values = jnp.array([1.0, 0.5, 0.0])
        t2 = t.with_time_profile(times, values)
        assert t2.time_profile_data is not None

    def test_with_id(self):
        t = Tracers(id="tracer1", metabolite="glc")
        assert t.id == "tracer1"
