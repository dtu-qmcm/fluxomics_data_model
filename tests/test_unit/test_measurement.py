import pytest
import jax.numpy as jnp

from fluxomics_data_model.experiment.measurement import (
    Group,
    Datum,
    MeasurementData,
    Measurement,
    MeasurementModel,
    LabelingMeasurement,
    FluxMeasurement,
    NetFlux,
    ExchangeFlux,
    MetaboliteSize,
    MetaboliteSizeMeasurement,
)
from fluxomics_data_model.core.common import TextualOrMath


class TestGroup:
    def test_creation(self):
        g = Group(id="Ala_158", expression=TextualOrMath(textual="Ala_158"))
        assert g.id == "Ala_158"

    def test_with_times(self):
        g = Group(
            id="Ala_158",
            times="0.5, 1.0, 2.0",
            expression=TextualOrMath(textual="Ala_158"),
        )
        assert g.times == "0.5, 1.0, 2.0"

    def test_time_points_property(self):
        g = Group(
            id="g1", times="0.5, 1.0", expression=TextualOrMath(textual="x")
        )
        tp = g.time_points
        assert tp is not None
        assert jnp.allclose(tp, jnp.array([0.5, 1.0]))

    def test_duplicate_group_ids_raises(self):
        with pytest.raises(ValueError, match="unique"):
            LabelingMeasurement(
                groups=[
                    Group(id="g1", expression=TextualOrMath(textual="x")),
                    Group(id="g1", expression=TextualOrMath(textual="y")),
                ]
            )


class TestNetFlux:
    def test_creation(self):
        nf = NetFlux(id="out_Ac", expression=TextualOrMath(textual="out_Ac"))
        assert nf.id == "out_Ac"


class TestExchangeFlux:
    def test_creation(self):
        ef = ExchangeFlux(
            id="PGI_xch", expression=TextualOrMath(textual="PGI_xch")
        )
        assert ef.id == "PGI_xch"


class TestMetaboliteSize:
    def test_creation(self):
        ms = MetaboliteSize(id="ATP", expression=TextualOrMath(textual="ATP"))
        assert ms.id == "ATP"


class TestDatum:
    def test_creation(self):
        d = Datum(id="Ala_158_M0", value=0.5, stddev=0.01)
        assert d.id == "Ala_158_M0"
        assert d.value == 0.5
        assert d.stddev == 0.01

    def test_with_optional_fields(self):
        d = Datum(id="x", value=1.0, stddev=0.1, time=2.0, pos=0)
        assert d.time == 2.0
        assert d.pos == 0


class TestMeasurementData:
    def test_values_property(self):
        data = MeasurementData(
            data=[
                Datum(id="a", value=0.5, stddev=0.01),
                Datum(id="b", value=0.3, stddev=0.02),
            ]
        )
        vals = data.values
        assert jnp.allclose(vals, jnp.array([0.5, 0.3]))

    def test_errors_property(self):
        data = MeasurementData(
            data=[
                Datum(id="a", value=0.5, stddev=0.01),
                Datum(id="b", value=0.3, stddev=0.02),
            ]
        )
        errs = data.errors
        assert jnp.allclose(errs, jnp.array([0.01, 0.02]))

    def test_times_property(self):
        data = MeasurementData(
            data=[
                Datum(id="a", value=0.5, stddev=0.01, time=1.0),
                Datum(id="b", value=0.3, stddev=0.02, time=2.0),
            ]
        )
        times = data.times
        assert jnp.allclose(times, jnp.array([1.0, 2.0]))

    def test_times_property_no_times(self):
        data = MeasurementData(
            data=[
                Datum(id="a", value=0.5, stddev=0.01),
            ]
        )
        assert data.times is None

    def test_get_data_for_id(self):
        data = MeasurementData(
            data=[
                Datum(id="a", value=0.5, stddev=0.01),
                Datum(id="b", value=0.3, stddev=0.02),
                Datum(id="a", value=0.4, stddev=0.01),
            ]
        )
        result = data.get_data_for_id("a")
        assert len(result) == 2


class TestMeasurementModel:
    def test_creation(self):
        mm = MeasurementModel()
        assert mm.labeling_measurement is None
        assert mm.flux_measurement is None

    def test_with_labeling(self):
        mm = MeasurementModel(
            labeling_measurement=LabelingMeasurement(
                groups=[
                    Group(id="g1", expression=TextualOrMath(textual="x")),
                ]
            )
        )
        assert len(mm.labeling_measurement.groups) == 1

    def test_with_flux_measurement(self):
        mm = MeasurementModel(
            flux_measurement=FluxMeasurement(
                net_fluxes=[
                    NetFlux(id="r1", expression=TextualOrMath(textual="r1"))
                ],
            )
        )
        assert len(mm.flux_measurement.net_fluxes) == 1


class TestMeasurement:
    def test_creation(self):
        m = Measurement(
            model=MeasurementModel(),
            data=MeasurementData(),
        )
        assert m.model is not None
        assert m.data is not None
