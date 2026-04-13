import pytest
import jax.numpy as jnp

from fluxomics_data_model.output.simulation import (
    FluxValue,
    MetaboliteSizeValue,
    Variables,
    Simulation,
)


class TestFluxValue:
    def test_creation(self):
        fv = FluxValue(flux="r1", type="net")
        assert fv.flux == "r1"
        assert fv.type == "net"
        assert fv.value is None

    def test_with_bounds(self):
        fv = FluxValue(flux="r1", type="net", lo=-10.0, hi=100.0)
        assert fv.bounds == (-10.0, 100.0)

    def test_bounds_array(self):
        fv = FluxValue(flux="r1", type="net", lo=-5.0, hi=50.0)
        arr = fv.bounds_array
        assert jnp.allclose(arr, jnp.array([-5.0, 50.0]))

    def test_default_bounds_array(self):
        fv = FluxValue(flux="r1", type="net")
        arr = fv.bounds_array
        assert arr[0] == float("-inf")
        assert arr[1] == float("inf")


class TestMetaboliteSizeValue:
    def test_creation(self):
        mv = MetaboliteSizeValue(metabolite="ATP")
        assert mv.metabolite == "ATP"

    def test_with_bounds(self):
        mv = MetaboliteSizeValue(metabolite="ATP", lo=0.1, hi=10.0)
        assert mv.bounds == (0.1, 10.0)

    def test_bounds_array_default_lo(self):
        mv = MetaboliteSizeValue(metabolite="ATP")
        arr = mv.bounds_array
        assert arr[0] == 0.0
        assert arr[1] == float("inf")


class TestVariables:
    def test_creation(self):
        v = Variables()
        assert v.flux_values == []
        assert v.metabolitesize_values == []

    def test_with_flux_values(self):
        v = Variables(
            flux_values=[
                FluxValue(flux="r1", type="net", value=1.0),
                FluxValue(flux="r2", type="net", value=2.0),
            ]
        )
        assert len(v.flux_values) == 2

    def test_duplicate_flux_ids_raises(self):
        with pytest.raises(ValueError, match="unique"):
            Variables(
                flux_values=[
                    FluxValue(flux="r1", type="net", value=1.0),
                    FluxValue(flux="r1", type="net", value=2.0),
                ]
            )

    def test_get_flux_bounds_matrix(self):
        v = Variables(
            flux_values=[
                FluxValue(flux="r1", type="net", lo=-10.0, hi=100.0),
            ]
        )
        mat = v.get_flux_bounds_matrix(["r1"])
        assert mat.shape == (1, 2)
        assert jnp.allclose(mat[0], jnp.array([-10.0, 100.0]))

    def test_get_flux_bounds_missing_reaction(self):
        v = Variables(
            flux_values=[
                FluxValue(flux="r1", type="net", lo=0.0, hi=10.0),
            ]
        )
        mat = v.get_flux_bounds_matrix(["r1", "r2"])
        assert mat.shape == (2, 2)
        assert mat[1, 0] == float("-inf")

    def test_duplicate_metabolite_size_ids_raises(self):
        with pytest.raises(ValueError, match="unique"):
            Variables(
                metabolitesize_values=[
                    MetaboliteSizeValue(metabolite="ATP"),
                    MetaboliteSizeValue(metabolite="ATP"),
                ]
            )


class TestSimulation:
    def test_creation(self):
        sim = Simulation()
        assert sim.type == "auto"
        assert sim.method == "auto"

    def test_with_variables(self):
        sim = Simulation(
            variables=Variables(
                flux_values=[
                    FluxValue(flux="r1", type="net", value=1.0),
                ]
            )
        )
        assert len(sim.variables.flux_values) == 1

    def test_get_optimization_bounds(self):
        sim = Simulation(
            variables=Variables(
                flux_values=[
                    FluxValue(flux="r1", type="net", lo=0.0, hi=10.0),
                ]
            )
        )
        fb, mb = sim.get_optimization_bounds(["r1"], ["ATP"])
        assert fb.shape == (1, 2)
        assert mb.shape == (1, 2)

    def test_get_optimization_bounds_no_variables(self):
        sim = Simulation()
        fb, mb = sim.get_optimization_bounds(["r1"], ["ATP"])
        assert fb.shape == (1, 2)
        assert mb.shape == (1, 2)
