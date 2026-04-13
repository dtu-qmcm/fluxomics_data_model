import pytest
import jax.numpy as jnp

from fluxomics_data_model.model.reaction import Reaction


class TestReactionCreation:
    def test_basic_creation(self):
        r = Reaction(id="PGI", reactants=["G6P"], products=["F6P"])
        assert r.id == "PGI"
        assert r.reactants == ["G6P"]
        assert r.products == ["F6P"]
        assert r.reversibility is True

    def test_irreversible(self):
        r = Reaction(
            id="PFK", reactants=["F6P"], products=["FBP"], reversibility=False
        )
        assert r.reversibility is False

    def test_with_name(self):
        r = Reaction(
            id="PGI",
            name="Phosphoglucose isomerase",
            reactants=["G6P"],
            products=["F6P"],
        )
        assert r.name == "Phosphoglucose isomerase"

    def test_multiple_reactants_products(self):
        r = Reaction(
            id="ALD",
            reactants=["FBP"],
            products=["DHAP", "GAP"],
        )
        assert r.reactants == ["FBP"]
        assert r.products == ["DHAP", "GAP"]

    def test_frozen(self):
        r = Reaction(id="PGI", reactants=["G6P"], products=["F6P"])
        with pytest.raises(Exception):
            r.id = "changed"


class TestReactionProperties:
    def test_reactant_ids(self):
        r = Reaction(id="r1", reactants=["A", "B"], products=["C"])
        assert r.reactant_ids == frozenset(["A", "B"])

    def test_product_ids(self):
        r = Reaction(id="r1", reactants=["A"], products=["B", "C"])
        assert r.product_ids == frozenset(["B", "C"])

    def test_metabolites(self):
        r = Reaction(id="r1", reactants=["A", "B"], products=["B", "C"])
        assert r.metabolites == frozenset(["A", "B", "C"])


class TestReactionEquation:
    def test_reversible_equation(self):
        r = Reaction(
            id="PGI", reactants=["G6P"], products=["F6P"], reversibility=True
        )
        assert r.equation == "G6P <=> F6P"

    def test_irreversible_equation(self):
        r = Reaction(
            id="PFK", reactants=["F6P"], products=["FBP"], reversibility=False
        )
        assert r.equation == "F6P => FBP"

    def test_multi_metabolite_equation(self):
        r = Reaction(
            id="r1",
            reactants=["A", "B"],
            products=["C", "D"],
            reversibility=True,
        )
        assert r.equation == "A + B <=> C + D"


class TestReactionVariants:
    def test_no_variants(self):
        r = Reaction(id="PGI", reactants=["G6P"], products=["F6P"])
        assert r.atom_mapping_ids is None
        assert r.is_variant_reaction is False
        assert r.n_variants == 1

    def test_with_variants(self):
        r = Reaction(
            id="SCS",
            reactants=["SUCCOA"],
            products=["SUC"],
            atom_mapping_ids=["SCS___1", "SCS___2"],
        )
        assert r.atom_mapping_ids == ["SCS___1", "SCS___2"]
        assert r.is_variant_reaction is True
        assert r.n_variants == 2

    def test_single_variant(self):
        r = Reaction(
            id="r1",
            reactants=["A"],
            products=["B"],
            atom_mapping_ids=["r1___1"],
        )
        assert r.is_variant_reaction is False
        assert r.n_variants == 1


class TestReactionStoichiometry:
    def test_default_stoichiometric_vector(self):
        r = Reaction(id="r1", reactants=["A"], products=["B"])
        vec = r.get_stoichiometric_vector(["A", "B", "C"])
        assert jnp.allclose(vec, jnp.array([-1.0, 1.0, 0.0]))

    def test_custom_stoichiometry(self):
        r = Reaction(
            id="r1",
            reactants=["A"],
            products=["B"],
            stoichiometry_dict={"A": -2.0, "B": 1.0, "C": 0.0},
        )
        vec = r.get_stoichiometric_vector(["A", "B", "C"])
        assert jnp.allclose(vec, jnp.array([-2.0, 1.0, 0.0]))

    def test_stoichiometric_vector_unrelated_metabolite(self):
        r = Reaction(id="r1", reactants=["A"], products=["B"])
        vec = r.get_stoichiometric_vector(["X", "A", "B"])
        assert jnp.allclose(vec, jnp.array([0.0, -1.0, 1.0]))


class TestReactionFluxBounds:
    def test_reversible_default_bounds(self):
        r = Reaction(
            id="r1", reactants=["A"], products=["B"], reversibility=True
        )
        bounds = r.flux_bounds
        assert bounds[0] < 0
        assert bounds[1] > 0

    def test_irreversible_default_bounds(self):
        r = Reaction(
            id="r1", reactants=["A"], products=["B"], reversibility=False
        )
        bounds = r.flux_bounds
        assert bounds[0] == 0.0

    def test_with_flux_bounds(self):
        r = Reaction(id="r1", reactants=["A"], products=["B"])
        r2 = r.with_flux_bounds(-5.0, 10.0)
        bounds = r2.flux_bounds
        assert jnp.allclose(bounds, jnp.array([-5.0, 10.0]))


class TestReactionWithStoichiometry:
    def test_with_stoichiometry(self):
        r = Reaction(id="r1", reactants=["A"], products=["B"])
        r2 = r.with_stoichiometry({"A": -2.0, "B": 2.0})
        assert r2.stoichiometry_dict == {"A": -2.0, "B": 2.0}
