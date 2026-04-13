import pytest
import jax.numpy as jnp

from fluxomics_data_model.model.metabolite import Metabolite
from fluxomics_data_model.core.common import Annotation


class TestMetaboliteCreation:
    def test_basic_creation(self):
        m = Metabolite(id="g6p", name="glucose 6-phosphate", atoms=6)
        assert m.id == "g6p"
        assert m.name == "glucose 6-phosphate"
        assert m.atoms == 6

    def test_minimal_creation(self):
        m = Metabolite(id="pyr")
        assert m.id == "pyr"
        assert m.name is None
        assert m.atoms == 0

    def test_all_fields(self):
        m = Metabolite(
            id="cit",
            name="citrate",
            atoms=6,
            weight=192.12,
            charge=-3,
            formula="C6H5O7",
            compartment="c",
        )
        assert m.weight == 192.12
        assert m.charge == -3
        assert m.formula == "C6H5O7"
        assert m.compartment == "c"

    def test_frozen(self):
        m = Metabolite(id="g6p")
        with pytest.raises(Exception):
            m.id = "changed"


class TestMetaboliteFormulaValidation:
    def test_valid_formula(self):
        m = Metabolite(id="glc", formula="C6H12O6")
        assert m.formula == "C6H12O6"

    def test_invalid_formula(self):
        with pytest.raises(Exception):
            Metabolite(id="bad", formula="not_a_formula")

    def test_empty_formula(self):
        m = Metabolite(id="x", formula="")
        assert m.formula == ""

    def test_zero_formula(self):
        m = Metabolite(id="x", formula="0")
        assert m.formula == "0"

    def test_none_formula(self):
        m = Metabolite(id="x")
        assert m.formula is None

    def test_formula_with_subscripts(self):
        m = Metabolite(id="atp", formula="C10H16N5O13P3")
        assert m.formula == "C10H16N5O13P3"


class TestMetaboliteMolecularWeight:
    def test_computed_from_formula(self):
        m = Metabolite(id="glc", formula="C6H12O6")
        mw = m.molecular_weight
        assert mw is not None
        assert abs(mw - 180.156) < 1.0

    def test_no_formula_no_weight(self):
        m = Metabolite(id="x")
        assert m.molecular_weight is None

    def test_empty_formula_no_weight(self):
        m = Metabolite(id="x", formula="")
        assert m.molecular_weight is None


class TestMetaboliteAtomVector:
    def test_atom_vector_with_atoms(self):
        m = Metabolite(id="glc", atoms=6)
        vec = m.atom_vector
        assert vec.shape == (6,)
        assert jnp.all(vec == 1.0)

    def test_atom_vector_zero_atoms(self):
        m = Metabolite(id="x", atoms=0)
        vec = m.atom_vector
        assert vec.shape == (0,)

    def test_with_jax_atoms(self):
        m = Metabolite(id="glc", atoms=3)
        new_atoms = jnp.array([1.0, 2.0, 3.0])
        m2 = m.with_jax_atoms(new_atoms)
        assert m2.id == "glc"
        result = m2.atom_vector
        assert jnp.allclose(result, new_atoms)


class TestMetaboliteAnnotations:
    def test_with_annotations(self):
        m = Metabolite(
            id="glc",
            annotations=[Annotation(name="inchi", content="InChI=1S/C6H12O6")],
        )
        assert len(m.annotations) == 1
        assert m.annotations[0].name == "inchi"

    def test_default_empty_annotations(self):
        m = Metabolite(id="glc")
        assert m.annotations == []
