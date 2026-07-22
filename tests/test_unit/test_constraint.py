import pytest

from fluxomics_data_converter.model.constraint import (
    ConstraintFormula,
    NetConstraints,
    ExchangeConstraints,
    MetaboliteSizeConstraints,
    Constraints,
)


class TestConstraintFormula:
    def test_textual_equality(self):
        f = ConstraintFormula(expression="r1 = 1.0")
        assert f.expression == "r1 = 1.0"
        assert f.is_mathml is False

    def test_named_constraint(self):
        f = ConstraintFormula(name="fix_r1", expression="r1 = 1.0")
        assert f.name == "fix_r1"

    def test_parse_sympy_equality(self):
        f = ConstraintFormula(expression="r1 = 1.0")
        lhs, op, rhs = f.parse_sympy()
        assert op == "="
        assert float(rhs) == 1.0

    def test_parse_sympy_gte(self):
        f = ConstraintFormula(expression="r1 >= 0.5")
        lhs, op, rhs = f.parse_sympy()
        assert op == ">="

    def test_parse_sympy_lte(self):
        f = ConstraintFormula(expression="r1 <= 10")
        lhs, op, rhs = f.parse_sympy()
        assert op == "<="

    def test_parse_sympy_no_operator_raises(self):
        f = ConstraintFormula(expression="r1 + 1")
        with pytest.raises(ValueError, match="No valid operator"):
            f.parse_sympy()

    def test_mathml_not_supported(self):
        f = ConstraintFormula(expression="<math>...</math>", is_mathml=True)
        with pytest.raises(NotImplementedError, match="MathML"):
            f.parse_sympy()

    def test_get_variable_names(self):
        f = ConstraintFormula(expression="r1 + r2 = 3.0")
        names = f.get_variable_names()
        assert names == {"r1", "r2"}

    def test_str_with_name(self):
        f = ConstraintFormula(name="fix", expression="r1 = 1.0")
        assert str(f) == "fix: r1 = 1.0"

    def test_str_without_name(self):
        f = ConstraintFormula(expression="r1 = 1.0")
        assert str(f) == "r1 = 1.0"


class TestNetConstraints:
    def test_empty(self):
        nc = NetConstraints()
        assert nc.formulas == []

    def test_with_formulas(self):
        nc = NetConstraints(
            formulas=[
                ConstraintFormula(expression="r1 = 1.0"),
                ConstraintFormula(expression="r2 >= 0.5"),
            ]
        )
        assert len(nc.formulas) == 2

    def test_summary(self):
        nc = NetConstraints(
            formulas=[
                ConstraintFormula(expression="r1 = 1.0"),
            ]
        )
        s = nc.summary()
        assert "Net Constraints" in s


class TestExchangeConstraints:
    def test_empty(self):
        ec = ExchangeConstraints()
        assert ec.formulas == []


class TestMetaboliteSizeConstraints:
    def test_empty(self):
        msc = MetaboliteSizeConstraints()
        assert msc.formulas == []


class TestConstraints:
    def test_empty(self):
        c = Constraints()
        assert c.net is None
        assert c.xch is None
        assert c.metabolitesize is None

    def test_with_net(self):
        c = Constraints(
            net=NetConstraints(
                formulas=[
                    ConstraintFormula(expression="r1 = 1.0"),
                ]
            )
        )
        assert c.net is not None
        assert len(c.net.formulas) == 1

    def test_get_all_variable_names(self):
        c = Constraints(
            net=NetConstraints(
                formulas=[
                    ConstraintFormula(expression="r1 = 1.0"),
                ]
            ),
            xch=ExchangeConstraints(
                formulas=[
                    ConstraintFormula(expression="r2 >= 0"),
                ]
            ),
        )
        names = c.get_all_variable_names()
        assert "r1" in names
        assert "r2" in names

    def test_summary(self):
        c = Constraints(
            net=NetConstraints(
                formulas=[
                    ConstraintFormula(expression="r1 = 1.0"),
                ]
            )
        )
        s = c.summary()
        assert "Constraints" in s
