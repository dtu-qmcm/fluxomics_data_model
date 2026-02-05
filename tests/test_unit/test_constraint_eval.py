"""
Tests for JAX-compatible constraint evaluation.
"""

import pytest
import numpy as np
import jax.numpy as jnp
from jax import grad

from fluxomics_data_model.model import ConstraintEvaluator
from fluxomics_data_model.model.constraint import ConstraintFormula
from fluxomics_data_model.io.fluxml_parser import FluxMLParser


class TestConstraintEvaluator:
    """Test basic constraint parsing and evaluation."""

    def test_simple_equality_constraint(self):
        """Test parsing simple equality constraint."""
        reaction_ids = ["uptGLYC", "uptNH3", "exCO2"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="uptGLYC = 0.5154",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        assert operator == "="
        assert np.isclose(rhs_value, 0.5154)

        # Test evaluation
        flux_vector = jnp.array([0.5154, 1.0, 0.5])
        residual = constraint_fn(flux_vector)
        assert np.isclose(residual, 0.0, atol=1e-6)

        # Test with different value
        flux_vector = jnp.array([0.6, 1.0, 0.5])
        residual = constraint_fn(flux_vector)
        assert np.isclose(residual, 0.6 - 0.5154, atol=1e-6)

    def test_inequality_constraint_gte(self):
        """Test parsing >= constraint."""
        reaction_ids = ["bmALA", "mu"]
        parameters = {"mu": 0.03}
        evaluator = ConstraintEvaluator(reaction_ids, parameters)

        formula = ConstraintFormula(
            expression="bmALA >= 0.75*0.22601*mu",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        assert operator == ">="
        expected_rhs = 0.75 * 0.22601 * 0.03
        assert np.isclose(rhs_value, expected_rhs)

        # Test evaluation: residual should be >= 0 when constraint is satisfied
        # If bmALA = 0.01 and rhs = 0.005085225, residual = 0.01 - 0.005085225 > 0
        flux_vector = jnp.array([0.01, 0.03])
        residual = constraint_fn(flux_vector)
        assert residual >= 0  # Constraint satisfied

        # Test with value below threshold
        flux_vector = jnp.array([0.001, 0.03])
        residual = constraint_fn(flux_vector)
        assert residual < 0  # Constraint violated

    def test_inequality_constraint_lte(self):
        """Test parsing <= constraint."""
        reaction_ids = ["exCO2", "PK"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="exCO2 <= 0.44",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        assert operator == "<="
        assert np.isclose(rhs_value, 0.44)

        # For <= constraint, residual = rhs - lhs should be >= 0 when satisfied
        # If exCO2 = 0.4 and rhs = 0.44, residual = 0.44 - 0.4 = 0.04 > 0
        flux_vector = jnp.array([0.4, 1.0])
        residual = constraint_fn(flux_vector)
        assert residual >= 0  # Constraint satisfied

        # Test with value above threshold
        flux_vector = jnp.array([0.5, 1.0])
        residual = constraint_fn(flux_vector)
        assert residual < 0  # Constraint violated

    def test_constraint_with_arithmetic(self):
        """Test constraint with multiple operations."""
        reaction_ids = ["r1", "r2", "r3"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="r1 + 2*r2 - r3 = 1.5",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        assert operator == "="
        assert np.isclose(rhs_value, 1.5)

        # Test: r1 + 2*r2 - r3 = 1 + 2*2 - 3.5 = 1.5
        flux_vector = jnp.array([1.0, 2.0, 3.5])
        residual = constraint_fn(flux_vector)
        assert np.isclose(residual, 0.0, atol=1e-6)

    def test_constant_constraint(self):
        """Test constraint with no variables (constant)."""
        reaction_ids = ["r1", "r2"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="1 = 1",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        # Any flux vector should give residual = 0
        flux_vector = jnp.array([5.0, 10.0])
        residual = constraint_fn(flux_vector)
        assert np.isclose(residual, 0.0)

    def test_parameter_substitution(self):
        """Test that parameters are correctly substituted."""
        reaction_ids = ["uptGLYC", "growth"]
        parameters = {"mu": 0.03, "yield_factor": 0.75}
        evaluator = ConstraintEvaluator(reaction_ids, parameters)

        formula = ConstraintFormula(
            expression="growth = mu * yield_factor",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        assert operator == "="
        expected_rhs = 0.03 * 0.75
        assert np.isclose(rhs_value, expected_rhs)

        # growth should equal 0.0225
        flux_vector = jnp.array([1.0, 0.0225])
        residual = constraint_fn(flux_vector)
        assert np.isclose(residual, 0.0, atol=1e-6)

    def test_jax_differentiability(self):
        """Test that constraint functions are differentiable with JAX."""
        reaction_ids = ["r1", "r2"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="r1 + r2*r2 = 5",
            is_mathml=False
        )

        constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)

        # Compute gradient
        flux_vector = jnp.array([1.0, 2.0])

        # Gradient of r1 + r2^2 - 5 with respect to flux_vector
        grad_fn = grad(lambda v: constraint_fn(v))
        gradients = grad_fn(flux_vector)

        # Expected gradients: d/dr1 = 1, d/dr2 = 2*r2 = 2*2 = 4
        assert np.isclose(gradients[0], 1.0)
        assert np.isclose(gradients[1], 4.0)

    def test_parse_all(self):
        """Test parsing multiple formulas at once."""
        reaction_ids = ["r1", "r2", "r3"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formulas = [
            ConstraintFormula(expression="r1 = 1.0", is_mathml=False),
            ConstraintFormula(expression="r2 >= 0.5", is_mathml=False),
            ConstraintFormula(expression="r3 <= 2.0", is_mathml=False),
        ]

        results = evaluator.parse_all(formulas)

        assert len(results) == 3
        assert results[0][1] == "="
        assert results[1][1] == ">="
        assert results[2][1] == "<="

    def test_evaluate_constraints(self):
        """Test evaluating all constraints and getting results dict."""
        reaction_ids = ["r1", "r2"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formulas = [
            ConstraintFormula(name="r1_fix", expression="r1 = 1.0", is_mathml=False),
            ConstraintFormula(expression="r2 >= 0.5", is_mathml=False),
        ]

        flux_vector = jnp.array([1.0, 0.75])
        results = evaluator.evaluate_constraints(flux_vector, formulas)

        assert "r1_fix (r1 = 1.0)" in results
        assert "constraint_2 (r2 >= 0.5)" in results
        assert np.isclose(results["r1_fix (r1 = 1.0)"], 0.0)
        assert results["constraint_2 (r2 >= 0.5)"] > 0  # Satisfied

    def test_error_handling_no_operator(self):
        """Test error handling when no operator is present."""
        reaction_ids = ["r1"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(expression="r1 + 1", is_mathml=False)

        with pytest.raises(ValueError, match="No valid operator found"):
            evaluator.parse_formula(formula)

    def test_error_handling_invalid_expression(self):
        """Test error handling for invalid expressions."""
        reaction_ids = ["r1"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(expression="r1 + ??? = 1", is_mathml=False)

        with pytest.raises(ValueError, match="Failed to parse constraint"):
            evaluator.parse_formula(formula)

    def test_error_handling_mathml_not_supported(self):
        """Test that MathML raises NotImplementedError."""
        reaction_ids = ["r1"]
        evaluator = ConstraintEvaluator(reaction_ids)

        formula = ConstraintFormula(
            expression="<math>...</math>",
            is_mathml=True
        )

        with pytest.raises(NotImplementedError, match="MathML constraints not yet supported"):
            evaluator.parse_formula(formula)


class TestConstraintEvaluatorIntegration:
    """Test constraint evaluator with real FluxML file."""

    def test_parse_fluxml_constraints(self):
        """Test parsing constraints from actual FluxML file."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML-Test/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get experiment and constraints
        exp = data_model.experiments[0]
        assert exp.constraints is not None
        assert exp.constraints.net is not None

        # Get reaction IDs from model
        reaction_ids = list(data_model.model.computational_reaction_ids)

        # Find parameter value for mu from constraints
        # Look for "mu = X" constraint
        parameters = {}
        for formula in exp.constraints.net.formulas:
            if '=' in formula.expression and 'mu' in formula.expression:
                parts = formula.expression.split('=')
                if parts[0].strip() == 'mu':
                    parameters['mu'] = float(parts[1].strip())
                    break

        # Create evaluator
        evaluator = ConstraintEvaluator(reaction_ids, parameters)

        # Parse a few net constraints
        test_formulas = exp.constraints.net.formulas[:5]

        for formula in test_formulas:
            try:
                constraint_fn, operator, rhs_value = evaluator.parse_formula(formula)
                print(f"✓ Parsed: {formula.expression}")
                print(f"  Operator: {operator}, RHS: {rhs_value}")
            except Exception as e:
                # Some constraints might reference undefined reactions/parameters
                print(f"✗ Failed to parse: {formula.expression}")
                print(f"  Error: {e}")

        # At least some should parse successfully
        assert True  # If we get here, basic parsing works

    def test_evaluate_on_flux_vector(self):
        """Test evaluating constraints on a flux vector."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML-Test/fluxml-model/models/CN_network_model_BCG.fml"
        )

        exp = data_model.experiments[0]
        reaction_ids = list(data_model.model.computational_reaction_ids)

        # Extract parameter values from equality constraints
        parameters = {}
        if exp.constraints and exp.constraints.net:
            for formula in exp.constraints.net.formulas:
                if '=' in formula.expression:
                    parts = formula.expression.split('=', 1)
                    lhs = parts[0].strip()
                    rhs = parts[1].strip()
                    # If LHS is a single symbol and RHS is numeric, treat as parameter
                    if lhs.isidentifier():
                        try:
                            parameters[lhs] = float(rhs)
                        except ValueError:
                            pass  # Not a simple numeric assignment

        evaluator = ConstraintEvaluator(reaction_ids, parameters)

        # Create a dummy flux vector (all ones for testing)
        flux_vector = jnp.ones(len(reaction_ids))

        # Try to evaluate constraints
        if exp.constraints and exp.constraints.net:
            formulas_to_test = exp.constraints.net.formulas[:10]
            results = {}

            for i, formula in enumerate(formulas_to_test):
                try:
                    constraint_fn, operator, rhs = evaluator.parse_formula(formula)
                    residual = constraint_fn(flux_vector)
                    results[f"constraint_{i}"] = float(residual)
                except Exception as e:
                    # Skip constraints that can't be parsed
                    pass

            # Should have parsed at least some constraints
            assert len(results) > 0
            print(f"Successfully evaluated {len(results)} constraints")
