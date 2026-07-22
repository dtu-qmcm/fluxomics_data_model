"""JAX-compatible constraint evaluation bridging SymPy and JAX.

:class:`ConstraintEvaluator` converts
:class:`~fluxomics_data_converter.model.constraint.ConstraintFormula` objects —
which store constraints as plain strings such as ``"bmALA >= 0.75*mu"`` —
into callable JAX functions that can be embedded in gradient-based
optimisation loops.

The bridge works as follows:

1. **Parsing** — :meth:`ConstraintFormula.parse_sympy` tokenises the
   expression and returns ``(lhs_expr, operator, rhs_expr)`` as SymPy
   objects.
2. **Parameter substitution** — known scalar parameters (e.g.
   ``{"mu": 0.03}``) are substituted into both sides, so the RHS reduces
   to a numeric constant.
3. **JAX compilation** — :func:`sympy.lambdify` converts the LHS expression
   into a Python callable with ``modules="jax"``, producing a JAX-traceable
   function.
4. **Index mapping** — the evaluator resolves each SymPy symbol to its
   position in the flux vector using the ``reaction_ids`` list supplied at
   construction time, then wraps the lambdified function in a closure that
   extracts the relevant elements.

Limitations
-----------
- The RHS must be numeric after parameter substitution.  Constraints of the
  form ``flux_A = flux_B`` (both sides symbolic) are not supported; rewrite
  them as ``flux_A - flux_B = 0``.
- MathML constraints raise ``NotImplementedError``.
"""

from typing import Dict, Callable, Tuple, List, Optional
import sympy as sp
import jax.numpy as jnp
from jax import Array

from .constraint import ConstraintFormula


class ConstraintEvaluator:
    """Evaluates constraint formulas as JAX-compatible functions.

    Parses textual constraint formulas using SymPy and converts them
    to JAX functions that can be used in optimization.
    """

    def __init__(
        self,
        reaction_ids: List[str],
        parameters: Optional[Dict[str, float]] = None,
    ):
        """Initialize constraint evaluator.

        Args:
            reaction_ids: List of reaction IDs (order determines flux
                vector indexing)
            parameters: Dictionary of parameter names and values
                (e.g., {'mu': 0.03})
        """
        self.reaction_ids = reaction_ids
        self.parameters = parameters or {}

        # Create mapping from reaction ID to flux vector index
        self.flux_map = {rid: idx for idx, rid in enumerate(reaction_ids)}

    def parse_formula(
        self, formula: ConstraintFormula
    ) -> Tuple[Callable[[Array], Array], str, float]:
        """Parse a constraint formula into a JAX-compatible function.

        The constraint is converted to the form: lhs - rhs [op] 0
        For inequality constraints (>= or <=), the returned function
        evaluates the residual that should be >= 0.
        For equality constraints (=), the returned function evaluates
        the residual that should be = 0.

        This function leverages ConstraintFormula.parse_sympy() for parsing,
        then adds:
        - Parameter substitution
        - Flux vector indexing
        - JAX function conversion

        Args:
            formula: ConstraintFormula to parse

        Returns:
            Tuple of (constraint_fn, operator, rhs_value) where:
            - constraint_fn: Function that takes flux vector and
              returns residual
            - operator: One of '>=', '<=', '='
            - rhs_value: Right-hand side value after parameter substitution

        Example:
            For "bmALA >= 0.75*0.22601*mu" with mu=0.03:
            - Returns (fn, '>=', 0.005085225) where fn(v) =
              v[idx_bmALA] - 0.005085225
            - In optimization: fn(v) >= 0 means bmALA >= 0.005085225
        """
        # Use ConstraintFormula's built-in SymPy parser
        try:
            lhs_expr, operator, rhs_expr = formula.parse_sympy()
        except NotImplementedError:
            # Re-raise NotImplementedError as-is (e.g., for MathML)
            raise
        except Exception as e:
            raise ValueError(f"Failed to parse constraint formula: {e}")

        # Substitute parameter values
        param_subs = {sp.Symbol(p): v for p, v in self.parameters.items()}
        lhs_expr = lhs_expr.subs(param_subs)
        rhs_expr = rhs_expr.subs(param_subs)

        # Check if RHS is numeric after substitution
        if not rhs_expr.is_number:
            raise ValueError(
                f"RHS is not numeric after parameter substitution: {rhs_expr}"
            )

        rhs_value = float(rhs_expr)

        # Create constraint function
        # For all operators, we compute: lhs - rhs
        # - For '>=': lhs - rhs >= 0
        # - For '<=': lhs - rhs <= 0 → -(lhs - rhs) >= 0 → rhs - lhs >= 0
        # - For '=': lhs - rhs = 0

        # Get free symbols in LHS (should be reaction IDs)
        free_symbols = lhs_expr.free_symbols

        if len(free_symbols) == 0:
            # Constant constraint (e.g., "1 >= 0")
            if operator == ">=":
                residual = float(lhs_expr) - rhs_value
            elif operator == "<=":
                residual = rhs_value - float(lhs_expr)
            else:  # '='
                residual = float(lhs_expr) - rhs_value

            def const_fn(flux_vector: Array) -> Array:
                return jnp.array(residual)

            return const_fn, operator, rhs_value

        # Build ordered list of symbols matching reaction order
        # Match free symbols to reaction IDs
        symbol_order = []
        for rid in self.reaction_ids:
            sym = sp.Symbol(rid)
            if sym in free_symbols:
                symbol_order.append(sym)

        # Convert to JAX function using lambdify
        if operator == ">=":
            # lhs >= rhs → lhs - rhs >= 0
            residual_expr = lhs_expr - rhs_value
        elif operator == "<=":
            # lhs <= rhs → rhs - lhs >= 0
            residual_expr = rhs_value - lhs_expr
        else:  # '='
            # lhs = rhs → lhs - rhs = 0
            residual_expr = lhs_expr - rhs_value

        # Create JAX function
        jax_fn = sp.lambdify(symbol_order, residual_expr, modules="jax")

        # Wrapper to extract relevant fluxes from flux vector
        def constraint_fn(flux_vector: Array) -> Array:
            # Extract flux values for symbols in the expression
            flux_values = [
                flux_vector[self.flux_map[str(sym)]] for sym in symbol_order
            ]
            return jax_fn(*flux_values) if flux_values else jax_fn()

        return constraint_fn, operator, rhs_value

    def parse_all(
        self, formulas: List[ConstraintFormula]
    ) -> List[Tuple[Callable[[Array], Array], str, float, Optional[str]]]:
        """Parse all constraint formulas.

        Args:
            formulas: List of ConstraintFormula objects

        Returns:
            List of tuples (constraint_fn, operator, rhs_value, name)
        """
        results = []
        for formula in formulas:
            constraint_fn, operator, rhs_value = self.parse_formula(formula)
            results.append((constraint_fn, operator, rhs_value, formula.name))
        return results

    def evaluate_constraints(
        self, flux_vector: Array, formulas: List[ConstraintFormula]
    ) -> Dict[str, float]:
        """Evaluate all constraints for a given flux vector.

        Args:
            flux_vector: JAX array of flux values
            formulas: List of constraint formulas

        Returns:
            Dictionary mapping constraint description to residual value
        """
        results = {}
        for i, formula in enumerate(formulas):
            constraint_fn, operator, rhs_value = self.parse_formula(formula)
            residual = float(constraint_fn(flux_vector))

            name = formula.name or f"constraint_{i + 1}"
            results[f"{name} ({formula.expression})"] = residual

        return results
