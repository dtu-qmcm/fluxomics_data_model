"""FluxML constraint definitions.

The three concrete constraint classes (:class:`NetConstraints`,
:class:`ExchangeConstraints`, :class:`MetaboliteSizeConstraints`) are
structurally identical — each holds a list of
:class:`ConstraintFormula` objects and produces a human-readable
``summary()``.  To avoid repeating the same body three times, the shared
behaviour lives in the private :class:`_ConstraintsMixin`.

Hierarchy::

    _ConstraintsMixin          ← shared summary() logic
        NetConstraints
        ExchangeConstraints
        MetaboliteSizeConstraints
    Constraints                ← container for all three
"""

from typing import ClassVar, Optional, List, Tuple, Any
from pydantic import BaseModel, Field, ConfigDict
import sympy as sp


class ConstraintFormula(BaseModel):
    """A single constraint formula.

    Can be either a textual mathematical expression or MathML.
    Textual format examples:
        - Equality: "uptGLYC = 0.5154"
        - Inequality: "bmALA >= 0.75*0.22601*mu"
        - Upper bound: "ENO <= 10"
        - Named: "ratio: uptUGlyc = 0.12*uptGLYC"

    Constraints can reference:
        - Reaction IDs (net or exchange fluxes): Variables
          representing flux through reactions
        - Metabolite IDs (pool sizes): Variables representing
          metabolite concentrations
        - Parameters: Named constants (e.g., 'mu' for growth rate)
    """

    name: Optional[str] = Field(
        default=None,
        description="Optional name for the constraint (for named constraints)",
    )
    expression: str = Field(
        description="Constraint expression (textual formula or MathML)"
    )
    is_mathml: bool = Field(
        default=False, description="Whether the expression is in MathML format"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def parse_sympy(self) -> Tuple[Any, str, Any]:
        """Parse expression into SymPy components.

        Returns:
            Tuple of (lhs_expr, operator, rhs_expr) where lhs_expr and rhs_expr
            are SymPy expressions and operator is one of '>=', '<=', '='

        Raises:
            ValueError: If expression cannot be parsed or has no valid operator
            NotImplementedError: If expression is in MathML format
        """
        if self.is_mathml:
            raise NotImplementedError("MathML constraints not yet supported")

        expr = self.expression.strip()

        # Detect operator
        operator = None
        for op in [">=", "<=", "="]:
            if op in expr:
                operator = op
                lhs_str, rhs_str = expr.split(op, 1)
                lhs_str = lhs_str.strip()
                rhs_str = rhs_str.strip()
                break

        if operator is None:
            raise ValueError(f"No valid operator found in constraint: {expr}")

        # Parse into SymPy expressions
        try:
            lhs_expr = sp.sympify(lhs_str)
            rhs_expr = sp.sympify(rhs_str)
        except Exception as e:
            raise ValueError(f"Failed to parse constraint '{expr}': {e}")

        return lhs_expr, operator, rhs_expr

    def get_variable_names(self) -> set[str]:
        """Extract all variable names from the constraint expression.

        Returns:
            Set of variable names (reaction IDs, metabolite IDs, or parameters)
        """
        if self.is_mathml:
            return set()  # Would need MathML parsing

        try:
            lhs_expr, _, rhs_expr = self.parse_sympy()
            lhs_vars = {str(s) for s in lhs_expr.free_symbols}
            rhs_vars = {str(s) for s in rhs_expr.free_symbols}
            return lhs_vars | rhs_vars
        except Exception:
            return set()

    def __str__(self) -> str:
        """String representation of the constraint."""
        if self.name:
            return f"{self.name}: {self.expression}"
        return self.expression


class _ConstraintsMixin:
    """Shared ``summary()`` implementation for single-type constraint classes.

    Each concrete subclass holds a ``formulas`` list and a ``_label`` class
    variable that names the constraint type (e.g. ``"Net"``).  This mixin
    provides a single, tested implementation of ``summary()`` so that
    :class:`NetConstraints`, :class:`ExchangeConstraints`, and
    :class:`MetaboliteSizeConstraints` do not duplicate the same body.
    """

    # Subclasses must declare these at class level:
    #   formulas: List[ConstraintFormula]
    #   _label:   ClassVar[str]

    _label: ClassVar[str] = "Constraint"  # overridden in each subclass

    def summary(self) -> str:
        """Return a human-readable summary of this constraint group.

        Counts equalities (``=``), lower-bound inequalities (``>=``), and
        upper-bound inequalities (``<=``), then lists up to five sample
        formulas.

        Returns:
            Multi-line string with counts and sample formula expressions.
        """
        if not self.formulas:  # type: ignore[attr-defined]
            return f"No {self._label.lower()} constraints"

        formulas = self.formulas  # type: ignore[attr-defined]
        lines = [f"{self._label} Constraints ({len(formulas)} formulas):"]

        equalities = sum(
            1
            for f in formulas
            if "=" in f.expression
            and ">=" not in f.expression
            and "<=" not in f.expression
        )
        inequalities_gte = sum(1 for f in formulas if ">=" in f.expression)
        inequalities_lte = sum(1 for f in formulas if "<=" in f.expression)

        lines.append(f"  Equalities (=): {equalities}")
        lines.append(f"  Lower bounds (>=): {inequalities_gte}")
        lines.append(f"  Upper bounds (<=): {inequalities_lte}")

        if len(formulas) <= 5:
            lines.append("  Formulas:")
            for f in formulas:
                lines.append(f"    - {str(f)}")
        else:
            lines.append("  Sample formulas:")
            for f in formulas[:3]:
                lines.append(f"    - {str(f)}")
            lines.append(f"    ... and {len(formulas) - 3} more")

        return "\n".join(lines)


class NetConstraints(_ConstraintsMixin, BaseModel):
    """Net flux constraints (``fluxml/constraints/net``).

    Apply to the net flux through each reaction, defined as
    ``forward_flux - reverse_flux``.  These typically fix or bound the
    net throughput of specific reactions in the model.
    """

    _label: ClassVar[str] = "Net"

    formulas: List[ConstraintFormula] = Field(
        default_factory=list, description="List of net flux constraint formulas"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")


class ExchangeConstraints(_ConstraintsMixin, BaseModel):
    """Exchange flux constraints (``fluxml/constraints/xch``).

    Apply to exchange fluxes — the individual forward or reverse component
    magnitudes of bidirectional reactions.  Exchange constraints allow
    separate control of the back-flux independently of the net flux.
    """

    _label: ClassVar[str] = "Exchange"

    formulas: List[ConstraintFormula] = Field(
        default_factory=list,
        description="List of exchange flux constraint formulas",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")


class MetaboliteSizeConstraints(_ConstraintsMixin, BaseModel):
    """Metabolite pool-size constraints (``fluxml/constraints/psize``).

    Apply to metabolite concentrations rather than reaction fluxes.
    Pool-size constraints are primarily used in non-stationary (INST) 13C
    MFA, where the pool sizes are free variables alongside the fluxes.
    """

    _label: ClassVar[str] = "Metabolite Size"

    formulas: List[ConstraintFormula] = Field(
        default_factory=list,
        description="List of metabolite size constraint formulas",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")


class Constraints(BaseModel):
    """FluxML constraints collection.

    Groups all constraint types:
    - Net constraints: Apply to net fluxes (forward - reverse) through reactions
    - Exchange constraints: Apply to individual forward/reverse exchange fluxes
    - Metabolite size constraints: Apply to metabolite pool sizes/concentrations

    The key difference is what variables they constrain:
    - Net/Exchange: Constrain reaction IDs (flux variables)
    - Metabolite size: Constrain metabolite IDs (pool size variables)

    Corresponds to fluxml/constraints
    """

    net: Optional[NetConstraints] = Field(
        default=None, description="Net flux constraints"
    )
    xch: Optional[ExchangeConstraints] = Field(
        default=None, description="Exchange flux constraints"
    )
    metabolitesize: Optional[MetaboliteSizeConstraints] = Field(
        default=None, description="Metabolite size constraints"
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def summary(self) -> str:
        """Generate a comprehensive summary of all constraints.

        Returns:
            Multi-line string with constraint statistics and samples
        """
        sections = []

        if self.net:
            sections.append(self.net.summary())

        if self.xch:
            sections.append(self.xch.summary())

        if self.metabolitesize:
            sections.append(self.metabolitesize.summary())

        if not sections:
            return "No constraints defined"

        # Add overall statistics at the top
        total_formulas = sum(
            [
                len(self.net.formulas) if self.net else 0,
                len(self.xch.formulas) if self.xch else 0,
                len(self.metabolitesize.formulas) if self.metabolitesize else 0,
            ]
        )

        header = f"=== Constraints ===\n  Total formulas: {total_formulas}\n\n"
        return header + "\n\n".join(sections)

    def __repr__(self) -> str:
        """String representation using summary."""
        return self.summary()

    def get_all_variable_names(self) -> set[str]:
        """Extract all unique variable names referenced in all constraints.

        Returns:
            Set of all variable names (reaction IDs, metabolite IDs, parameters)
        """
        all_vars = set()

        if self.net:
            for formula in self.net.formulas:
                all_vars.update(formula.get_variable_names())

        if self.xch:
            for formula in self.xch.formulas:
                all_vars.update(formula.get_variable_names())

        if self.metabolitesize:
            for formula in self.metabolitesize.formulas:
                all_vars.update(formula.get_variable_names())

        return all_vars
