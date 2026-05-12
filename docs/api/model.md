# Model

The model module defines the metabolic network components.

::: fluxomics_data_converter.model.metabolite.Metabolite
    options:
      show_root_heading: true

::: fluxomics_data_converter.model.reaction.Reaction
    options:
      show_root_heading: true
      members:
        - is_variant_reaction
        - n_variants

::: fluxomics_data_converter.model.atom_mapping.AtomTransition
    options:
      show_root_heading: true
      members:
        - to_letter_notation
        - transform_isotopomers
        - weights

::: fluxomics_data_converter.model.constraint.Constraints
    options:
      show_root_heading: true

::: fluxomics_data_converter.model.constraint_eval.ConstraintEvaluator
    options:
      show_root_heading: true
      members:
        - parse_formula
        - evaluate_constraints
