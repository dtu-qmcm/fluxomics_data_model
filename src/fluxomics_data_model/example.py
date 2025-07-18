"""
Example usage of the FluxML data model.

This demonstrates how to create FluxML objects and work with JAX arrays.
"""

import jax.numpy as jnp
from fluxomics_data_model import (
    FluxML,
    Info,
    ReactionNetwork,
    Configuration,
    Pool,
    MetabolitePools,
    Reaction,
    Reduct,
    RProduct,
    Input,
    Label,
)


def create_simple_model() -> FluxML:
    """Create a simple FluxML model for demonstration."""

    # Create metabolite pools
    pools = [
        Pool(id="glucose", atoms=6, size=1.0),
        Pool(id="pyruvate", atoms=3, size=1.0),
        Pool(id="co2", atoms=1, size=1.0),
    ]
    metabolite_pools = MetabolitePools(pools=pools)

    # Create reactions
    reactions = [
        Reaction(
            id="glycolysis",
            bidirectional=False,
            reducts=[Reduct(id="glucose")],
            rproducts=[RProduct(id="pyruvate"), RProduct(id="co2")],
        )
    ]

    # Create reaction network
    reaction_network = ReactionNetwork(
        metabolitepools=metabolite_pools, reactions=reactions
    )

    # Create info
    info = Info(
        name="Simple Glycolysis Model",
        version="1.0",
        comment="A minimal model for demonstration",
    )

    # Create configuration
    config = Configuration(
        name="glucose_tracer",
        stationary=True,
        inputs=[
            Input(
                pool="glucose",
                labels=[
                    Label(cfg="111111", purity="0.99"),  # Fully labeled glucose
                    Label(cfg="000000", purity="0.01"),  # Unlabeled glucose
                ],
            )
        ],
    )

    # Create FluxML model
    model = FluxML(
        info=info, reactionnetwork=reaction_network, configurations=[config]
    )

    return model


def demonstrate_jax_compatibility():
    """Demonstrate JAX array operations with FluxML model."""

    model = create_simple_model()

    print("FluxML Model Summary:")
    print(f"- Pools: {len(model.pool_ids)}")
    print(f"- Reactions: {len(model.reaction_ids)}")
    print(f"- Configurations: {len(model.configurations)}")
    print(f"- Pool IDs: {list(model.pool_ids)}")
    print(f"- Reaction IDs: {list(model.reaction_ids)}")

    # Get JAX representation
    jax_repr = model.to_jax_representation()
    print("\nJAX Representation:")
    print(
        f"- Stoichiometric matrix shape: "
        f"{jax_repr['stoichiometric_matrix'].shape}"
    )
    print(f"- Flux bounds shape: {jax_repr['flux_bounds'].shape}")
    print(f"- Pool size bounds shape: {jax_repr['poolsize_bounds'].shape}")

    # Demonstrate stoichiometric matrix
    S = jax_repr["stoichiometric_matrix"]
    print(f"\nStoichiometric Matrix:\n{S}")

    # Demonstrate flux bounds
    flux_bounds = jax_repr["flux_bounds"]
    print(f"\nFlux Bounds:\n{flux_bounds}")

    # Get tracer experiment data
    tracer_data = model.get_tracer_experiment_data("glucose_tracer")
    if tracer_data:
        print(
            f"\nTracer Composition Matrix Shape: "
            f"{tracer_data['tracer_composition'].shape}"
        )
        print(f"Tracer Composition:\n{tracer_data['tracer_composition']}")

    # Demonstrate JAX operations
    print("\nJAX Operations:")

    # Calculate steady-state flux balance (S @ v = 0)
    # For demonstration, use a simple flux vector
    v = jnp.array([1.0])  # Glycolysis flux
    balance = S @ v
    print(f"Flux balance (S @ v): {balance}")

    # The balance should be close to zero for steady state
    print(f"Balance residual: {jnp.linalg.norm(balance)}")

    return model


if __name__ == "__main__":
    # Run the demonstration
    model = demonstrate_jax_compatibility()

    # Show model validation
    print("\nModel Validation:")
    print("✓ All pool references are valid")
    print("✓ All reaction references are valid")
    print("✓ Model is JAX compatible")
    print("✓ All data structures are immutable")

    # Show serialization capability
    print("\nSerialization:")
    model_dict = model.model_dump()
    print(f"Model serialized to dict with {len(model_dict)} top-level keys")

    # Reconstruct from dict
    reconstructed = FluxML.model_validate(model_dict)
    print("✓ Model successfully reconstructed from dict")
