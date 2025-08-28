"""
Atom mapping functionality for fluxomics reactions.

This module provides classes and utilities for handling atom mappings in
metabolic reactions, including support for:
- Multiple mapping formats (letter notation, FluxML cfg strings, RDKit SMILES)
- Symmetric reactions with multiple equivalent mappings
- Isotopomer transformations
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pydantic import BaseModel, Field, ConfigDict
import re
import numpy as np


@dataclass(frozen=True)
class AtomAddress:
    """Address for a specific atom in a reaction.

    Attributes:
        mol: Molecule/compound ID
        index: Atom index within the molecule (1-based)
        element: Element type (e.g., 'C', 'N', 'O')
        instance: Instance number when multiple same compounds exist (1-based)
    """

    mol: str
    index: int
    element: Optional[str] = None
    instance: int = 1


class AtomMap(BaseModel):
    """A single atom mapping for a reaction.

    Maps product atoms to their source reactant atoms.
    Keys and values are AtomAddress objects with 1-based indices.
    """

    # We need to use tuple representation for Pydantic compatibility
    mapping: Dict[AtomAddress, AtomAddress] = Field(
        default_factory=dict,
        description="Mapping from product atoms to source reactant atoms",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    def add(
        self,
        product_cpd: str,
        product_atom: int,
        src_cpd: str,
        src_atom: int,
        atom_type: Optional[str] = "C",
        product_instance: int = 1,
        src_instance: int = 1,
    ) -> "AtomMap":
        """Add a mapping from product atom to source atom."""
        # Create new dict to maintain immutability
        new_mapping = dict(self.mapping)
        new_mapping[
            AtomAddress(product_cpd, product_atom, atom_type, product_instance)
        ] = AtomAddress(src_cpd, src_atom, atom_type, src_instance)
        return AtomMap(mapping=new_mapping)

    def product_atoms(self) -> AtomAddress:
        """Return number of atoms per product compound instance."""
        result: Dict[Tuple[str, int], int] = {}
        for cpd, atom, _, instance in self.mapping.keys():
            key = (cpd, instance)
            result[key] = max(result.get(key, 0), atom)
        return AtomAddress(result)

    def reactant_atoms(self) -> AtomAddress:
        """Return number of atoms per reactant compound instance."""
        result: Dict[Tuple[str, int], int] = {}
        for cpd, atom, _, instance in self.mapping.values():
            key = (cpd, instance)
            result[key] = max(result.get(key, 0), atom)
        return AtomAddress(result)


class AtomMapping(BaseModel):
    """
    Complete atom mapping information for a reaction.

    Supports multiple alternative mappings (variants/symmetries).
    """

    maps: List[AtomMap] = Field(
        default_factory=lambda: [AtomMap()],
        description="List of alternative atom mappings",
    )
    weights: Optional[List[float]] = Field(
        default=None,
        description="Weights for each mapping variant (defaults to uniform)",
    )

    model_config = ConfigDict(frozen=True, extra="forbid")

    @classmethod
    def from_letter_notation(
        cls,
        reactant_items: List[Tuple[str, str]],
        product_items: List[Tuple[str, str]],
    ) -> "AtomMapping":
        """Parse atom mapping from letter notation data without string building.

        Args:
            reactant_items: List of (compound_id, letter_cfg) tuples
            product_items: List of (product_id, letter_cfg) tuples
        """
        # Build a global list of all atoms from all reactants
        # (cpd, atom_idx, instance)
        all_reactant_atoms: List[Tuple[str, int, int]] = []
        # Sequence of all letters from reactants
        letter_sequence: List[str] = []

        # Count instances for each compound on reactant side
        reactant_instance_counts: Dict[str, int] = {}
        for cpd, letters in reactant_items:
            reactant_instance_counts[cpd] = (
                reactant_instance_counts.get(cpd, 0) + 1
            )
            instance = reactant_instance_counts[cpd]
            for i, ch in enumerate(letters, start=1):
                all_reactant_atoms.append((cpd, i, instance))
                letter_sequence.append(ch)

        # Build mapping for products
        atom_map = AtomMap()
        letter_usage_count: Dict[str, int] = {}

        # Count instances for each compound on product side
        product_instance_counts: Dict[str, int] = {}
        for cpd, letters in product_items:
            product_instance_counts[cpd] = (
                product_instance_counts.get(cpd, 0) + 1
            )
            instance = product_instance_counts[cpd]
            for i, ch in enumerate(letters, start=1):
                # Count how many times we've used this letter
                usage_count = letter_usage_count.get(ch, 0)

                # Find the nth occurrence of this letter in the sequence
                occurrences = [
                    idx
                    for idx, letter in enumerate(letter_sequence)
                    if letter == ch
                ]

                if usage_count >= len(occurrences):
                    raise ValueError(
                        f"Letter '{ch}' used more times in products "
                        f"than available in reactants. Available: "
                        f"{len(occurrences)}, Used: {usage_count + 1}"
                    )

                atom_idx = occurrences[usage_count]
                src_cpd, src_atom, src_inst = all_reactant_atoms[atom_idx]

                atom_map = atom_map.add(
                    cpd, i, src_cpd, src_atom, None, instance, src_inst
                )

                letter_usage_count[ch] = usage_count + 1

        return cls(maps=[atom_map])

    def to_letter_notation(self) -> str:
        """Convert atom mapping to letter notation string.

        Returns a string representation of the atom mapping in letter
        notation format. If multiple elements are present, they are shown
        on separate lines. If multiple mapping variants exist, they are
        shown as alternatives.

        Example: "A(abc) + B(de) -> C(abcde)"
        With elements: "C: A(ab) -> B(ab)\nN: A(c) -> B(c)"

        Returns:
            Letter notation string representation of the mapping(s)
        """
        from collections import defaultdict

        if not self.maps:
            return ""

        # Handle multiple mapping variants
        if len(self.maps) > 1:
            # For multiple variants, show each on a separate line
            results = []
            for i, atom_map in enumerate(self.maps):
                weight_str = f" [{self.weights[i]:.3f}]" if self.weights else ""
                single_mapping = AtomMapping(maps=[atom_map])
                result = single_mapping.to_letter_notation()
                if result:
                    results.append(f"Variant {i+1}{weight_str}: {result}")
            return "\n".join(results)

        # Single mapping
        atom_map = self.maps[0]
        if not atom_map.mapping:
            return ""

        # Group mappings by element type
        element_groups = defaultdict(list)
        for prod_addr, react_addr in atom_map.mapping.items():
            element = prod_addr.element or "C"
            element_groups[element].append((prod_addr, react_addr))

        # Process each element separately
        element_results = []

        for element in sorted(element_groups.keys()):
            mappings = element_groups[element]

            # Collect compounds and their atoms for this element
            reactant_atoms = defaultdict(lambda: defaultdict(list))
            product_atoms = defaultdict(lambda: defaultdict(list))

            # Assign letters sequentially based on reactant order
            letters = (
                "abcdefghijklmnopqrstuvwxyz"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
            )
            letter_idx = 0
            atom_to_letter = {}

            # Sort mappings by reactant position for consistent
            # letter assignment
            sorted_mappings = sorted(
                mappings, key=lambda x: (x[1].mol, x[1].instance, x[1].index)
            )

            for prod_addr, react_addr in sorted_mappings:
                key = (react_addr.mol, react_addr.instance, react_addr.index)
                if key not in atom_to_letter:
                    if letter_idx >= len(letters):
                        # Generate extended letters like aa, ab, ac...
                        div, mod = divmod(letter_idx - len(letters), 26)
                        atom_to_letter[key] = letters[div % 26] + letters[mod]
                    else:
                        atom_to_letter[key] = letters[letter_idx]
                    letter_idx += 1

                letter = atom_to_letter[key]
                reactant_atoms[react_addr.mol][react_addr.instance].append(
                    (react_addr.index, letter)
                )
                product_atoms[prod_addr.mol][prod_addr.instance].append(
                    (prod_addr.index, letter)
                )

            # Build reaction string for this element
            def build_side(atoms_dict):
                """Build one side of the reaction string."""
                parts = []
                for cpd in sorted(atoms_dict.keys()):
                    instances = atoms_dict[cpd]
                    for instance in sorted(instances.keys()):
                        atoms = sorted(instances[instance])
                        letter_str = "".join(ltr for _, ltr in atoms)
                        parts.append(f"{cpd}({letter_str})")
                return parts

            reactant_parts = build_side(reactant_atoms)
            product_parts = build_side(product_atoms)

            if reactant_parts and product_parts:
                reaction = (
                    f"{' + '.join(reactant_parts)} -> "
                    f"{' + '.join(product_parts)}"
                )
                # Only add element prefix if there are multiple elements
                if len(element_groups) > 1:
                    reaction = f"{element}: {reaction}"
                element_results.append(reaction)

        # Join results
        if len(element_results) == 0:
            return ""
        elif len(element_results) == 1:
            return element_results[0]
        else:
            return "\n".join(element_results)

    @classmethod
    def from_fluxml_cfg(
        cls,
        reactant_cfgs: Dict[str, str],
        product_cfgs: List[Tuple[str, str]],
        reactant_order: List[str],
    ) -> "AtomMapping":
        """Parse FluxML cfg strings.

        Args:
            reactant_cfgs: Dict of reactant_id -> cfg string
            product_cfgs: List of (product_id, cfg string) tuples
            reactant_order: Ordered list of reactant IDs
                (with repetitions for instances)

        In product cfg "C#1@2 C#3@1", each token means:
        - C#k@r: take atom k from reactant at position r
        """
        atom_map = AtomMap()

        # Build instance mapping: position -> (compound, instance)
        # This handles repeated compounds in reactant_order
        position_to_instance = {}
        compound_counts = {}
        for idx, cpd in enumerate(reactant_order):
            compound_counts[cpd] = compound_counts.get(cpd, 0) + 1
            position_to_instance[idx] = (cpd, compound_counts[cpd])

        # Track product instances
        product_instances = {}
        product_atom_indices = {}

        # Process product cfgs
        for pcid, cfg in product_cfgs:
            # Track instance number for this product
            product_instances[pcid] = product_instances.get(pcid, 0) + 1
            prod_instance = product_instances[pcid]

            # Track atom index for this product instance
            key = (pcid, prod_instance)
            if key not in product_atom_indices:
                product_atom_indices[key] = 0

            for tok in cfg.split():
                m = re.fullmatch(r"([A-Z])#(\d+)@(\d+)", tok)
                if not m:
                    raise ValueError(f"Bad product cfg token: {tok}")
                atom_type = m.group(1)
                src_atom = int(m.group(2))
                src_reactant_idx = int(m.group(3)) - 1  # to 0-based

                if src_reactant_idx not in position_to_instance:
                    raise IndexError(
                        f"Reactant index @{src_reactant_idx+1} out of range"
                    )

                src_cpd, src_instance = position_to_instance[src_reactant_idx]
                product_atom_indices[key] += 1
                prod_atom_idx = product_atom_indices[key]

                atom_map = atom_map.add(
                    pcid,
                    prod_atom_idx,
                    src_cpd,
                    src_atom,
                    atom_type,
                    prod_instance,
                    src_instance,
                )

        return cls(maps=[atom_map])

    def to_fluxml_string(
        self,
        reactants: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
    ) -> Optional[str]:
        """Convert to FluxML-style atom mapping string.

        Args:
            reactants: Optional list of reactant IDs. If None,
                derives from maps.
            products: Optional list of product IDs. If None,
                derives from maps.

        Returns:
            FluxML-style string representation of the mapping
        """
        if not self.maps:
            return None

        # Use first mapping (for variants, caller should iterate)
        atom_map = self.maps[0]

        # Derive reactants and products from mapping if not provided
        if reactants is None:
            reactants = sorted(
                set(addr.mol for addr in atom_map.mapping.values())
            )
        if products is None:
            products = sorted(set(addr.mol for addr in atom_map.mapping.keys()))

        # Build reactant and product strings
        reactant_parts = []
        product_parts = []

        # Group mappings by compound and instance
        reactant_instances: Dict[str, set] = {}
        for addr in atom_map.mapping.values():
            if addr.mol not in reactant_instances:
                reactant_instances[addr.mol] = set()
            reactant_instances[addr.mol].add(addr.instance)

        for cpd in reactants:
            if cpd in reactant_instances:
                for instance in sorted(reactant_instances[cpd]):
                    atoms = []
                    for addr in atom_map.mapping.values():
                        if addr.mol == cpd and addr.instance == instance:
                            atoms.append((addr.index, addr.element or "C"))
                    if atoms:
                        atom_str = " ".join(
                            f"{at}#{idx}" for idx, at in sorted(atoms)
                        )
                        num_inst = len(reactant_instances[cpd])
                        prefix = f"{num_inst}" if num_inst > 1 else ""
                        reactant_parts.append(f"{prefix}{cpd}({atom_str})")

        product_instances: Dict[str, set] = {}
        for addr in atom_map.mapping.keys():
            if addr.mol not in product_instances:
                product_instances[addr.mol] = set()
            product_instances[addr.mol].add(addr.instance)

        for cpd in products:
            if cpd in product_instances:
                for instance in sorted(product_instances[cpd]):
                    atoms = []
                    for prod_addr in atom_map.mapping.keys():
                        if (
                            prod_addr.mol == cpd
                            and prod_addr.instance == instance
                        ):
                            src_addr = atom_map.mapping[prod_addr]
                            atoms.append(
                                (
                                    prod_addr.index,
                                    prod_addr.element or "C",
                                    src_addr.mol,
                                    src_addr.index,
                                )
                            )
                    if atoms:
                        atom_str = " ".join(
                            f"{at}#{pa}" for pa, at, _, _ in sorted(atoms)
                        )
                        num_inst = len(product_instances[cpd])
                        prefix = f"{num_inst}" if num_inst > 1 else ""
                        product_parts.append(f"{prefix}{cpd}({atom_str})")

        if not reactant_parts or not product_parts:
            return None

        reactant_side = " + ".join(reactant_parts)
        product_side = " + ".join(product_parts)
        return f'"{reactant_side} => {product_side}"'

    def merge_symmetric_mappings(self, other: "AtomMapping") -> "AtomMapping":
        """Merge with another mapping to handle symmetric reactions."""
        new_maps = list(self.maps) + list(other.maps)
        # If weights exist, combine them proportionally
        if self.weights and other.weights:
            total = len(self.maps) + len(other.maps)
            self_weight = len(self.maps) / total
            other_weight = len(other.maps) / total
            new_weights = [w * self_weight for w in self.weights] + [
                w * other_weight for w in other.weights
            ]
        else:
            new_weights = None

        return AtomMapping(maps=new_maps, weights=new_weights)

    def transform_isotopomers(
        self,
        reactant_distributions: Dict[Tuple[str, int], np.ndarray],
        reactant_order: List[str],
        product_order: List[str],
        which_map: int = 0,
    ) -> Dict[Tuple[str, int], np.ndarray]:
        """Transform isotopomer distributions through the reaction.

        Args:
            reactant_distributions: Dict of (compound_id, instance) ->
                isotopomer distribution
            reactant_order: Order of reactants (with repetitions for multiple
                instances)
            product_order: Order of products (with repetitions for multiple
                instances)
            which_map: Which mapping variant to use

        Returns:
            Dict of (product_id, instance) -> isotopomer distribution
        """
        if which_map >= len(self.maps):
            raise ValueError(f"Map index {which_map} out of range")

        atom_map = self.maps[which_map]

        # Build concatenated reactant atom order with instances
        src_atoms: List[Tuple[str, int, str, int]] = []
        # Count instances per compound
        reactant_instances: Dict[str, int] = {}
        for cpd in reactant_order:
            reactant_instances[cpd] = reactant_instances.get(cpd, 0) + 1
            instance = reactant_instances[cpd]
            key = (cpd, instance)
            if key not in reactant_distributions:
                raise ValueError(f"Missing distribution for {key}")
            n_atoms = len(reactant_distributions[key]).bit_length() - 1
            for i in range(1, n_atoms + 1):
                src_atoms.append((cpd, i, "C", instance))

        # Build product atom order with instances
        prod_atoms: List[Tuple[str, int, str, int]] = []
        product_instances: Dict[str, int] = {}
        for cpd in product_order:
            product_instances[cpd] = product_instances.get(cpd, 0) + 1
            instance = product_instances[cpd]
            # Find max atom index for this product instance
            max_idx = 0
            for p_cpd, p_atom, _, p_inst in atom_map.mapping.keys():
                if p_cpd == cpd and p_inst == instance:
                    max_idx = max(max_idx, p_atom)
            for i in range(1, max_idx + 1):
                prod_atoms.append((cpd, i, "C", instance))

        # Create mapping from product atom positions to source positions
        src_pos_lookup = {atom: i for i, atom in enumerate(src_atoms)}
        src_positions = []
        for prod_atom in prod_atoms:
            if prod_atom not in atom_map.mapping:
                raise ValueError(f"Product atom {prod_atom} not in mapping")
            src_atom = atom_map.mapping[prod_atom]
            src_positions.append(src_pos_lookup[src_atom])

        # Build joint reactant distribution
        joint_dist = np.array([1.0])
        instance_counts: Dict[str, int] = {}
        for cpd in reactant_order:
            instance_counts[cpd] = instance_counts.get(cpd, 0) + 1
            instance = instance_counts[cpd]
            key = (cpd, instance)
            joint_dist = np.kron(joint_dist, reactant_distributions[key])

        # Transform to product distribution
        n_src = len(src_atoms)
        n_prod = len(prod_atoms)

        if joint_dist.size != (1 << n_src):
            raise ValueError("Reactant distribution size mismatch")

        # Compute product distribution
        prod_dist = np.zeros(1 << n_prod)
        src_indices = np.arange(1 << n_src, dtype=np.uint32)

        # Map bits from source to product positions
        prod_indices = np.zeros_like(src_indices)
        for k, src_pos in enumerate(src_positions):
            bit = (src_indices >> src_pos) & 1
            prod_indices |= bit << k

        np.add.at(prod_dist, prod_indices, joint_dist)

        # Split into per-product distributions
        result = {}
        offset = 0
        product_instance_counts: Dict[str, int] = {}
        for cpd in product_order:
            product_instance_counts[cpd] = (
                product_instance_counts.get(cpd, 0) + 1
            )
            instance = product_instance_counts[cpd]
            n_atoms = sum(
                1
                for (c, _, _, inst) in prod_atoms
                if c == cpd and inst == instance
            )
            if n_atoms > 0:
                size = 1 << n_atoms
                result[(cpd, instance)] = prod_dist[offset : offset + size]
                offset += size

        return result

    def transform_with_symmetry(
        self,
        reactant_distributions: Dict[Tuple[str, int], np.ndarray],
        reactant_order: List[str],
        product_order: List[str],
    ) -> Dict[Tuple[str, int], np.ndarray]:
        """Average isotopomer transformation over all mapping variants."""
        if not self.maps:
            raise ValueError("No atom maps available")

        weights = (
            self.weights
            if self.weights
            else [1.0 / len(self.maps)] * len(self.maps)
        )

        result: Dict[Tuple[str, int], np.ndarray] = {}
        for i, weight in enumerate(weights):
            prod_dists = self.transform_isotopomers(
                reactant_distributions, reactant_order, product_order, i
            )
            for key, dist in prod_dists.items():
                if key not in result:
                    result[key] = weight * dist
                else:
                    result[key] += weight * dist

        return result
