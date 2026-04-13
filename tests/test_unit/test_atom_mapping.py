"""
Test atom mapping functionality including variant handling.
"""

import pytest
import numpy as np
from fluxomics_data_model.model.atom_mapping import (
    AtomMapping,
    AtomMap,
    AtomAddress,
)
from fluxomics_data_model.io.fluxml_parser import FluxMLParser


class TestAtomMappingBasic:
    """Test basic atom mapping functionality."""

    def test_parse_letter_notation_simple(self):
        """Test parsing simple letter notation."""
        reactant_items = [("A", "abc")]
        product_items = [("B", "abc")]

        atom_map = AtomMapping.parse_letter_notation(
            reactant_items, product_items
        )

        # Check that mapping exists
        assert len(atom_map.mapping) == 3

        # Verify mappings
        assert atom_map.mapping[AtomAddress("B", 1, None, 1)] == AtomAddress(
            "A", 1, None, 1
        )
        assert atom_map.mapping[AtomAddress("B", 2, None, 1)] == AtomAddress(
            "A", 2, None, 1
        )
        assert atom_map.mapping[AtomAddress("B", 3, None, 1)] == AtomAddress(
            "A", 3, None, 1
        )

    def test_parse_letter_notation_complex(self):
        """Test parsing complex letter notation with reordering."""
        reactant_items = [("A", "abc"), ("B", "de")]
        product_items = [("C", "edc"), ("D", "ba")]

        atom_map = AtomMapping.parse_letter_notation(
            reactant_items, product_items
        )

        # Check that mapping exists
        assert len(atom_map.mapping) == 5

        # Verify mappings - e, d, c from second and first reactant
        assert atom_map.mapping[AtomAddress("C", 1, None, 1)] == AtomAddress(
            "B", 2, None, 1
        )  # e
        assert atom_map.mapping[AtomAddress("C", 2, None, 1)] == AtomAddress(
            "B", 1, None, 1
        )  # d
        assert atom_map.mapping[AtomAddress("C", 3, None, 1)] == AtomAddress(
            "A", 3, None, 1
        )  # c

    def test_to_letter_notation(self):
        """Test converting atom mapping to letter notation."""
        reactant_items = [("A", "abc")]
        product_items = [("B", "abc")]

        atom_map = AtomMapping.parse_letter_notation(
            reactant_items, product_items
        )
        atom_mapping = AtomMapping(
            reaction_id="test",
            reactants=["A"],
            products=["B"],
            maps={"test": atom_map},
        )

        notation = atom_mapping.to_letter_notation()
        assert notation == "A(abc) -> B(abc)"

    def test_parse_fluxml_cfg_simple(self):
        """Test parsing FluxML cfg strings."""
        reactant_cfgs = {"A": "C#1 C#2 C#3"}
        product_cfgs = [("B", "C#1@1 C#2@1 C#3@1")]
        reactant_order = ["A"]

        atom_map = AtomMapping.parse_fluxml_cfg(
            reactant_cfgs=reactant_cfgs,
            product_cfgs=product_cfgs,
            reactant_order=reactant_order,
        )

        # Check that mapping exists
        assert len(atom_map.mapping) == 3

        # Verify mappings
        assert atom_map.mapping[AtomAddress("B", 1, "C", 1)] == AtomAddress(
            "A", 1, "C", 1
        )
        assert atom_map.mapping[AtomAddress("B", 2, "C", 1)] == AtomAddress(
            "A", 2, "C", 1
        )
        assert atom_map.mapping[AtomAddress("B", 3, "C", 1)] == AtomAddress(
            "A", 3, "C", 1
        )

    def test_parse_fluxml_cfg_reordering(self):
        """Test parsing FluxML cfg with atom reordering."""
        reactant_cfgs = {"A": "C#1 C#2 C#3"}
        product_cfgs = [("B", "C#3@1 C#2@1 C#1@1")]
        reactant_order = ["A"]

        atom_map = AtomMapping.parse_fluxml_cfg(
            reactant_cfgs=reactant_cfgs,
            product_cfgs=product_cfgs,
            reactant_order=reactant_order,
        )

        # Check that mapping exists
        assert len(atom_map.mapping) == 3

        # Verify reversed mappings
        assert atom_map.mapping[AtomAddress("B", 1, "C", 1)] == AtomAddress(
            "A", 3, "C", 1
        )
        assert atom_map.mapping[AtomAddress("B", 2, "C", 1)] == AtomAddress(
            "A", 2, "C", 1
        )
        assert atom_map.mapping[AtomAddress("B", 3, "C", 1)] == AtomAddress(
            "A", 1, "C", 1
        )


class TestAtomMappingVariants:
    """Test atom mapping with variants."""

    def test_atom_mapping_with_multiple_variants(self):
        """Test AtomMapping with multiple variants."""
        # Create two simple atom maps
        map1 = AtomMap()
        map1 = map1.add("B", 1, "A", 1, "C", 1, 1)
        map1 = map1.add("B", 2, "A", 2, "C", 1, 1)

        map2 = AtomMap()
        map2 = map2.add("B", 1, "A", 2, "C", 1, 1)
        map2 = map2.add("B", 2, "A", 1, "C", 1, 1)

        # Create AtomMapping with variants
        atom_mapping = AtomMapping(
            reaction_id="test",
            reactants=["A"],
            products=["B"],
            maps={"test___1": map1, "test___2": map2},
            weights={"test___1": 0.5, "test___2": 0.5},
        )

        # Check that we have 2 variants
        assert len(atom_mapping.maps) == 2
        assert "test___1" in atom_mapping.maps
        assert "test___2" in atom_mapping.maps

        # Check weights
        assert atom_mapping.weights["test___1"] == 0.5
        assert atom_mapping.weights["test___2"] == 0.5

    def test_to_letter_notation_multiple_variants(self):
        """Test converting multiple variants to letter notation."""
        # Create two simple atom maps
        map1 = AtomMap()
        map1 = map1.add("B", 1, "A", 1, "C", 1, 1)
        map1 = map1.add("B", 2, "A", 2, "C", 1, 1)

        map2 = AtomMap()
        map2 = map2.add("B", 1, "A", 2, "C", 1, 1)
        map2 = map2.add("B", 2, "A", 1, "C", 1, 1)

        atom_mapping = AtomMapping(
            reaction_id="test",
            reactants=["A"],
            products=["B"],
            maps={"test___1": map1, "test___2": map2},
            weights={"test___1": 0.5, "test___2": 0.5},
        )

        notation = atom_mapping.to_letter_notation()

        # Should show both variants
        assert "test___1" in notation
        assert "test___2" in notation
        assert "[0.500]" in notation  # weight display


class TestIsotopomerTransformation:
    """Test isotopomer transformation functionality."""

    def test_transform_isotopomers_simple(self):
        """Test simple isotopomer transformation."""
        # Create a simple mapping: A(abc) -> B(abc)
        atom_map = AtomMap()
        atom_map = atom_map.add("B", 1, "A", 1, "C", 1, 1)
        atom_map = atom_map.add("B", 2, "A", 2, "C", 1, 1)
        atom_map = atom_map.add("B", 3, "A", 3, "C", 1, 1)

        atom_mapping = AtomMapping(
            reaction_id="test",
            reactants=["A"],
            products=["B"],
            maps={"test": atom_map},
        )

        # Create input distribution: all unlabeled
        reactant_dist = {("A", 1): np.array([1.0, 0, 0, 0, 0, 0, 0, 0])}

        # Transform
        product_dist = atom_mapping.transform_isotopomers(
            reactant_dist, atom_map_id="test"
        )

        # Output should be identical
        assert ("B", 1) in product_dist
        assert product_dist[("B", 1)][0] == 1.0
        assert sum(product_dist[("B", 1)][1:]) == 0.0

    def test_transform_with_symmetry(self):
        """Test transformation with multiple symmetric variants."""
        # Create two variants: normal and reversed
        map1 = AtomMap()
        map1 = map1.add("B", 1, "A", 1, "C", 1, 1)
        map1 = map1.add("B", 2, "A", 2, "C", 1, 1)

        map2 = AtomMap()
        map2 = map2.add("B", 1, "A", 2, "C", 1, 1)
        map2 = map2.add("B", 2, "A", 1, "C", 1, 1)

        atom_mapping = AtomMapping(
            reaction_id="test",
            reactants=["A"],
            products=["B"],
            maps={"test___1": map1, "test___2": map2},
            weights={"test___1": 0.5, "test___2": 0.5},
        )

        # Input: 50% labeled at position 1, 50% at position 2
        # Using 2-carbon compound (4 isotopomers: 00, 01, 10, 11)
        reactant_dist = {("A", 1): np.array([0.0, 0.5, 0.5, 0.0])}

        # Transform with symmetry
        product_dist = atom_mapping.transform_with_symmetry(reactant_dist)

        # With symmetry, both variants should average out
        assert ("B", 1) in product_dist
        # Due to symmetry, distribution should be symmetric
        assert np.isclose(product_dist[("B", 1)][1], product_dist[("B", 1)][2])


class TestFluxMLParsingVariants:
    """Test parsing FluxML files with variants."""

    def test_parse_bsDAP_reaction(self):
        """Test parsing bsDAP reaction with 4 variants (2×2)."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get bsDAP reaction
        reaction = data_model.model.reactions.get_by_id("bsDAP")
        assert reaction is not None

        # Should have 4 variant IDs
        assert reaction.atom_mapping_ids is not None
        assert len(reaction.atom_mapping_ids) == 4
        assert "bsDAP___1" in reaction.atom_mapping_ids
        assert "bsDAP___2" in reaction.atom_mapping_ids
        assert "bsDAP___3" in reaction.atom_mapping_ids
        assert "bsDAP___4" in reaction.atom_mapping_ids

        # Get atom mapping
        atom_mapping = data_model.model.atom_mappings.get("bsDAP")
        assert atom_mapping is not None

        # Should have 4 maps
        assert len(atom_mapping.maps) == 4
        assert "bsDAP___1" in atom_mapping.maps
        assert "bsDAP___2" in atom_mapping.maps
        assert "bsDAP___3" in atom_mapping.maps
        assert "bsDAP___4" in atom_mapping.maps

        # All weights should be equal (0.25 each)
        assert atom_mapping.weights is not None
        for weight in atom_mapping.weights.values():
            assert np.isclose(weight, 0.25)

    def test_parse_bsMET_reaction(self):
        """Test parsing bsMET reaction with 2 variants (1×2)."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get bsMET reaction
        reaction = data_model.model.reactions.get_by_id("bsMET")
        assert reaction is not None

        # Should have 2 variant IDs
        assert reaction.atom_mapping_ids is not None
        assert len(reaction.atom_mapping_ids) == 2
        assert "bsMET___1" in reaction.atom_mapping_ids
        assert "bsMET___2" in reaction.atom_mapping_ids

        # Get atom mapping
        atom_mapping = data_model.model.atom_mappings.get("bsMET")
        assert atom_mapping is not None

        # Should have 2 maps
        assert len(atom_mapping.maps) == 2
        assert "bsMET___1" in atom_mapping.maps
        assert "bsMET___2" in atom_mapping.maps

        # All weights should be equal (0.5 each)
        assert atom_mapping.weights is not None
        for weight in atom_mapping.weights.values():
            assert np.isclose(weight, 0.5)

    def test_to_fluxml_string_variant(self):
        """Test that variant atom mappings have different FluxML string outputs."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get bsMET atom mapping
        atom_mapping = data_model.model.atom_mappings.get("bsMET")
        assert atom_mapping is not None

        # Get both variant maps
        map1 = atom_mapping.maps["bsMET___1"]
        map2 = atom_mapping.maps["bsMET___2"]

        # Check that SUC mappings are different between variants
        # Variant 1: SUC atom 1 should come from SUCCOA atom 2
        from fluxomics_data_model.model.atom_mapping import AtomAddress

        assert map1.mapping[AtomAddress("SUC", 1, "C", 1)] == AtomAddress(
            "SUCCOA", 2, "C", 1
        )

        # Variant 2: SUC atom 1 should come from SUCCOA atom 1
        assert map2.mapping[AtomAddress("SUC", 1, "C", 1)] == AtomAddress(
            "SUCCOA", 1, "C", 1
        )

        # Verify they have different mappings for SUC
        suc_map1 = {k: v for k, v in map1.mapping.items() if k.mol == "SUC"}
        suc_map2 = {k: v for k, v in map2.mapping.items() if k.mol == "SUC"}
        assert suc_map1 != suc_map2

        # Test that FluxML strings are now different
        fluxml_str1 = atom_mapping.to_fluxml_string(atom_map_id="bsMET___1")
        fluxml_str2 = atom_mapping.to_fluxml_string(atom_map_id="bsMET___2")
        assert fluxml_str1 != fluxml_str2

        # Verify that SUC atoms have different sources in the FluxML strings
        # Variant 1 should have SUC with swapped atoms
        assert "SUC(C#2@3 C#1@3 C#4@3 C#3@3)" in fluxml_str1
        # Variant 2 should have SUC in order
        assert "SUC(C#1@3 C#2@3 C#3@3 C#4@3)" in fluxml_str2

    def test_atom_mapping_summary(self):
        """Test atom mapping summary functionality."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Test bsMET with 2 variants
        atom_mapping = data_model.model.atom_mappings.get("bsMET")
        summary = atom_mapping.summary()
        assert "bsMET" in summary
        assert "2 variants" in summary
        assert "bsMET___1" in summary
        assert "bsMET___2" in summary
        assert "0.500" in summary  # weight

        # Test str representation
        str_repr = str(atom_mapping)
        assert str_repr == summary

        # Test bsDAP with 4 variants
        atom_mapping_dap = data_model.model.atom_mappings.get("bsDAP")
        summary_dap = atom_mapping_dap.summary()
        assert "bsDAP" in summary_dap
        assert "4 variants" in summary_dap
        assert "0.250" in summary_dap  # weight


class TestBackwardCompatibility:
    """Test backward compatibility with non-variant reactions."""

    def test_parse_simple_reaction(self):
        """Test parsing simple reaction without variants."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get a simple reaction (e.g., bsGLY)
        reaction = data_model.model.reactions.get_by_id("bsGLY")
        assert reaction is not None

        # Should not have variant IDs (or None)
        assert reaction.atom_mapping_ids is None

        # Get atom mapping
        atom_mapping = data_model.model.atom_mappings.get("bsGLY")
        assert atom_mapping is not None

        # Should have 1 map with reaction ID as key
        assert len(atom_mapping.maps) == 1
        assert "bsGLY" in atom_mapping.maps

    def test_to_letter_notation_simple(self):
        """Test converting simple reaction to letter notation."""
        parser = FluxMLParser()
        data_model = parser.parse_file(
            "/home/te/Projects/data_model/fluxomics_data_model/data/FluxML_tests/fluxml-model/models/CN_network_model_BCG.fml"
        )

        # Get a simple reaction
        atom_mapping = data_model.model.atom_mappings.get("bsGLY")
        assert atom_mapping is not None

        # Convert to letter notation
        notation = atom_mapping.to_letter_notation()
        assert notation is not None
        assert "->" in notation
