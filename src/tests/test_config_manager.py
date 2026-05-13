"""Tests for the ConfigManager module."""

from __future__ import annotations

import os
import re
import unittest

import yaml

from genome_plotter.functions.ConfigManager import Config


class TestConfigManager(unittest.TestCase):
    """Test cases for Config manager."""

    CONFIG_YAML = os.path.join(
        os.path.dirname(__file__), "..", "genome_plotter", "assets", "config.yaml"
    )

    with open(CONFIG_YAML) as f:
        config_obj = Config.model_validate(yaml.safe_load(f))
    hex_color_match = re.compile(r"^#[0-9A-F]{6}$", re.IGNORECASE)

    def test_get_chromosome_colors(self) -> None:
        """Test chromosome colors retrieval."""
        chromosome_colors = self.config_obj.color_schema.chromosome_colors

        # Test return types:
        self.assertIsInstance(chromosome_colors, dict)

        # Fields are checked:
        chromosome_features = [
            "centromere",
            "heterochromatin",
            "intergenic",
            "exon",
            "gene",
            "dummy",
        ]
        for feature in chromosome_features:
            self.assertIn(feature, chromosome_colors)

        for feature in chromosome_features:
            color = chromosome_colors[feature]
            self.assertIsInstance(color, str)
            self.assertTrue(self.hex_color_match.match(color))

    def test_get_cytobanc_colors(self) -> None:
        """Test cytoband colors retrieval."""
        cytobanc_colors = self.config_obj.color_schema.cytoband_colors

        # Test return types:
        self.assertIsInstance(cytobanc_colors, dict)

        # Fields are checked:
        chromosome_features = [
            "border",
            "stalk",
            "gvar",
            "gneg",
            "gpos25",
            "gpos50",
            "gpos75",
            "gpos100",
            "acen",
        ]
        for feature in chromosome_features:
            self.assertIn(feature, cytobanc_colors)

        for feature in chromosome_features:
            color = cytobanc_colors[feature]
            self.assertIsInstance(color, str)
            self.assertTrue(self.hex_color_match.match(color))

    def test_get_arrow_colors(self) -> None:
        """Test arrow colors retrieval."""
        arrow_colors = self.config_obj.color_schema.arrow_colors

        # Test return types:
        self.assertIsInstance(arrow_colors, dict)

        # Fields are checked:
        chromosome_features = ["line_color", "utr_color", "cds_color"]
        for feature in chromosome_features:
            self.assertIn(feature, arrow_colors)

        for feature in chromosome_features:
            color = arrow_colors[feature]
            self.assertIsInstance(color, str)
            self.assertTrue(self.hex_color_match.match(color))

    def test_get_gwas_color(self) -> None:
        """Test GWAS color retrieval."""
        gwas_color = self.config_obj.color_schema.gwas_point
        self.assertIsInstance(gwas_color, str)
        self.assertTrue(self.hex_color_match.match(gwas_color))

    def test_get_custom_gene_window(self) -> None:
        """Test custom gene window retrieval."""
        gene_window = self.config_obj.plot_parameters.custom_gene_window
        self.assertIsInstance(gene_window, int)


def _fresh_config() -> Config:
    """Load a fresh Config instance from the bundled YAML for each test."""
    _yaml = os.path.join(
        os.path.dirname(__file__), "..", "genome_plotter", "assets", "config.yaml"
    )
    with open(_yaml) as f:
        return Config.model_validate(yaml.safe_load(f))


class TestUpdateBasicParameters(unittest.TestCase):
    """Tests for Config.update_basic_parameters."""

    def test_known_field_is_updated(self) -> None:
        """A known BasicParameters field is updated to the supplied value."""
        config = _fresh_config()
        config.update_basic_parameters(chunk_size=900)
        self.assertEqual(config.basic_parameters.chunk_size, 900)

    def test_unknown_key_is_silently_ignored(self) -> None:
        """An unrecognised key is silently ignored and nothing else changes."""
        config = _fresh_config()
        original = config.basic_parameters.chunk_size
        config.update_basic_parameters(nonexistent_field=42)
        self.assertEqual(config.basic_parameters.chunk_size, original)

    def test_multiple_fields_updated_at_once(self) -> None:
        """Multiple known fields can be updated in a single call."""
        config = _fresh_config()
        config.update_basic_parameters(chunk_size=200, missing_tolerance=0.1)
        self.assertEqual(config.basic_parameters.chunk_size, 200)
        tolerance = config.basic_parameters.missing_tolerance
        self.assertIsNotNone(tolerance)
        assert tolerance is not None
        self.assertAlmostEqual(tolerance, 0.1)

    def test_optional_field_accepts_string_value(self) -> None:
        """An optional string field accepts a non-None string value."""
        config = _fresh_config()
        config.update_basic_parameters(data_folder="/tmp/data")
        self.assertEqual(config.basic_parameters.data_folder, "/tmp/data")

    def test_optional_field_can_be_reset_to_none(self) -> None:
        """An optional field can be reset back to None."""
        config = _fresh_config()
        config.update_basic_parameters(data_folder="/tmp/data")
        config.update_basic_parameters(data_folder=None)
        self.assertIsNone(config.basic_parameters.data_folder)

    def test_mix_of_known_and_unknown_keys(self) -> None:
        """Known keys are updated while unknown keys are ignored."""
        config = _fresh_config()
        config.update_basic_parameters(chunk_size=333, bogus=True)
        self.assertEqual(config.basic_parameters.chunk_size, 333)

    def test_zero_tolerance_is_accepted(self) -> None:
        """0.0 is a valid tolerance value and must not be treated as falsy."""
        config = _fresh_config()
        config.update_basic_parameters(missing_tolerance=0.0)
        self.assertEqual(config.basic_parameters.missing_tolerance, 0.0)


if __name__ == "__main__":
    unittest.main()
