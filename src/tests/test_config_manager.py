"""Tests for the ConfigManager module."""

from __future__ import annotations

import json
import os
import re
import unittest

from genome_plotter.functions.ConfigManager import Config


class TestConfigManager(unittest.TestCase):
    """Test cases for Config manager."""

    CONFIG_JSON = os.path.join(
        os.path.dirname(__file__), "..", "genome_plotter", "assets", "config.json"
    )

    with open(CONFIG_JSON) as f:
        config_obj = Config(**json.load(f))
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


if __name__ == "__main__":
    unittest.main()
