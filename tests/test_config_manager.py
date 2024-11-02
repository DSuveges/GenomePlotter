from __future__ import annotations

import json
import re
import unittest

from functions.ConfigManager import Config


class TestConfigManager(unittest.TestCase):
    CONFIG_JSON = "config.json"

    with open(CONFIG_JSON, "w") as f:
        config_obj = Config(**json.load(f))
    hex_color_match = re.compile(r"^#[0-9A-F]{6}$", re.IGNORECASE)

    # def save_config(self, filename=None):
    # def set_data_folder(self, data_folder):
    # def set_width(self, width):
    # def set_pixel_size(self, pixel_size):
    # def set_dark_start(self, dark_start):
    # def set_dark_max(self, dark_max):
    # def set_plot_folder(self, plot_folder):
    # -> def get_chromosome_colors(self):
    # -> def get_cytobanc_colors(self):
    # -> def get_gwas_color(self):
    # -> def get_arrow_colors(self):
    # def get_data_folder(self):
    # def get_source(self, resource, key, chromosome=None):
    # def get_cytoband_file(self):
    # def get_chromosome_file(self, chromsome):
    # def get_gwas_file(self):
    # def get_gencode_file(self):
    # def get_pixel(self):
    # def get_chunk_size(self):
    # def get_width(self):
    # def get_dark_max(self):
    # def get_dark_start(self):

    def test_get_chromosome_colors(self):
        chromosome_colors = self.config_obj.get_chromosome_colors()

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

    def test_get_cytobanc_colors(self):
        cytobanc_colors = self.config_obj.get_cytobanc_colors()

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

    def test_get_arrow_colors(self):
        arrow_colors = self.config_obj.get_arrow_colors()

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

    def test_get_gwas_color(self):
        gwas_color = self.config_obj.get_gwas_color()
        self.assertIsInstance(gwas_color, str)
        self.assertTrue(self.hex_color_match.match(gwas_color))

    def test_get_custom_gene_window(self):
        gene_window = self.config_obj.get_custom_gene_window()
        self.assertIsInstance(gene_window, int)


if __name__ == "__main__":
    unittest.main()
