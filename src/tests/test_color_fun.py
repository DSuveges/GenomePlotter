"""Tests for the ColorFunctions module."""

from __future__ import annotations

import re
import unittest

import pandas as pd

from genome_plotter.functions.ColorFunctions import (
    ColorPicker,
    color_darkener,
    hex_to_rgb,
    linear_gradient,
    rgb_to_hex,
)


class TestColorFunctions(unittest.TestCase):
    """Test cases for color functions."""

    def test_linear_gradient(self) -> None:
        """Test linear gradient generation."""
        # Testing the default lenght of the gradient:
        gradient = linear_gradient("#000000", "#FFFFFF")
        self.assertEqual(len(gradient), 10)

        # Testing a custom lenght of the gradient:
        gradient = linear_gradient("#000000", "#FFFFFF", 15)
        self.assertEqual(len(gradient), 15)

        # Testing a custom lenght of the gradient:
        gradient = linear_gradient("#000000", "#FFFFFF", 0)
        self.assertEqual(len(gradient), 0)

        # Testing a custom lenght of the gradient:
        with self.assertRaises(ValueError):
            gradient = linear_gradient("#000000", "#FFFFFF", "cica")  # type: ignore[arg-type]

        # Testing if some weird stuff is going on with the input:
        with self.assertRaises(ValueError):
            gradient = linear_gradient("#cica", "#FFFFFF")

        with self.assertRaises(ValueError):
            gradient = linear_gradient("#000000", "#cica")

    def test_hex_to_rgb(self) -> None:
        """Test hexadecimal to RGB conversion."""
        # Test for good output:
        hex_col = "#000000"
        rgb_col = hex_to_rgb(hex_col)

        self.assertIsInstance(rgb_col, list)
        self.assertEqual(len(rgb_col), 3)
        for i in rgb_col:
            self.assertIsInstance(i, int)
        self.assertEqual([0, 0, 0], rgb_col)

        # Test for another good output:
        hex_col = "#ffffff"
        self.assertEqual([255, 255, 255], hex_to_rgb(hex_col))

        # Testing for bad output:
        for bad_input in ["cica", True, 13, "#209345209", "ffffff"]:
            with self.assertRaises(ValueError):
                hex_to_rgb(bad_input)  # type: ignore[arg-type]

    def test_rgb_to_hex(self) -> None:
        """Test RGB to hexadecimal conversion."""
        # Test for good output:
        rgb_col = [0, 0, 0]
        hex_col = rgb_to_hex(rgb_col)

        self.assertIsInstance(hex_col, str)
        self.assertEqual(hex_col[0], "#")

    def test_color_darkener(self) -> None:
        """Test color darkener function."""
        # Basic correct input:
        color = "#DDDDDD"
        x = 234
        width = 200
        threshold = 0.5
        max_diff_value = 0.9

        # Testing for input type:
        with self.assertRaises(TypeError):
            color_darkener(
                "cicaful",
                x,
                width=width,
                threshold=threshold,
                max_diff_value=max_diff_value,
            )
        with self.assertRaises(TypeError):
            color_darkener(
                color,
                x,
                width="pocok",  # type: ignore[arg-type]
                threshold=threshold,
                max_diff_value=max_diff_value,
            )
        with self.assertRaises(TypeError):
            color_darkener(
                color, x, width=width, threshold=3.0, max_diff_value=max_diff_value
            )
        with self.assertRaises(TypeError):
            color_darkener(
                color,
                x,
                width=width,
                threshold="pocok",  # type: ignore[arg-type]
                max_diff_value=max_diff_value,
            )
        with self.assertRaises(TypeError):
            color_darkener(
                color, x, width=width, threshold=threshold, max_diff_value=1232
            )

        # Testing if the threshold is appreciated:
        self.assertEqual(
            color_darkener(color, 50, width, threshold, max_diff_value), color
        )
        self.assertNotEqual(
            color_darkener(color, x, width, threshold, max_diff_value), color
        )

    def test_color_picker(self) -> None:
        """Test ColorPicker class."""
        # Good set of parameters:
        color_map = {
            "centromere": "#9393FF",
            "heterochromatin": "#F9D2C2",
            "intergenic": "#A3E0D1",
            "exon": "#FFD326",
            "gene": "#6CB8CC",
            "dummy": "#B3F29D",
        }
        width = 200
        count = 20
        dark_max = 0.15
        dark_start = 0.75

        with self.assertRaises(TypeError):
            ColorPicker("Pocok", dark_max, dark_start, count, width)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ColorPicker({"color": "color"}, dark_max, dark_start, count, width)
        with self.assertRaises(ValueError):
            color_map_wrong = color_map.copy()
            color_map_wrong["exon"] = "cica"
            ColorPicker(
                colors=color_map_wrong,
                dark_max=dark_max,
                dark_threshold=dark_start,
                count=count,
                width=width,
            )
        with self.assertRaises(ValueError):
            ColorPicker(color_map, "aaaa", dark_start, count, width)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ColorPicker(color_map, 12.0, dark_start, count, width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, -0.2, count, width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, dark_start, "pocok", width)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, dark_start, count, "foo")  # type: ignore[arg-type]

        # Get correct object initialized:
        cp = ColorPicker(color_map, dark_max, dark_start, count, width)

        self.assertIsInstance(cp.color_map, dict)
        self.assertTrue(
            pd.Series(list(color_map.keys())).isin(list(cp.color_map.keys())).all()
        )

        for feature in cp.color_map.keys():
            self.assertIsInstance(cp.color_map[feature], list)
            self.assertEqual(len(cp.color_map[feature]), count)

        # So far, so good. Testing functionality:
        self.assertIsInstance(cp.map_color("dummy", 0.5), str)
        color = cp.map_color("dummy", 0.5)
        self.assertTrue(re.match(r"#[0-9a-f]{6}", color))
        self.assertEqual(cp.map_color("dummy", 0.5), color_map["dummy"].lower())

        # Do we get heterochromatin:
        self.assertEqual(
            cp.map_color("exon", None), color_map["heterochromatin"].lower()
        )

        # No keys for unknown feature:
        self.assertEqual(cp.map_color("cicaful", 0.5), "#000000")

        # Let's test mapper:
        with self.assertRaises(TypeError):
            cp.pick_color("pocok")
        with self.assertRaises(TypeError):
            cp.pick_color(pd.Series({"GC_ratio": 0.3, "GENCODE": "exon"}))

        self.assertEqual(
            cp.pick_color(pd.Series({"GC_ratio": 0.3, "GENCODE": "dummy", "x": 150})),
            color_map["dummy"].lower(),
        )

        self.assertEqual(
            cp.pick_color(pd.Series({"GC_ratio": None, "GENCODE": "exon", "x": 150})),
            color_map["heterochromatin"].lower(),
        )


class TestPickColorsVectorized(unittest.TestCase):
    """Tests for ColorPicker.pick_colors_vectorized."""

    _COLOR_MAP = {
        "centromere": "#9393FF",
        "heterochromatin": "#F9D2C2",
        "intergenic": "#A3E0D1",
        "exon": "#FFD326",
        "gene": "#6CB8CC",
        "dummy": "#B3F29D",
    }

    def _picker(self, **kwargs: object) -> ColorPicker:
        """Return a ColorPicker with darkening enabled unless overridden."""
        defaults: dict[str, object] = {
            "dark_max": 0.15,
            "dark_threshold": 0.75,
            "count": 20,
            "width": 200,
        }
        defaults.update(kwargs)
        return ColorPicker(self._COLOR_MAP, **defaults)  # type: ignore[arg-type]

    def _df(self, rows: list[dict[str, object]]) -> pd.DataFrame:
        """Build a minimal DataFrame for color-picking tests."""
        return pd.DataFrame(rows)

    def _reference(self, picker: ColorPicker, df: pd.DataFrame) -> pd.Series:
        """Compute colors the slow way for comparison."""
        return df.apply(picker.pick_color, axis=1)

    # --- output shape / format ---

    def test_output_is_series_with_matching_index(self) -> None:
        """Output is a pd.Series with the same index as the input DataFrame."""
        cp = self._picker()
        df = self._df([{"GC_ratio": 0.4, "GENCODE": "intergenic", "x": 10}])
        result = cp.pick_colors_vectorized(df)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(list(result.index), list(df.index))

    def test_output_length_matches_input(self) -> None:
        """Output length equals the number of input rows."""
        cp = self._picker()
        df = self._df([
            {"GC_ratio": 0.2, "GENCODE": "exon", "x": 50},
            {"GC_ratio": 0.8, "GENCODE": "gene", "x": 190},
            {"GC_ratio": None, "GENCODE": "heterochromatin", "x": 0},
        ])
        self.assertEqual(len(cp.pick_colors_vectorized(df)), 3)

    def test_all_outputs_are_valid_hex_strings(self) -> None:
        """Every returned value is a six-digit hex color string."""
        cp = self._picker()
        features = ["exon", "gene", "intergenic", "centromere", "heterochromatin", "dummy"]
        rows = [{"GC_ratio": 0.5, "GENCODE": f, "x": 100} for f in features]
        result = cp.pick_colors_vectorized(self._df(rows))
        pattern = re.compile(r"^#[0-9a-f]{6}$")
        for color in result:
            self.assertRegex(color, pattern)

    # --- correctness vs apply(pick_color) ---

    def test_matches_apply_for_all_features(self) -> None:
        """Vectorised output matches apply(pick_color) for every feature type."""
        cp = self._picker()
        features = ["exon", "gene", "intergenic", "centromere", "heterochromatin", "dummy"]
        rows = [{"GC_ratio": 0.5, "GENCODE": f, "x": 50} for f in features]
        df = self._df(rows)
        pd.testing.assert_series_equal(
            cp.pick_colors_vectorized(df).reset_index(drop=True),
            self._reference(cp, df).reset_index(drop=True),
            check_dtype=False,
        )

    def test_matches_apply_across_gc_values(self) -> None:
        """Vectorised output matches apply(pick_color) at low, mid, and high GC ratios."""
        cp = self._picker()
        rows = [
            {"GC_ratio": gc, "GENCODE": "exon", "x": 10}
            for gc in [0.0, 0.25, 0.5, 0.75, 1.0]
        ]
        df = self._df(rows)
        pd.testing.assert_series_equal(
            cp.pick_colors_vectorized(df).reset_index(drop=True),
            self._reference(cp, df).reset_index(drop=True),
            check_dtype=False,
        )

    def test_matches_apply_below_and_above_dark_threshold(self) -> None:
        """Vectorised output matches apply(pick_color) on both sides of the darkening threshold."""
        cp = self._picker(dark_threshold=0.75, width=200)
        rows = [
            {"GC_ratio": 0.5, "GENCODE": "exon", "x": x}
            for x in [0, 50, 100, 149, 150, 151, 175, 199, 200]
        ]
        df = self._df(rows)
        pd.testing.assert_series_equal(
            cp.pick_colors_vectorized(df).reset_index(drop=True),
            self._reference(cp, df).reset_index(drop=True),
            check_dtype=False,
        )

    def test_nan_gc_always_returns_heterochromatin_color(self) -> None:
        """A NaN GC_ratio returns heterochromatin[0] regardless of the GENCODE label."""
        cp = self._picker()
        df = self._df([{"GC_ratio": None, "GENCODE": "exon", "x": 50}])
        result = cp.pick_colors_vectorized(df).iloc[0]
        expected = cp.pick_color(pd.Series({"GC_ratio": None, "GENCODE": "exon", "x": 50}))
        self.assertEqual(result, expected)

    def test_dummy_feature_never_darkened(self) -> None:
        """Dummy chunks return the same color regardless of x position."""
        cp = self._picker()
        rows = [{"GC_ratio": 0.5, "GENCODE": "dummy", "x": x} for x in [0, 100, 200]]
        result = cp.pick_colors_vectorized(self._df(rows))
        self.assertEqual(result.iloc[0], result.iloc[1])
        self.assertEqual(result.iloc[0], result.iloc[2])

    def test_works_without_darkening(self) -> None:
        """Vectorised lookup works when darkening parameters are not set."""
        cp = ColorPicker(self._COLOR_MAP, count=20)
        df = self._df([
            {"GC_ratio": 0.3, "GENCODE": "exon", "x": 0},
            {"GC_ratio": 0.7, "GENCODE": "gene", "x": 0},
        ])
        result = cp.pick_colors_vectorized(df)
        self.assertEqual(len(result), 2)
        pattern = re.compile(r"^#[0-9a-f]{6}$")
        for color in result:
            self.assertRegex(color, pattern)


if __name__ == "__main__":
    unittest.main()
