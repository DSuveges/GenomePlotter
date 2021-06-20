import unittest
import re

import pandas as pd

from functions.ColorFunctions import linear_gradient, hex_to_rgb, rgb_to_hex, color_darkener, ColorPicker

class TestColorFunctions(unittest.TestCase):
    def test_linear_gradient(self):

        # Testing the default lenght of the gradient:
        gradient = linear_gradient('#000000', '#FFFFFF')
        self.assertEqual(len(gradient), 10)

        # Testing a custom lenght of the gradient:
        gradient = linear_gradient('#000000', '#FFFFFF', 15)
        self.assertEqual(len(gradient), 15)

        # Testing a custom lenght of the gradient:
        gradient = linear_gradient('#000000', '#FFFFFF', 0)
        self.assertEqual(len(gradient), 0)

        # Testing a custom lenght of the gradient:
        with self.assertRaises(ValueError):
            gradient = linear_gradient('#000000', '#FFFFFF', 'cica')

        # Testing if some weird stuff is going on with the input:
        with self.assertRaises(ValueError):
            gradient = linear_gradient('#cica', '#FFFFFF')

        with self.assertRaises(ValueError):
            gradient = linear_gradient('#000000', '#cica')

    def test_hex_to_rgb(self):

        # Test for good output:
        hex_col = '#000000'
        rgb_col = hex_to_rgb(hex_col)

        self.assertIsInstance(rgb_col, list)
        self.assertEqual(len(rgb_col), 3)
        for i in rgb_col:
            self.assertIsInstance(i, int)
        self.assertEqual([0, 0, 0], rgb_col)

        # Test for another good output:
        hex_col = '#ffffff'
        self.assertEqual([255, 255, 255], hex_to_rgb(hex_col))

        # Testing for bad output:
        for bad_input in ['cica', True, 13, '#209345209', 'ffffff']:
            with self.assertRaises(ValueError):
                hex_to_rgb(bad_input)

    def test_rgb_to_hex(self):

        # Test for good output:
        rgb_col = [0, 0, 0]
        hex_col = rgb_to_hex(rgb_col)

        self.assertIsInstance(hex_col, str)
        self.assertEqual(hex_col[0], '#')

    def test_color_darkener(self):

        # Basic correct input:
        color = '#DDDDDD'
        x = 234
        width = 200
        threshold = 0.5
        max_diff_value = 0.9

        # Testing for input type:
        with self.assertRaises(TypeError):
            color_darkener('cicaful', x, width=width, threshold=threshold, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(color, x, width='pocok', threshold=threshold, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(color, x, width=width, threshold=3.0, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(color, x, width=width, threshold='pocok', max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(color, x, width=width, threshold=threshold, max_diff_value=1232)

        # Testing if the threshold is appreciated:
        self.assertEqual(color_darkener(color, 50, width, threshold, max_diff_value), color)
        self.assertNotEqual(color_darkener(color, x, width, threshold, max_diff_value), color)

    def test_color_picker(self):

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
            ColorPicker('Pocok', dark_max, dark_start, count, width)
        with self.assertRaises(ValueError):
            ColorPicker({'color': 'color'}, dark_max, dark_start, count, width)
        with self.assertRaises(ValueError):
            color_map_wrong = color_map.copy()
            color_map_wrong['exon'] = 'cica'
            ColorPicker(
                colors=color_map_wrong, dark_max=dark_max, dark_threshold=dark_start, count=count, width=width
            )
        with self.assertRaises(ValueError):
            ColorPicker(color_map, 'aaaa', dark_start, count, width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, 12.0, dark_start, count, width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, -0.2, count, width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, dark_start, 'pocok', width)
        with self.assertRaises(ValueError):
            ColorPicker(color_map, dark_max, dark_start, count, 'foo')

        # Get correct object initialized:
        cp = ColorPicker(color_map, dark_max, dark_start, count, width)

        self.assertIsInstance(cp.color_map, dict)
        self.assertTrue(pd.Series(list(color_map.keys())).isin(list(cp.color_map.keys())).all())

        for feature in cp.color_map.keys():
            self.assertIsInstance(cp.color_map[feature], list)
            self.assertEqual(len(cp.color_map[feature]), count)

        # So far, so good. Testing functionality:
        self.assertIsInstance(cp.map_color('dummy', 0.5), str)
        color = cp.map_color('dummy', 0.5)
        self.assertTrue(re.match(r'#[0-9a-f]{6}', color))
        self.assertEqual(cp.map_color('dummy', 0.5), color_map['dummy'].lower())

        # Do we get heterochromatin:
        self.assertEqual(cp.map_color('exon', None), color_map['heterochromatin'].lower())

        # No keys for unknown feature:
        self.assertEqual(cp.map_color('cicaful', 0.5), '#000000')

        # Let's test mapper:
        with self.assertRaises(TypeError):
            cp.pick_color('pocok')
        with self.assertRaises(TypeError):
            cp.pick_color(pd.Series({'GC_ratio': 0.3, 'GENCODE': 'exon'}))

        self.assertEqual(
            cp.pick_color(pd.Series({'GC_ratio': 0.3, 'GENCODE': 'dummy', 'x': 150})),
            color_map['dummy'].lower()
        )

        self.assertEqual(
            cp.pick_color(pd.Series({'GC_ratio': None, 'GENCODE': 'exon', 'x': 150})),
            color_map['heterochromatin'].lower()
        )

if __name__ == '__main__':
    unittest.main()
