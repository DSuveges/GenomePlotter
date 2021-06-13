import unittest

import pandas as pd

from functions.color_fun import linear_gradient, hex_to_RGB, RGB_to_hex, color_darkener

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

    def test_hex_to_RGB(self):

        # Test for good output:
        hex_col = '#000000'
        rgb_col = hex_to_RGB(hex_col)

        self.assertIsInstance(rgb_col, list)
        self.assertEqual(len(rgb_col), 3)
        for i in rgb_col:
            self.assertIsInstance(i, int)
        self.assertEqual([0, 0, 0], rgb_col)

        # Test for another good output:
        hex_col = '#ffffff'
        self.assertEqual([255, 255, 255], hex_to_RGB(hex_col))

        # Testing for bad output:
        for bad_input in ['cica', True, 13, '#209345209', 'ffffff']:
            with self.assertRaises(ValueError):
                hex_to_RGB(bad_input)

    def test_rgb_to_hex(self):

        # Test for good output:
        rgb_col = [0, 0, 0]
        hex_col = RGB_to_hex(rgb_col)

        self.assertIsInstance(hex_col, str)
        self.assertEqual(hex_col[0], '#')

    def test_color_darkener(self):

        # Basic correct input:
        row = pd.DataFrame({
            'x': [145, 100],
            'color': ['#FFFFFF', '#123456'],
        })
        width = 200
        threshold = 0.5
        max_diff_value = 0.9

        # Testing for input type:
        with self.assertRaises(TypeError):
            color_darkener('cicaful', width=width, threshold=threshold, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(row, width='pocok', threshold=threshold, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(row, width=width, threshold=3.0, max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(row, width=width, threshold='pocok', max_diff_value=max_diff_value)
        with self.assertRaises(TypeError):
            color_darkener(row, width=width, threshold=threshold, max_diff_value=1232)

        df = pd.DataFrame({
            'x': [45, 130, 200],
            'color': ['#FFFFFF', '#123456', '#234566'],
        })

        # Testing if the threshold is appreciated:
        row = df.iloc[0]
        self.assertEqual(color_darkener(row, width, threshold, max_diff_value), row['color'])
        row = df.iloc[1]
        self.assertNotEqual(color_darkener(row, width, threshold, max_diff_value), row['color'])

        # Testing if the diff value is appreciated:
        row = df.iloc[2]
        self.assertEqual(color_darkener(row, width, threshold, 1.0), '#000000')

if __name__ == '__main__':
    unittest.main()
