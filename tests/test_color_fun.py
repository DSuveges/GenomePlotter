import unittest

from functions.color_fun import linear_gradient, hex_to_RGB, RGB_to_hex

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

if __name__ == '__main__':
    unittest.main()
