import colorsys
import pandas as pd
import numpy as np
import logging
import re

def hex_to_rgb(hex_color: str) -> list:
    """Converting hexadecimal color to rgb

    params:
        hex_color (str): color represented in hexadecimal format eg '#FFFFFF'

    returns:
        list: RGB color representation, list with 3 integers eg [255, 255, 255]
    """

    if not isinstance(hex_color, str):
        raise ValueError('The provided hexadecimal definition has to be string eg #000000.')
    elif len(hex_color) != 7 or not hex_color.startswith('#'):
        raise ValueError('The provided hexadecimal definition has to starts with #')

    hex_color = hex_color.lower()

    # Pass 16 to the integer function for change of base
    return [int(hex_color[i:i + 2], 16) for i in range(1, 6, 2)]


def rgb_to_hex(rgb_color: list) -> str:
    """Converting rgb color to hexadecimal representation

    params:
        rgb_color (list): RGB color representation, list with 3 integers

    returns:
        str: color represented in hexadecimal values eg. '#ffffff'
    """
    # Components need to be integers for hex to make sense
    rgb_color = [int(x) for x in rgb_color]
    return "#" + "".join(
        ["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in rgb_color]
    )


def linear_gradient(start_hex: str, finish_hex: str = "#FFFFFF", length: int = 10) -> list:
    """Generating color gradient between two hexadecimal color of a given length

    Params:
        start_hex (str): starting color in hexadecimal format, requried
        finish_hex (str): ending color in hexadecimal format, default: '#FFFFFF'
        length (int): number of colors in the gradient

    Returns:
        list: 'length' number of colors in hexadecimal format
    """
    # Starting and ending colors in RGB form
    start_rgb = hex_to_rgb(start_hex)
    finish_rgb = hex_to_rgb(finish_hex)

    # Initilize a list of the output colors with the starting color
    rgb_list = [start_hex.lower()]

    if not isinstance(length, int):
        raise ValueError('The number of returned colors have to be specified by an integer.')
    elif length == 0:
        return []

    # Calcuate a color at each evenly spaced value of t from 1 to n
    for step in range(1, length):

        # Interpolate RGB vector for color at the current value of t
        curr_vector = [
            int(start_rgb[j] + (float(step) / (length - 1)) * (finish_rgb[j] - start_rgb[j]))
            for j in range(3)
        ]

        # Add it to our list of output colors
        rgb_list.append(rgb_to_hex(curr_vector))

    return rgb_list


def color_darkener(color: str, x: int, width: int, threshold: float, max_diff_value: float) -> str:
    """Function to decrease the luminosity of a hex color

    Params:
        color (str): primary color assigned based on primary annotation eg. intergenic, exon etc,
        x (int): x position of the chunk
        width (str): how many chunks do we have in one line. (always >= row['x'])
        threshold (float): fraction of the width, where the darkening starts (<= 1.0)
        max_diff_value (float): the max value of darkening (<=1)

    Returns:
        str: the darkness adjusted color in hex
    """

    if not isinstance(color, str) or not re.match(r'#[0-9A-F]{6}', color):
        raise TypeError('Color is expected to be given as a hexadecimal value (eg. "#F12AC4").')

    if not isinstance(x, int):
        raise TypeError('x has to be integer giving the numeric position of the chunk.')

    if not isinstance(width, int):
        raise TypeError('width has to be an integer.')

    if not isinstance(threshold, float) or (threshold > 1):
        raise TypeError('threshold has to be a float below 1. The fraction of the width where the darkening starts')

    if not isinstance(max_diff_value, float) or (max_diff_value > 1):
        raise TypeError('The darkening has to be a float below 1.')

    col_frac = x / width

    # We have the color, now based on the column we make it a bit darker:
    if col_frac > threshold:
        diff = (col_frac - threshold) / (1 - threshold)
        factor = 1 - max_diff_value * diff

        # Get rgb code of the hexa code:
        rgb_code = hex_to_rgb(color)

        # Get the hls code of the rgb:
        hls_code = colorsys.rgb_to_hls(
            rgb_code[0] / 255,
            rgb_code[1] / 255,
            rgb_code[2] / 255
        )

        # Scaling luminosity then convert to RGB:
        new_rgb = colorsys.hls_to_rgb(hls_code[0],
                                      hls_code[1] * factor,
                                      hls_code[2])

        # Get the modifed hexacode:
        color = rgb_to_hex([x * 255 for x in new_rgb])

    return(color)


class ColorPicker(object):

    # These are the supported and expected features:
    features = ['exon', 'gene', 'intergenic', 'centromere', 'heterochromatin', 'dummy']

    def __init__(self, colors: dict, dark_max: float = None, dark_threshold: float = None,
                 count: int = 20, width: int = None) -> None:
        """Initializing ColorPicker

        Params:
            colors (dict):
            dark_max (float): How much darker color should be reached.
            dark_threshold (float): Where the darkenin should be started.
            count (int): Number of steps in the gradient.
            width (int): Number of chunks in a row.

        Returns:
            None
        """

        # Checking if colors is of the good type:
        if not isinstance(colors, dict):
            raise TypeError(f'The color object has to be a dictionary. {type(colors)} is given.')

        # Checking if all features can be found in the color set:
        if not pd.Series(self.features).isin(list(colors.keys())).all():
            raise ValueError(f'The following keys must be defined in the color sets: {",".join(self.features)}')

        # Checking if all the values are good:
        for hex_color in colors.values():
            if not re.match(r'#[0-9A-F]{6}', hex_color) or not isinstance(hex_color, str):
                raise ValueError('All colors should be in hexadecimal format eg "#1ED5FA"')

        # Checking dark max and the threshold:
        for val in [dark_max, dark_threshold]:
            if (not isinstance(val, float)) or (val > 1) or (val < 0):
                raise ValueError(f'dark_max and dark_threshold has to be float between 0 and 1 instead: {val}')

        # Checking dark max and the threshold:
        for val in [count, width]:
            if not isinstance(val, int):
                raise ValueError(f'Count and width have to be integers. Got: {val}')

        # Generating color gradients for a given length:
        self.color_map = {x: linear_gradient(colors[x], length=count) for x in self.features}

        self.dark_max = dark_max
        self.dark_threshold = dark_threshold
        self.width = width
        self.count = count

    def map_color(self, feature: str, gc_content: float) -> str:

        if feature == 'dummy':
            color = self.color_map['dummy'][0]
        elif gc_content is None or np.isnan(gc_content):
            color = self.color_map['heterochromatin'][0]
        else:
            try:
                color = self.color_map[feature][int(gc_content * (self.count - 1))]
            except KeyError:
                logging.error(f'Feature {feature} was not found in color mapper. Returning black.')
                color = '#000000'

        return color

    def pick_color(self, row: pd.Series) -> str:

        expected_columns = ['GC_ratio', 'GENCODE', 'x']
        if not isinstance(row, pd.Series) or not pd.Series(expected_columns).isin(row.index).all():
            raise TypeError(f'The row has to be a pd.Series with the following keys: {",".join(expected_columns)}')

        # Get the base color:
        color = self.map_color(row['GENCODE'], row['GC_ratio'])

        # Darken the color:
        if (row['GENCODE'] != 'dummy') and self.width is not None and (row['x'] / self.width) > self.dark_threshold:
            color = color_darkener(color, row['x'], self.width, self.dark_threshold, self.dark_max)

        return color
