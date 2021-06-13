import colorsys
import pandas as pd

def hex_to_RGB(hex_color: str) -> list:
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

def RGB_to_hex(rgb_color: list) -> str:
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
    start_rgb = hex_to_RGB(start_hex)
    finish_rgb = hex_to_RGB(finish_hex)

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
        rgb_list.append(RGB_to_hex(curr_vector))

    return rgb_list

def color_darkener(row: pd.Series, width: int, threshold: float, max_diff_value: float) -> str:
    """Function to decrease the luminosity of a hex color

    Params:
        row (pd.Series): one row of the genome dataframe with the following indexes:
            'color': primari color assigned based on primary annotation eg. intergenic, exon etc,
            'x': x position of the given chunk.
        width (str): how many chunks do we have in one line. (always >= row['x'])
        threshold (float): fraction of the width, where the darkening starts (<= 1.0)
        max_diff_value (float): the max value of darkening (<=1)

    Returns:
        str: the darkness adjusted color in hex
    """

    if not isinstance(row, pd.Series) or ('x' not in row) or ('color' not in row):
        raise TypeError('row has to be a ps.Series object with indices "color" and "x".')

    if not isinstance(width, int):
        raise TypeError('width has to be an integer.')

    if not isinstance(threshold, float) or (threshold > 1):
        raise TypeError('threshold has to be a float below 1. The fraction of the width where the darkening starts')

    if not isinstance(max_diff_value, float) or (max_diff_value > 1):
        raise TypeError('The darkening has to be a float below 1.')

    color = row['color']
    col_frac = row['x'] / width

    # Ok, we have the color, now based on the column we make it a bit darker:
    if col_frac > threshold:
        diff = (col_frac - threshold) / (1 - threshold)
        factor = 1 - max_diff_value * diff

        # Get rgb code of the hexa code:
        rgb_code = hex_to_RGB(color)

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
        color = RGB_to_hex([x * 255 for x in new_rgb])

    return(color)
