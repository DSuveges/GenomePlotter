import colorsys

def hex_to_RGB(hex_color):
    """ "#FFFFFF" -> [255,255,255] """

    if not isinstance(hex_color, str):
        raise ValueError('The provided hexadecimal definition has to be string eg #000000.')
    elif len(hex_color) != 7 or not hex_color.startswith('#'):
        raise ValueError('The provided hexadecimal definition has to starts with #')

    hex_color = hex_color.lower()

    # Pass 16 to the integer function for change of base
    return [int(hex_color[i:i + 2], 16) for i in range(1, 6, 2)]

def RGB_to_hex(RGB):
    """ [255,255,255] -> "#FFFFFF" """
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#" + "".join(
        ["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB]
    )

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    """
    returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the # sign (eg "#FFFFFF")
    """
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)

    # Initilize a list of the output colors with the starting color
    RGB_list = [start_hex.lower()]

    if not isinstance(n, int):
        raise ValueError('The number of returned colors have to be specified by an integer.')
    elif n == 0:
        return []

    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):

        # Interpolate RGB vector for color at the current value of t
        curr_vector = [
            int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j]))
            for j in range(3)
        ]

        # Add it to our list of output colors
        RGB_list.append(RGB_to_hex(curr_vector))

    return RGB_list

def color_darkener(row, width, threshold, max_diff_value):
    """
    width - how many chunks do we have in one line.
    threshold - at which column the darkening starts
    max_diff_value - the max value of darkening

    Once the colors are assigned, we make the colors darker for those columns
    that are at the end of the plot.
    """

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
