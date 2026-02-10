"""Module containing SVG manipulation utilities."""

from __future__ import annotations

import os
from typing import Any

os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/lib"
import cairosvg


class svg_handler:
    """Functions to manipulate SVG elements and save to files."""

    __svg_rect__ = '<rect x="{}" y="{}" width="{}" height="{}" style="stroke-width:1;stroke:{}; fill: {}" />\n'
    __svg_label__ = '<text x="{}" y="{}" text-anchor="{}" font-family="sans-serif" \
        font-size="{}px" fill="{}">{}</text>\n'
    __svg_line__ = (
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}" {} />\n'
    )

    def __init__(
        self: svg_handler,
        svg_string: str,
        width: float,
        height: float,
        background: str | None = None,
    ) -> None:
        """Initialize SVG handler with content and dimensions.

        Args:
            svg_string (str): SVG content string.
            width (float): Width of the SVG.
            height (float): Height of the SVG.
            background (str | None): Background color.
        """
        self.__svg__ = svg_string
        self.__width__ = width
        self.__height__ = height
        self.__background = background

    def group(self: svg_handler, translate: tuple[float, float] = (0, 0)) -> None:
        """Group and transform SVG elements with translation.

        Args:
            translate (tuple[float, float]): Translation coordinates (x, y).
        """
        self.__svg__ = '<g transform="translate(%s %s)">\n%s\n</g>\n' % (
            translate[0],
            translate[1],
            self.__svg__,
        )

        # Updating coordinates:
        self.__width__ += abs(translate[0])
        self.__height__ += abs(translate[1])

    def appendSvg(self: svg_handler, svg_string: str) -> None:
        """Append SVG string to the current content.

        Args:
            svg_string (str): SVG content to append.
        """
        self.__svg__ += svg_string

    def mergeSvg(self: svg_handler, svg_obj: svg_handler) -> None:
        """Merge another SVG handler object into this one.

        Args:
            svg_obj (svg_handler): SVG handler to merge.
        """
        self.__svg__ += svg_obj.getSvg()

        self.__width__ = max(self.__width__, svg_obj.getWidth())
        self.__height__ = max(self.__height__, svg_obj.getHeight())

    def __closeSvg(self: svg_handler) -> None:
        """Close the SVG document by adding header and footer."""
        svg_header = f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" \
            xmlns:xlink="http://www.w3.org/1999/xlink" \
            xml:space="preserve"  width = "{self.__width__}" height = "{self.__height__}" >\n'

        # If background color is set, we use it:
        if self.__background is not None:
            svg_header += (
                f'<rect width="100%" height="100%" fill="{self.__background}" />'
            )

        self.__closedSVG__ = svg_header + self.__svg__ + "\n</svg>\n"

    def savePng(self: svg_handler, filename: str = "test.png") -> None:
        """Save SVG as PNG file.

        Args:
            filename (str): Output filename.
        """
        self.__closeSvg()

        cairosvg.svg2png(bytestring=self.__closedSVG__, write_to=filename)

    def saveSvg(self: svg_handler, filename: str = "test.svg") -> None:
        """Save SVG to file.

        Args:
            filename (str): Output filename.
        """
        self.__closeSvg()

        with open(filename, "w") as f:
            f.write(self.__closedSVG__)

    def getSvg(self: svg_handler) -> str:
        """Return the SVG content string.

        Returns:
            str: SVG content.
        """
        return self.__svg__

    def getWidth(self: svg_handler) -> float:
        """Return the width of the SVG.

        Returns:
            float: Width value.
        """
        return self.__width__

    def getHeight(self: svg_handler) -> float:
        """Return the height of the SVG.

        Returns:
            float: Height value.
        """
        return self.__height__

    def draw_rectangle(
        self: svg_handler,
        x: float,
        y: float,
        width: float,
        height: float,
        stroke: str,
        fill: str,
    ) -> None:
        """Draw a rectangle on the SVG.

        Args:
            x (float): X coordinate.
            y (float): Y coordinate.
            width (float): Rectangle width.
            height (float): Rectangle height.
            stroke (str): Stroke color.
            fill (str): Fill color.
        """
        self.__svg__ += self.__svg_rect__.format(x, y, width, height, stroke, fill)

    def draw_line(
        self: svg_handler,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        stroke: str = "#000000",
        stroke_width: int = 3,
        **kwargs: Any,
    ) -> None:
        """Draw a line on the SVG.

        Args:
            x1 (float): Start X coordinate.
            y1 (float): Start Y coordinate.
            x2 (float): End X coordinate.
            y2 (float): End Y coordinate.
            stroke (str): Stroke color.
            stroke_width (int): Stroke width.
            **kwargs (Any): Additional SVG attributes.
        """
        extra_args = ""

        if kwargs:
            for key, value in kwargs.items():
                extra_args += f' {key.replace("_", "-")}="{value}"'
        print(extra_args)
        self.__svg__ += self.__svg_line__.format(
            x1, y1, x2, y2, stroke, stroke_width, extra_args
        )

    def add_text(
        self: svg_handler,
        x: float,
        y: float,
        text: str,
        size: int | float = 10,
        fill: str = "#000000",
        anchor: str = "start",
    ) -> None:
        """Add text to the SVG.

        Args:
            x (float): X coordinate.
            y (float): Y coordinate.
            text (str): Text content.
            size (int | float): Font size.
            fill (str): Text color.
            anchor (str): Text anchor position.
        """
        self.__svg__ += self.__svg_label__.format(x, y, anchor, size, fill, text)
