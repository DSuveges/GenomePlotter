"""Module containing SVG manipulation utilities."""

from __future__ import annotations

import os
from typing import Any

os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/lib"
import cairosvg


class SvgHandler:
    """Build, compose, and export SVG documents."""

    _RECT = '<rect x="{}" y="{}" width="{}" height="{}" style="stroke-width:1;stroke:{}; fill: {}" />\n'
    _LABEL = '<text x="{}" y="{}" text-anchor="{}" font-family="sans-serif" font-size="{}px" fill="{}">{}</text>\n'
    _LINE = '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}" {} />\n'

    def __init__(
        self,
        svg_string: str,
        width: float,
        height: float,
        background: str | None = None,
    ) -> None:
        """Initialize SVG handler with content and dimensions.

        Args:
            svg_string (str): Initial SVG content.
            width (float): Canvas width.
            height (float): Canvas height.
            background (str | None): Optional background fill color.
        """
        self._svg = svg_string
        self._width = width
        self._height = height
        self._background = background

    def group(self, translate: tuple[float, float] = (0, 0)) -> None:
        """Wrap current content in a <g> translate transform.

        Args:
            translate (tuple[float, float]): (x, y) translation offset.
        """
        self._svg = '<g transform="translate(%s %s)">\n%s\n</g>\n' % (
            translate[0],
            translate[1],
            self._svg,
        )
        self._width += abs(translate[0])
        self._height += abs(translate[1])

    def append_svg(self, svg_string: str) -> None:
        """Append raw SVG content to the current canvas.

        Args:
            svg_string (str): SVG content to append.
        """
        self._svg += svg_string

    def merge_svg(self, other: SvgHandler) -> None:
        """Merge another SvgHandler's content into this one.

        Canvas dimensions are expanded to fit both objects. Callers are
        responsible for positioning other via group() before merging so
        that the max-dimension logic yields the correct total canvas size.

        Args:
            other (SvgHandler): Handler whose content is merged in.
        """
        self._svg += other.get_svg()
        self._width = max(self._width, other.get_width())
        self._height = max(self._height, other.get_height())

    def _build_svg(self) -> str:
        """Return the complete SVG document as a string."""
        header = (
            f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink" '
            f'xml:space="preserve" width="{self._width}" height="{self._height}">\n'
        )
        if self._background is not None:
            header += f'<rect width="100%" height="100%" fill="{self._background}" />'
        return header + self._svg + "\n</svg>\n"

    def save_png(self, filename: str = "test.png") -> None:
        """Render and save the SVG as a PNG file.

        Args:
            filename (str): Output file path.
        """
        cairosvg.svg2png(bytestring=self._build_svg(), write_to=filename)

    def save_svg(self, filename: str = "test.svg") -> None:
        """Save the SVG document to a file.

        Args:
            filename (str): Output file path.
        """
        with open(filename, "w") as f:
            f.write(self._build_svg())

    def get_svg(self) -> str:
        """Return the raw (unwrapped) SVG content string."""
        return self._svg

    def get_width(self) -> float:
        """Return the current canvas width."""
        return self._width

    def get_height(self) -> float:
        """Return the current canvas height."""
        return self._height

    def draw_rectangle(
        self,
        x: float,
        y: float,
        width: float,
        height: float,
        stroke: str,
        fill: str,
    ) -> None:
        """Append a rectangle element.

        Args:
            x (float): X coordinate.
            y (float): Y coordinate.
            width (float): Rectangle width.
            height (float): Rectangle height.
            stroke (str): Border color.
            fill (str): Fill color.
        """
        self._svg += self._RECT.format(x, y, width, height, stroke, fill)

    def draw_line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        stroke: str = "#000000",
        stroke_width: int = 3,
        **kwargs: Any,
    ) -> None:
        """Append a line element.

        Args:
            x1 (float): Start X.
            y1 (float): Start Y.
            x2 (float): End X.
            y2 (float): End Y.
            stroke (str): Line color.
            stroke_width (int): Line width.
            **kwargs: Additional SVG attributes (underscores become hyphens).
        """
        extra = "".join(
            f' {key.replace("_", "-")}="{value}"' for key, value in kwargs.items()
        )
        self._svg += self._LINE.format(x1, y1, x2, y2, stroke, stroke_width, extra)

    def add_text(
        self,
        x: float,
        y: float,
        text: str,
        size: int | float = 10,
        fill: str = "#000000",
        anchor: str = "start",
    ) -> None:
        """Append a text element.

        Args:
            x (float): X coordinate.
            y (float): Y coordinate.
            text (str): Text content.
            size (int | float): Font size in pixels.
            fill (str): Text color.
            anchor (str): SVG text-anchor value.
        """
        self._svg += self._LABEL.format(x, y, anchor, size, fill, text)
