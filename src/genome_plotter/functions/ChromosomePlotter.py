"""Module for plotting chromosome data as SVG."""

from __future__ import annotations

import os

os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/lib"

import cairosvg
import pandas as pd
from loguru import logger


class ChromosomePlotter:
    """Functions to store, process and save SVG object of a single chromosome."""

    chunk_svg = '<rect x="{}" y="{}" width="{}" height="{}" style="stroke-width:1;stroke:{}; fill: {}" />'

    def __init__(self, input_data: pd.DataFrame, pixel: int) -> None:
        """Initialize plotter object.

        Args:
            input_data (pd.DataFrame): DataFrame with chromosome data.
            pixel (int): Size of pixel (one unit of genetic information).
        """
        self._pixel = pixel
        self._chromosome_data = input_data
        self._chromosome_name = input_data["chr"].iloc[0]
        self._chunk_size = input_data.iloc[0].end

        self._width = pixel * (input_data.x.max() + 1)
        self._height = pixel * (input_data.y.max() + 1)

        self._plot_string = ""

    def _add_centromere(self) -> None:
        """Add centromere visualization to the chromosome plot."""
        if "centromere" not in self._chromosome_data.GENCODE.to_list():
            return

        centromere_rows = self._chromosome_data.loc[
            self._chromosome_data.GENCODE == "centromere"
        ]
        centromere_start = centromere_rows.y.min() * self._pixel
        centromere_end = centromere_rows.y.max() * self._pixel

        centromere_midpoint = (centromere_end - centromere_start) / 2
        centromere_height = centromere_end - centromere_start
        centromere_x = self._width / 3

        centromere_string = (
            f'<path d="M -1 0 '
            f'C 0 {centromere_midpoint}, {centromere_x / 2} {centromere_midpoint}, {centromere_x / 2} {centromere_midpoint} '
            f'C {centromere_x / 2} {centromere_midpoint}, 0 {centromere_midpoint}, -1 {centromere_height} Z" fill="white"/>\n'
        )
        half_centromere = f'<g transform="translate(0, {centromere_start})">\n\t{centromere_string}\n</g>\n'
        other_half = (
            f'<g transform="rotate(180 0 {centromere_start + centromere_midpoint}) '
            f'translate(-{self._width}, 0)">\n\t{half_centromere}\n</g>\n'
        )
        self._plot_string += f'\n<g id="centromere">\n\t{half_centromere}\t{other_half}</g>\n'

    def draw_dummy(self) -> None:
        """Draw a dummy chromosome representation."""
        width = self._width
        height = self._height
        data = self._chromosome_data

        dummy_color = data.loc[data.GENCODE != "centromere"].color.tolist()[0]
        centromere_color = data.loc[data.GENCODE == "centromere"].color.tolist()[0]

        centromere_rows = data.loc[data.GENCODE == "centromere"]
        centromere_start = centromere_rows.y.min() * self._pixel
        centromere_end = centromere_rows.y.max() * self._pixel - centromere_start

        logger.info(f"centromere_start: {centromere_start}, centromere_end: {centromere_end}")

        self._plot_string += self.chunk_svg.format(0, 0, width, height, dummy_color, dummy_color)
        self._plot_string += self.chunk_svg.format(
            0, centromere_start, width, centromere_end, centromere_color, centromere_color
        )
        self._add_centromere()

    def draw_chromosome(self) -> None:
        """Draw the chromosome with colored chunks."""
        pixel = self._pixel
        data = self._chromosome_data
        xs = data["x"].to_numpy() * pixel
        ys = data["y"].to_numpy() * pixel
        colors = data["color"].to_numpy()
        self._plot_string = "\n".join(
            self.chunk_svg.format(x, y, pixel, pixel, c, c)
            for x, y, c in zip(xs, ys, colors)
        )
        self._add_centromere()

    def get_plot_with(self) -> int:
        """Return the plot width in pixels."""
        return self._width

    def get_plot_height(self) -> int:
        """Return the plot height in pixels."""
        return self._height

    def return_svg(self) -> str:
        """Return the SVG plot string."""
        return self._plot_string

    def _wrap_svg_string(self) -> str:
        return (
            '<svg width="%s" height="%s" version="1.1" xmlns="http://www.w3.org/2000/svg" '
            'xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">\n%s</svg>'
            % (self._width, self._height, self._plot_string)
        )

    def save_png(self, file_name: str) -> None:
        """Save the plot as a PNG file.

        Args:
            file_name (str): Output file path.
        """
        cairosvg.svg2png(bytestring=self._wrap_svg_string(), write_to=file_name)

    def wrap_svg(self, file_name: str) -> None:
        """Wrap the SVG content and save to file.

        Args:
            file_name (str): Output file path.
        """
        with open(file_name, "w") as f:
            f.write(self._wrap_svg_string())
