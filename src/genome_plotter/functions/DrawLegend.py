#!/usr/bin/env python
"""Module for drawing legends for chromosome plots."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict
from typing import TYPE_CHECKING

from genome_plotter.functions.ColorFunctions import linear_gradient
from genome_plotter.functions.svg_handler import svg_handler

if TYPE_CHECKING:
    from genome_plotter.functions.ConfigManager import Config


class DrawLegend:
    """Class for drawing legends for chromosome plots.

    Note:
        This class is under active development.
        TODO: Integrate datasource versions.
        TODO: Add example gene.
    """

    def __init__(self: DrawLegend, config_manager: Config) -> None:
        """Initialize the legend drawer.

        Args:
            config_manager (Config): Configuration manager object.
        """
        self.__colors__ = asdict(config_manager.color_schema)
        self.__pixel__ = config_manager.plot_parameters.pixel_size
        self.__config_manager__ = config_manager

    def draw_colors_legend(self: DrawLegend) -> None:
        """Draw a legend describing colors of the chromosome plot."""
        # These values and the order is hardcoded:
        regions = OrderedDict(
            {
                "exon": "Exons",
                "gene": "Introns",
                "centromere": "Centromere",
                "heterochromatin": "Heterochromatic region",
            }
        )

        # Extract variables:
        colors = self.__colors__
        pixel = self.__pixel__ * 5

        # Initialize vertical position:
        y_position: float = 0

        # Initialize svg object:
        self.__svg_color__ = svg_handler("", 750, 300)

        # Generating the color gradients:
        for name, label in regions.items():
            x_position: float = 0
            start_color = colors[name]
            end_color = colors[name] if name == "heterochromatin" else "#FFFFFF"
            gradient = linear_gradient(start_color, end_color, length=20)

            # Drawing boxes:
            for color in gradient:
                self.__svg_color__.draw_rectangle(
                    x_position, y_position, pixel / 2, pixel, color, color
                )
                x_position += pixel / 2

            # Adding label for the box:
            self.__svg_color__.add_text(
                x_position + pixel * 0.5, y_position + pixel * 0.6, label, pixel / 2
            )
            y_position += pixel * 1.1

        # Adding percent values:
        for x_perc in [0, 50, 100]:
            x1 = 0 + x_perc * pixel / 10
            y1 = 0
            y2 = 4 * pixel * 1.1 + 10

            # Adding vertical line:
            self.__svg_color__.draw_line(
                x1, y1, x1, y2, stroke_width=2, stroke_dasharray="4"
            )

            # Adding percent mark:
            self.__svg_color__.add_text(
                x1 - pixel * 0.5, y2 + pixel * 0.5, f"{x_perc}%", pixel * 0.5
            )

        # Adding GC content label:
        self.__svg_color__.add_text(4 * pixel, 5.8 * pixel, "GC-content", pixel * 0.5)

        # Group together:
        self.__svg_color__.group(translate=(30, 100))

        # Adding label:
        self.__svg_color__.add_text(200, 50, "Chromosome plot colors", pixel * 0.8)

    def save(self: DrawLegend, filename: str) -> None:
        """Save the legend as PNG.

        Args:
            filename (str): Output file name.
        """
        self.__svg_color__.savePng(filename)
