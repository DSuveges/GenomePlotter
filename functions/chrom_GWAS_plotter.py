"""Module for processing chromosome GWAS plots."""

from __future__ import annotations

import gzip

import pandas as pd


class process_chrom:
    """This class plots the chromosomes and saves a gzipped text file."""

    defs: dict[str, str]
    plot: str

    def __init__(
        self: process_chrom, width: int, height: int, pixel: int, defs: int = 0
    ) -> None:
        """Initialize the chromosome processor.

        Args:
            width (int): Width of the plot.
            height (int): Height of the plot.
            pixel (int): Pixel size.
            defs (int): Whether to include defs tag.
        """
        self.width = width
        self.height = height
        self.pixel = pixel
        self.plot = ""

        # If we would like, we can include the defs tag as well:
        if defs == 1:
            self.__add_def(
                "chunk",
                '<rect x="0" y="0" width="%s" height="%s" style="stroke-width:0" />\n'
                % (pixel, pixel),
            )

    def __add_def(self: process_chrom, ID: str, svg_string: str) -> None:
        """Add a <g> object to the <defs> section of the SVG.

        Args:
            ID (str): Definition identifier.
            svg_string (str): SVG definition string.
        """
        if hasattr(self, "defs"):
            self.defs[ID] = svg_string
        else:
            self.defs = {}
            self.defs[ID] = svg_string

    def draw_chunk(self: process_chrom, row: pd.Series) -> None:
        """Draw a chunk on the plot.

        Args:
            row (pd.Series): Row with x, y, and color columns.
        """
        x = row["x"]
        y = row["y"]
        color = row["color"]
        self.plot += f'<use x="{x*self.pixel}" y="{y*self.pixel}" href="#chunk" style="fill: {color};"/>\n'

    def draw_GWAS(self: process_chrom, row: pd.Series) -> None:
        """Draw a GWAS signal as a dot.

        Args:
            row (pd.Series): Row with x and y columns.
        """
        x = row["x"]
        y = row["y"]
        self.plot += f'<circle cx="{x*self.pixel}" cy="{y*self.pixel}" r="2" stroke="black" stroke-width="1" fill="black" />\n'

    def mark_centromere(
        self: process_chrom, centr_start: int, centr_end: int
    ) -> None:
        """Compile SVG path for centromere and add to plot.

        Args:
            centr_start (int): Centromere start position.
            centr_end (int): Centromere end position.
        """
        # The y coordinate of the centromere start and end points 0 adjusted:
        centr_y_start = 0
        centr_y_end = self.pixel * (centr_end - centr_start) / self.width
        centr_y_midpoint = centr_y_end / 2
        centr_x_midpoint = centr_y_end * 2
        print(centr_y_end, centr_y_start)

        # Creating the centromere path:
        half_centromere = (
            '<path d="M %s %s C %s %s, %s %s, %s %s C %s %s, %s %s, %s %s Z" fill="white"/>'
            % (
                0,
                centr_y_start,
                0,
                centr_y_midpoint,
                centr_x_midpoint / 2,
                centr_y_midpoint,
                centr_x_midpoint,
                centr_y_midpoint,
                centr_x_midpoint / 2,
                centr_y_midpoint,
                0,
                centr_y_midpoint,
                0,
                centr_y_end,
            )
        )

        # Set position for the left side:
        left_centromere_side = '<g transform="translate(0, %s)">\n\t%s\n</g>\n' % (
            (centr_start * self.pixel) / self.width,
            half_centromere,
        )

        # Set position for the right side:
        right_centromere_side = (
            '<g transform="rotate(180 0 %s) translate(-%s, -%s)">\n\t%s\n</g>\n'
            % (
                centr_y_midpoint,
                self.pixel * self.width,
                (centr_start * self.pixel) / self.width,
                half_centromere,
            )
        )

        # Adding centromeres to the plot:
        self.plot += '<g id="centromere">\n%s %s </g>\n' % (
            left_centromere_side,
            right_centromere_side,
        )

    def save(self: process_chrom, chr_name: str) -> None:
        """Save the assembled lines into a gzipped text file.

        Args:
            chr_name (str): Chromosome name for the output file.
        """
        with gzip.open(f"chr{chr_name}_genome_chunks.txt.gz", "wt") as f:
            f.write(self.plot)
