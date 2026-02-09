"""Module for annotating genes on chromosome plots."""

from __future__ import annotations

import logging

import pandas as pd


class PositionConverter:
    """Helper class to convert genomic position to plot y coordinate."""

    def __init__(
        self: PositionConverter, width: int, chunkSize: int, pixel: int
    ) -> None:
        """Initialize the position converter.

        Args:
            width (int): Number of chunks per row.
            chunkSize (int): Size of each chunk.
            pixel (int): Pixel size.
        """
        self.width = width
        self.chunkSize = chunkSize
        self.pixel = pixel

    def convert(self: PositionConverter, position: int) -> float:
        """Convert genomic position to plot coordinate.

        Args:
            position (int): Genomic position.

        Returns:
            float: Plot y coordinate.
        """
        return position / (self.width * self.chunkSize) * self.pixel


class GeneAnnotator:
    """Adds gene annotations to chromosomes."""

    # Text line:
    text = '<text x="{x}" y="{y}" alignment-baseline="baseline" \
        text-anchor="start" font-family="sans-serif" font-size="{font_size}px" fill="{font_color}">{gene_name}</text>\n'

    # Segment line:
    segment = '<line x1="{x1}" x2="{x2}" y1="{y1}" y2="{y2}" stroke="black" stroke-width="2"/>\n'

    # Font color:
    __fontColor = "black"

    def __init__(
        self: GeneAnnotator,
        gene_file: str,
        centromerePosition: int,
        chromosome: str,
        chunk_size: int,
        pixel: int,
        height: int | float,
        width: int,
    ) -> None:
        """Initialize the gene annotator.

        Args:
            gene_file (str): Path to compressed tab-delimited gene file.
            centromerePosition (int): Position of the centromere.
            chromosome (str): Chromosome identifier.
            chunk_size (int): Size of each chunk.
            pixel (int): Pixel size for plotting.
            height (int): Plot height.
            width (int): Number of chunks per row.

        Note:
            Expected columns in gene_file: chr, start, end, name.
        """
        # Open gene file and store data:
        df = pd.read_csv(gene_file, sep="\t", compression="infer")
        self.__df = df.loc[df.chr == chromosome]
        logging.info(
            "Adding annotation for {} genes on this chromosome.".format(len(self.__df))
        )

        # Store other values:
        self.__chunkSize = chunk_size
        self.__pixel = pixel
        self.__svg_top: float = 0
        self.__svg_bottom: float = height
        self.__width = width

        self.__fontSize = pixel * 10
        self.__centromerePosition = centromerePosition
        self.__chromosome = chromosome

        # Initialize return values:
        self.__gene_annot = ""

    def genes_on_chromosome_arms(
        self: GeneAnnotator, df: pd.DataFrame, arm: str
    ) -> None:
        """Annotate genes on a chromosome arm.

        Args:
            df (pd.DataFrame): DataFrame with gene data.
            arm (str): Chromosome arm ('p' or 'q').
        """
        # Initialize position converter:
        converter = PositionConverter(self.__width, self.__chunkSize, self.__pixel)

        # Get centromere position:
        limit = converter.convert(self.__centromerePosition)

        # Extract values:
        font_size = self.__fontSize
        font_color = self.__fontColor

        for _, row in df.iterrows():
            pos = converter.convert(row["start"])

            if arm == "p":
                if limit <= pos:
                    text_pos = limit - 10
                    limit = text_pos - font_size - 10
                else:
                    text_pos = pos
                    limit = text_pos - font_size - 10
            if arm == "q":
                if limit > pos:
                    text_pos = limit + 10
                    limit = text_pos + font_size + 10
                else:
                    text_pos = pos
                    limit = text_pos + font_size + 10

            # Update svg_top and bottom if required:
            if limit - 50 < self.__svg_top:
                self.__svg_top = limit - font_size - 50

            if limit + 50 > self.__svg_bottom:
                self.__svg_bottom = limit + font_size + 50

            # Add segments:
            self.__gene_annot += self.segment.format(
                **{"x1": 0, "x2": 100, "y1": pos, "y2": pos}
            )
            self.__gene_annot += self.segment.format(
                **{"x1": 100, "x2": 200, "y1": pos, "y2": text_pos}
            )
            self.__gene_annot += self.segment.format(
                **{"x1": 200, "x2": 300, "y1": text_pos, "y2": text_pos}
            )

            # Add text:
            self.__gene_annot += self.text.format(
                **{
                    "x": 200,
                    "y": text_pos - 10,
                    "font_size": font_size,
                    "font_color": font_color,
                    "gene_name": row["name"],
                }
            )

    def generate_gene_annotation(self: GeneAnnotator) -> None:
        """Generate SVG annotations for all genes."""
        # Initialize gene annotation:
        self.__gene_annot = ""

        # Generate text for the p-arm:
        self.genes_on_chromosome_arms(
            self.__df.loc[self.__df.start <= self.__centromerePosition].sort_values(
                by=["start"], ascending=False
            ),
            "p",
        )

        # Generate text for the q-arm:
        self.genes_on_chromosome_arms(
            self.__df.loc[self.__df.start > self.__centromerePosition].sort_values(
                by=["start"], ascending=True
            ),
            "q",
        )

    def get_annotation(self: GeneAnnotator) -> str:
        """Return the gene annotation SVG string.

        Returns:
            str: SVG content string.
        """
        return self.__gene_annot

    def get_dimensions(self: GeneAnnotator) -> tuple[float, float]:
        """Get the dimensions of the annotation plot.

        Returns:
            tuple[float, float]: Top and bottom coordinates.
        """
        return (self.__svg_top, self.__svg_bottom)
