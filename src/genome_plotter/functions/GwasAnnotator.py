"""Module for annotating chromosomes with GWAS associations."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


class GwasAnnotator:
    """Adds GWAS associations to the chromosome."""

    # GWAS hit svg definition:
    gwas_hit = (
        '<circle cx="{}" cy="{}" r="{}" stroke="{}" stroke-width="1" fill="{}" />'
    )

    # Unit circle:
    circle_unit = 63

    # GWAS hit count cap:
    gwas_cap = 10

    def __init__(
        self: GwasAnnotator,
        pixel: int,
        chromosome: str,
        gwas_file: str,
        chunk_size: int,
        width: int,
        xoffset: int = 0,
        yoffset: int = 0,
        gwas_color: str = "black",
    ) -> None:
        """Initialize GWAS annotator.

        Args:
            pixel (int): Pixel size for plotting.
            chromosome (str): Chromosome identifier.
            gwas_file (str): Path to GWAS data file.
            chunk_size (int): Size of each genomic chunk.
            width (int): Number of chunks per row.
            xoffset (int): X offset for positioning.
            yoffset (int): Y offset for positioning.
            gwas_color (str): Color for GWAS points.
        """
        gwas_df = pd.read_csv(
            gwas_file,
            compression="gzip",
            sep="\t",
            quotechar='"',
            header=0,
            dtype={"#chr": str, "start": int, "end": int, "rsID": str, "trait": str},
        )

        filtered_locations = gwas_df.loc[gwas_df["#chr"] == chromosome]

        # Count hits per chunk index, capped at gwas_cap:
        chunk_counts = (
            filtered_locations.start.apply(lambda x: int(x / chunk_size))
            .value_counts()
            .clip(upper=self.gwas_cap)
        )

        self.__positions = pd.DataFrame(
            {
                "counts": chunk_counts.values,
                "x": (chunk_counts.index % width).astype(int),
                "y": (chunk_counts.index // width).astype(int),
            }
        )

        self.__pixel = pixel
        self.__xoffset = xoffset
        self.__yoffset = yoffset
        self.__gwas_color = gwas_color

    def generate_gwas(self: GwasAnnotator) -> str:
        """Generate SVG for GWAS hits.

        Returns:
            str: SVG string containing GWAS hit circles.
        """
        pos = self.__positions
        gwas_color = self.__gwas_color

        radius = np.sqrt(pos["counts"] ** 2 * self.circle_unit / math.pi)
        cx = pos["x"] * self.__pixel + radius / 2 + self.__xoffset
        cy = pos["y"] * self.__pixel + radius / 2 + self.__yoffset

        return "\n".join(
            self.gwas_hit.format(cx_val, cy_val, r_val, gwas_color, gwas_color)
            for cx_val, cy_val, r_val in zip(cx, cy, radius)
        )
