"""Module for annotating chromosomes with GWAS associations."""

from __future__ import annotations

import math

import pandas as pd


class gwas_annotator:
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
        self: gwas_annotator,
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
        # Reading gwas file:
        gwas_df = pd.read_csv(
            gwas_file,
            compression="gzip",
            sep="\t",
            quotechar='"',
            header=0,
            dtype={"#chr": str, "start": int, "end": int, "rsID": str, "trait": str},
        )

        # Filtering dataframe for the given chromosome:
        filtered_locations = gwas_df.loc[gwas_df["#chr"] == chromosome]

        # Looping through all GWAS hits on the chromosome and calculate coordinates and scale:
        data = []
        for index, value in (
            filtered_locations.start.apply(lambda x: int(x / chunk_size))
            .value_counts()
            .items()
        ):
            value = value if value < 10 else 10
            data.append(
                {"counts": value, "x": int(index % width), "y": int(index / width)}
            )

        self.__positions = pd.DataFrame(data)

        self.__pixel = pixel
        self.__xoffset = xoffset
        self.__yoffset = yoffset
        self.__gwas_color = gwas_color

    def generate_gwas(self: gwas_annotator) -> str:
        """Generate SVG for GWAS hits.

        Returns:
            str: SVG string containing GWAS hit circles.
        """
        pixel = self.__pixel
        positions = self.__positions
        xoffset = self.__xoffset
        yoffset = self.__yoffset
        gwas_color = self.__gwas_color

        # GWAS Points:
        gwas_points = []

        # Based on the x/y coordinates, let's draw the point:
        for _, row in positions.iterrows():
            # The radius of the circle is proportional to the number of GWAS hits in the given chunk:
            radius = math.sqrt(row["counts"] ** 2 * self.circle_unit / math.pi)

            # Adding point{}
            gwas_points.append(
                self.gwas_hit.format(
                    (row["x"] * pixel) + radius / 2 + xoffset,
                    (row["y"] * pixel) + radius / 2 + yoffset,
                    radius,
                    gwas_color,
                    gwas_color,
                )
            )

        return "\n".join(gwas_points)
