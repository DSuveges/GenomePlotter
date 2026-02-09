"""A method to fetch cytological band information from Ensembl REST API."""

from __future__ import annotations

import logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)


# get cytoband data:
class FetchCytobands:
    """Function to retrieve cytogenic bands from Ensembl."""

    def __init__(self: FetchCytobands, url: str) -> None:
        """Initialize the object with the URL to fetch the cytobands.

        Args:
            url (str): URL to fetch the cytoband data from.
        """
        response = requests.get(url)
        data = response.json()

        logger.info("Cytobands successfully fetched. Parsing.")

        # Saving assembly:
        self.assembly = data["default_coord_system_version"]
        logger.info(f"Current genome assembly: {self.assembly}")

        bands = []
        for region in data["top_level_region"]:
            if "bands" in region:
                bands += region["bands"]

        df = pd.DataFrame(bands)
        df.rename(
            columns={"id": "name", "seq_region_name": "chr", "stain": "type"},
            inplace=True,
        )

        logger.info(f"Number of bands in the genome: {len(df):,}.")

        df = df[["chr", "start", "end", "name", "type"]]
        self.cytobands = df.sort_values(by=["chr", "start"])

    def save_cytoband_data(self: FetchCytobands, outfile: str) -> None:
        """Save the cytoband data to a file.

        Args:
            outfile (str): The file to save the data to.
        """
        logger.info(f"Saving cytoband file: {outfile}.")
        self.cytobands.to_csv(outfile, sep="\t", index=False, compression="gzip")

    def get_assembly_build(self: FetchCytobands) -> str:
        """Return the assembly build.

        Returns:
            str: The assembly build.
        """
        return self.assembly
