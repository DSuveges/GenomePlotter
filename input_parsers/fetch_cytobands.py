"""A method to fetch cytological band information from Ensembl REST API."""

from __future__ import annotations

import logging

import pandas as pd
import requests

logger = logging.getLogger(__name__)


# get cytoband data:
class FetchCytobands(object):
    """Function to retrieve cytogenic bands from Ensembl"""

    def __init__(self, url):
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

    def save_cytoband_data(self, outfile):
        logger.info(f"Saving cytoband file: {outfile}.")
        self.cytobands.to_csv(outfile, sep="\t", index=False, compression="gzip")

    def get_assembly_build(self):
        return self.assembly
