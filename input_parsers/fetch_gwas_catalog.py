"""Module to fetch and parse the human genome sequence data from Ensembl FTP."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from functions.FetchFromFtp import FetchFromFtp

if TYPE_CHECKING:
    from functions.ConfigManager import SourcePrototype

logger = logging.getLogger(__name__)


# Fetch GWAS Catalog data and parse:
class FetchGwas(FetchFromFtp):
    """Function to fetch GWAS data from ftp and format as bed"""

    def __init__(self, gwas_parameters: SourcePrototype):
        """
        Based on the gwas source related parameters, fetches, parses and saves data
        """

        # Extract parameters:
        self.host = gwas_parameters.host
        self.path = gwas_parameters.path
        self.source_file = gwas_parameters.source_file
        self.processed_file = gwas_parameters.processed_file
        self.release_date = None

        # Initialize host:
        FetchFromFtp.__init__(self, self.host)

    def retrieve_data(self):
        """
        Retrieving release data and association table.
        Then the connection is closed.
        """

        # Get release date
        self.release_date = self.fetch_last_update_date(self.path)

        # Parse data
        self.fetch_tsv(self.path, self.source_file)

        logger.info(f"Successfully fetched {len(self.tsv_data):,} GWAS associations.")

        # Close connection:
        self.close_connection()

    # Processing gwas data:
    def process_gwas_data(self):
        gwas_df = self.tsv_data
        gwas_df.CHR_ID = gwas_df.CHR_ID.astype(str)

        # Filtering columns:
        filt = gwas_df.loc[
            (~gwas_df.CHR_ID.isna()) & (~gwas_df.CHR_POS.isna()),
            ["CHR_ID", "CHR_POS", "SNPS"],
        ]

        filt = filt.loc[
            (~filt.CHR_ID.str.contains("x"))
            & (~filt.CHR_ID.str.contains(";"))
            & (filt.SNPS.str.contains("rs", case=False))
        ]

        logger.info(f"Number of filtered associations:  {len(filt):,}.")
        logger.info("Formatting data...")

        # Set proper types again:
        filt.CHR_POS = filt.CHR_POS.astype(int)

        # rename columns:
        filt["end"] = filt.CHR_POS
        filt.rename(
            columns={"CHR_POS": "start", "CHR_ID": "#chr", "SNPS": "rsID"}, inplace=True
        )

        # Order columns and getting rid of duplicates:
        filt = filt[["#chr", "start", "end", "rsID"]].drop_duplicates()

        # Save dataframe:
        self.gwas_df = filt

    # Save data
    def save_gwas_data(self, data_dir):
        gwas_output_filename = f"{data_dir}/{self.processed_file}"
        self.gwas_df.to_csv(
            gwas_output_filename, sep="\t", compression="infer", index=False
        )

    # Extract release date:
    def get_release_date(self):
        return self.release_date
