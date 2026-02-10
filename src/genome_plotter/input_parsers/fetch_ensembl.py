"""Module to fetch and parse the human genome sequence data from Ensembl FTP."""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import pandas as pd
import requests

from genome_plotter.functions.FetchFromFtp import FetchFromFtp

if TYPE_CHECKING:
    from genome_plotter.functions.ConfigManager import SourcePrototype

logger = logging.getLogger(__name__)


# get ensembl version
def fetch_ensembl_version(url: str) -> int:
    """Fetch the current Ensembl release version.

    Args:
        url (str): URL to the Ensembl REST API.

    Returns:
        int: Current Ensembl release number.
    """
    response = requests.get(url)
    data = response.json()
    return data["releases"][0]


# Fetch and process Gencode data:
class FetchGenome(FetchFromFtp):
    """Function to retrieve genome sequence from Ensembl."""

    chromosomes: list[str] = []

    def __init__(self: FetchGenome, ensembl_parameters: SourcePrototype) -> None:
        """Initialize the class with the parameters required to fetch the data.

        Args:
            ensembl_parameters (SourcePrototype): Parameters to fetch the data.
        """
        assert ensembl_parameters.release and ensembl_parameters.path, (
            "Ensembl release and path needs to be provided."
        )

        # These are the parameters required to fetch the data:
        self.host = ensembl_parameters.host
        self.path = ensembl_parameters.path.format(ensembl_parameters.release)
        self.source_file = ensembl_parameters.source_file
        self.parsed_file = ensembl_parameters.processed_file

        # Initialize host:
        assert self.host is not None, "Host is required for genome data fetch."
        FetchFromFtp.__init__(self, self.host)

    def retrieve_data(self: FetchGenome) -> None:
        """Retrieve the genome data from the FTP server."""
        assert self.host is not None, "Host is required."
        assert self.source_file is not None, "Source file is required."
        ftp = FetchFromFtp(self.host)
        self.resp = ftp.fetch_file(self.path, self.source_file)

        logger.info("Sequence data successfully fetched. Parsing...")

    def parse_genome(
        self: FetchGenome, chunk_size: int, threshold: float, data_folder: str
    ) -> None:
        """Parse the genome data and save the processed data.

        Args:
            chunk_size (int): Size of the chunk to process.
            threshold (float): Threshold to skip chunks.
            data_folder (str): Folder to save the processed data.
        """
        self.chunk_size = chunk_size
        self.threshold = threshold
        self.data_folder = data_folder

        # This variable contains sequence for one chromosome:
        chrom_data = ""
        chrom_name: str | None = None

        for line_bytes in self.resp:
            line: str = line_bytes.decode("utf-8")

            # Process header:
            if re.match(">", line):
                # If there's data in the buffer, save it:
                if chrom_data:
                    # We are skipping non-canonical chromosomes:
                    self.process_chromosome(chrom_data, chrom_name)

                    # Empty chromosome data buffer:
                    chrom_data = ""

                # Extract chromosome name from header:
                x = re.match(r">(\S+) ", line)
                try:
                    chrom_name = x.group(1) if x else None
                except AttributeError:
                    logger.error(f"Error parsing chromosome name: {line}")
                    raise ValueError(f"Error parsing chromosome name: {line}")

            # Append the sequence:
            chrom_data += line.strip()

        # The last chunk is passed:
        if chrom_name:
            self.chromosomes.append(chrom_name)

        logger.info(f"Parsing chromosome {chrom_name} is done.")
        self.process_chromosome(chrom_data, chrom_name)

    def process_chromosome(
        self: FetchGenome, chrom_data: pd.DataFrame, chr_name: str | None
    ) -> None:
        """Process the chromosome sequence data into defined chunks. Save resulting dataset into tsv.

        Args:
            chrom_data (pd.DataFrame): The chromosome sequence data.
            chr_name (str): The name of the chromosome.
        """
        # Test chromosome name if we want to process or not:
        if chr_name is None or len(chr_name) > 3:
            logger.warning(f"Non-canoncial chromosome ({chr_name}) is skipped.")
            return None

        logger.info(f"Processing chromosome: {chr_name}")
        raw_data = []
        file_name = f"{self.data_folder}/{self.parsed_file.format(chr_name)}"
        chunk_size = self.chunk_size
        threshold = self.threshold

        for i in range(0, len(chrom_data), chunk_size):
            # Removing Ns:
            chunk = chrom_data[i : i + chunk_size].replace("N", "")

            # Skipping chunks where the Ns are above the threshold:
            if len(chunk) < chunk_size * threshold:
                gc_content = None
            else:
                # Calculate GC content:
                gc_content = (chunk.count("G") + chunk.count("C")) / len(chunk)

            raw_data.append(
                {
                    "chr": chr_name,
                    "start": i,
                    "end": i + chunk_size,
                    "GC_ratio": gc_content,
                }
            )

        # Save data:
        df = pd.DataFrame(raw_data)
        df.to_csv(file_name, sep="\t", compression="infer", index=False, na_rep="NA")
