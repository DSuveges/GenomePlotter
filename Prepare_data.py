"""This script fetches and prepares the input data for the genome plotter project."""
from __future__ import annotations

import argparse
import json
import logging.config
import os
from dataclasses import asdict

import yaml

from functions.ConfigManager import Config
from functions.InputParsers import (
    FetchCytobands,
    FetchGenome,
    FetchGwas,
    fetch_ensembl_version,
)
from input_parsers.fetch_gencode import FetchGencode


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="This script fetches and parses input data for the genome plotter project"
    )
    parser.add_argument(
        "-d",
        "--dataDir",
        help="Folder into which the input data and the temporary files will be saved",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="JSON file with configuration data",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-l", "--logfile", help="Name of the logfile", required=False, type=str
    )
    parser.add_argument(
        "-s",
        "--chunk_size",
        help="Chunk size to pool genomic sequence",
        required=True,
        type=int,
    )
    parser.add_argument(
        "-t",
        "--tolerance",
        help="Fraction of a chunk that cannot be N.",
        required=True,
        type=float,
    )

    return parser.parse_args()


def get_cytoband_data(cytoband_url: str, cytoband_output_file: str) -> str:
    """Fetch cytoband data from the given URL and save it to the output file.
    
    Args:
        cytoband_url (str): URL to fetch the cytoband data from.
        cytoband_output_file (str): File to save the cytoband data to.
    
    Returns:
        str: The genome build of the cytoband data.
    """
    cytoband_retrieve = FetchCytobands(cytoband_url)
    cytoband_retrieve.save_cytoband_data(cytoband_output_file)
    return cytoband_retrieve.get_assembly_build()


def main(configuration: Config) -> None:
    """Main function to fetch and prepare the input data for the genome plotter project.
    
    Args:
        configuration (Config): The configuration object containing the input data.
    """
    # Extracting relevant parameters:
    basic_parameters = configuration.basic_parameters
    data_dir = basic_parameters.data_folder
    chunk_size = basic_parameters.chunk_size
    tolerance = basic_parameters.missing_tolerance

    # Report the other command line parameters:
    logger.info(f"Chunk size: {chunk_size}")
    logger.info(f"Tolerance for unsequenced bases: {tolerance}")

    # Fetching GWAS Catalog data:
    logger.info("Fetching GWAS data...")
    gwas_retrieve = FetchGwas(asdict(configuration.source_data.gwas_data))
    gwas_retrieve.retrieve_data()
    gwas_retrieve.process_gwas_data()
    gwas_retrieve.save_gwas_data(data_dir)
    configuration.source_data.gwas_data.release_date = gwas_retrieve.get_release_date()

    # Fetching cytological bands:
    logger.info("Fetching cytoband information...")
    configuration.source_data.cytoband_data.genome_build = get_cytoband_data(
        configuration.source_data.cytoband_data.url,
        f"{data_dir}/{configuration.source_data.cytoband_data.processed_file}",
    )

    # Fetching GENCODE data:
    logging.info("Fetching GENCODE data.")
    gencode_retrieve = FetchGencode(asdict(configuration.source_data.gencode_data))
    gencode_retrieve.retrieve_data()
    gencode_retrieve.process_gencode_data()
    gencode_retrieve.save_gencode_data(data_dir)
    configuration.source_data.gencode_data.release_date = (
        gencode_retrieve.get_release_date()
    )
    configuration.source_data.gencode_data.version = gencode_retrieve.get_release()

    # Fetching Ensembl version and genome build:
    logger.info("Fetching Ensembl release...")
    ensembl_release = fetch_ensembl_version(
        configuration.source_data.ensembl_data.version_url
    )
    configuration.source_data.ensembl_data.release = ensembl_release
    logger.info(f"Current Ensembl release: {ensembl_release}")

    # Fetching the human genome:
    logger.info("Fetching the human genome sequence...")
    genome_retrieve = FetchGenome(asdict(configuration.source_data.ensembl_data))
    genome_retrieve.retrieve_data()
    genome_retrieve.parse_genome(chunk_size, tolerance, data_dir)

    # Save config file:
    updated_config_file = "config_updated.json"
    logger.info(f"Saving updated configuration as {updated_config_file}.")
    configuration.save(updated_config_file)


def validate_input(data_dir: str, config_file: str) -> None:
    """Validate the input parameters.
    
    Args:
        data_dir (str): The directory to save the data to.
        config_file (str): The configuration file.
        
    Raises:
        ValueError: If the data directory or the configuration file do not exist.
    """
    # Checking if output dir exists:
    if not os.path.abspath(data_dir):
        raise ValueError(f"The provided folder ({data_dir}) does not exist")

    # Reading configuration file:
    if not os.path.isfile(config_file):
        raise ValueError(f"The provided config file ({config_file}) does not exist.")


if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()

    # Initialise logger:
    with open("logger_config.yaml", "r") as stream:
        logger_config = yaml.safe_load(stream, Loader=yaml.FullLoader)

    logging.config.dictConfig(logger_config)
    logger = logging.getLogger(__name__)

    # Validate input parameters:
    logger.info("Validating input parameters...")
    validate_input(args.dataDir, args.config)

    logger.info(f"Pre-processed data is saved to {args.dataDir}")
    logger.info(f"Configuration file: {args.config}")

    # Initilise configuration:
    with open(args.config) as f:
        try:
            configuration = Config(**json.load(f))
        except json.decoder.JSONDecodeError:
            raise ValueError(
                f"The provided config file ({args.config}) is not a valid JSON file."
            )

    # Update configuration with command line options:
    configuration.update_basic_parameters(
        data_folder=os.path.abspath(args.dataDir),
        chunk_size=args.chunk_size,
        missing_tolerance=args.tolerance,
    )

    main(configuration)
