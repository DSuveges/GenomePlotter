"""This script fetches and prepares the input data for the genome plotter project."""

from __future__ import annotations

import argparse
import os

import yaml
from loguru import logger

from genome_plotter import LOG_FORMAT
from genome_plotter.functions.ConfigManager import Config
from genome_plotter.input_parsers.data_integrator import integrate_data
from genome_plotter.input_parsers.fetch_cytobands import FetchCytobands
from genome_plotter.input_parsers.fetch_ensembl import (
    FetchGenome,
    fetch_ensembl_version,
)
from genome_plotter.input_parsers.fetch_gencode import FetchGencode
from genome_plotter.input_parsers.fetch_gwas_catalog import FetchGwas


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
        help="JSON file with configuration data. If not provided, the default config bundled with the package is used.",
        required=False,
        default=None,
        type=str,
    )
    parser.add_argument(
        "-s",
        "--chunk_size",
        help="Chunk size to pool genomic sequence in base pairs.",
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


def _validate_run_params(
    data_dir: str | None,
    chunk_size: int | None,
    tolerance: float | None,
) -> tuple[str, int, float]:
    """Validate parameters required by run() and return them as non-optional types.

    Raises:
        ValueError: If any required parameter is missing or invalid.
    """
    if not data_dir:
        raise ValueError("Data directory is not provided.")
    if not chunk_size:
        raise ValueError("Chunk size is not provided.")
    if tolerance is None:
        raise ValueError("Tolerance for unsequenced bases is not provided.")
    return data_dir, chunk_size, tolerance


def _fetch_gwas(configuration: Config, data_dir: str) -> None:
    logger.info("Fetching GWAS data...")
    gwas_retrieve = FetchGwas(configuration.source_data.gwas_data)
    gwas_retrieve.retrieve_data()
    gwas_retrieve.process_gwas_data()
    gwas_retrieve.save_gwas_data(data_dir)
    configuration.source_data.gwas_data.release_date = gwas_retrieve.get_release_date()


def _fetch_cytobands(configuration: Config, data_dir: str) -> None:
    logger.info("Fetching cytoband information...")
    configuration.source_data.cytoband_data.genome_build = get_cytoband_data(
        configuration.source_data.cytoband_data.url,
        f"{data_dir}/{configuration.source_data.cytoband_data.processed_file}",
    )


def _fetch_gencode(configuration: Config, data_dir: str) -> None:
    logger.info("Fetching GENCODE data.")
    gencode_retrieve = FetchGencode(configuration.source_data.gencode_data)
    gencode_retrieve.retrieve_data()
    gencode_retrieve.process_gencode_data()
    gencode_retrieve.save_gencode_data(data_dir)
    configuration.source_data.gencode_data.release_date = gencode_retrieve.get_release_date()
    configuration.source_data.gencode_data.version = gencode_retrieve.get_release()


def _fetch_genome(
    configuration: Config, data_dir: str, chunk_size: int, tolerance: float
) -> None:
    logger.info("Fetching Ensembl release...")
    if configuration.source_data.ensembl_data.version_url is None:
        raise ValueError("Could not pull Ensembl version.")
    ensembl_release = fetch_ensembl_version(
        configuration.source_data.ensembl_data.version_url
    )
    configuration.source_data.ensembl_data.release = ensembl_release
    logger.info(f"Current Ensembl release: {ensembl_release}")

    logger.info("Fetching the human genome sequence...")
    genome_retrieve = FetchGenome(configuration.source_data.ensembl_data)
    genome_retrieve.retrieve_data()
    genome_retrieve.parse_genome(chunk_size, tolerance, data_dir)


def run(configuration: Config) -> None:
    """Run the data preparation pipeline.

    Args:
        configuration (Config): The configuration object containing the input data.

    Raises:
        ValueError: If required parameters are missing or invalid.
    """
    basic_parameters = configuration.basic_parameters
    data_dir, chunk_size, tolerance = _validate_run_params(
        basic_parameters.data_folder,
        basic_parameters.chunk_size,
        basic_parameters.missing_tolerance,
    )

    logger.info(f"Chunk size: {chunk_size}")
    logger.info(f"Tolerance for unsequenced bases: {tolerance}")

    _fetch_gwas(configuration, data_dir)
    _fetch_cytobands(configuration, data_dir)
    _fetch_gencode(configuration, data_dir)
    _fetch_genome(configuration, data_dir, chunk_size, tolerance)

    logger.info("Integrating parsed data...")
    integrate_data(
        output_dir=data_dir,
        chromosomes=[str(i + 1) for i in range(22)] + ["Y", "X", "MT"],
        cytoband_file=configuration.get_cytoband_file(),
        gencode_file=configuration.get_gencode_file(),
    )

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
    if not os.path.exists(os.path.join(os.getcwd(), data_dir)):
        raise ValueError(f"The provided folder ({data_dir}) does not exist")

    if not os.path.isfile(config_file):
        raise ValueError(f"The provided config file ({config_file}) does not exist.")


def main() -> None:
    """Entry point for the prepare-data CLI command."""
    args = parse_args()

    logger.add("genome_plotter.log", level="DEBUG", format=LOG_FORMAT)

    if args.config is None:
        args.config = os.path.join(os.path.dirname(__file__), "assets", "config.yaml")
        logger.info("No config file provided, using default bundled config.")

    logger.info("Validating input parameters...")
    validate_input(args.dataDir, args.config)

    logger.info(f"Pre-processed data is saved to {args.dataDir}")
    logger.info(f"Configuration file: {args.config}")

    with open(args.config) as f:
        try:
            configuration = Config.model_validate(yaml.safe_load(f))
        except yaml.YAMLError as e:
            raise ValueError(
                f"The provided config file ({args.config}) is not a valid YAML file."
            ) from e

    configuration.update_basic_parameters(
        data_folder=os.path.abspath(args.dataDir),
        chunk_size=args.chunk_size,
        missing_tolerance=args.tolerance,
    )

    run(configuration)


if __name__ == "__main__":
    main()
