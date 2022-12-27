"""Steps to fetch and process required data for genome plotter."""
from __future__ import annotations

from typing import TYPE_CHECKING

import hydra

if TYPE_CHECKING:
    from omegaconf import DictConfig

import os
import json
import logging
import sys

from modules.InputParsers import (
    FetchGwas,
    FetchGenome,
    FetchCytobands,
    FetchGencode,
    fetch_ensembl_version
)

@hydra.main(version_base=None, config_path=".", config_name="config")
def main(cfg: DictConfig) -> None:
    """Downloading data."""

    # Initialize logger:
    handlers = [logging.StreamHandler(sys.stdout)]

    # Initialize logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

    logging.info('Fetching data for Genome Plotter started....')

    # Checking if output dir exists:
    if not os.path.isdir(cfg.shared_parameters.source_data_folder):
        raise ValueError(f'The provided folder ({cfg.shared_parameters.source_data_folder}) does not exist')

    logging.info(f'Pre-processed data is saved to {cfg.shared_parameters.source_data_folder}')

    # Report the other command line parameters:
    logging.info(f'Chunk size: {cfg.basic_parameters.chunk_size}')
    logging.info(f'Tolerance for unsequenced bases: {tolerance}')

    # Fetching GWAS Catalog data:
    logging.info('Fetching GWAS data...')
    gwas_retrieve = FetchGwas(configuration['source_data']['gwas_data'])
    gwas_retrieve.retrieve_data()
    gwas_retrieve.process_gwas_data()
    gwas_retrieve.save_gwas_data(cfg.shared_parameters.source_data_folder)
    configuration['source_data']['gwas_data']['release_date'] = gwas_retrieve.get_release_date()

    # Fetching cytological bands:
    cytoband_url = configuration['source_data']['cytoband_data']['url']
    cytoband_output_file = f"{cfg.shared_parameters.source_data_folder}/{configuration['source_data']['cytoband_data']['processed_file']}"
    logging.info('Fetching cytoband information.')
    cytoband_retrieve = FetchCytobands(cytoband_url)
    cytoband_retrieve.save_cytoband_data(cytoband_output_file)
    configuration['source_data']['cytoband_data']['genome_build'] = cytoband_retrieve.get_assembly_build()

    # Fetching GENCODE data:
    logging.info('Fetching GENCODE data.')
    gencode_retrieve = FetchGencode(configuration['source_data']['gencode_data'])
    gencode_retrieve.retrieve_data()
    gencode_retrieve.process_gencode_data()
    gencode_retrieve.save_gencode_data(data_dir)
    configuration['source_data']['gencode_data']['release_date'] = gencode_retrieve.get_release_date()
    configuration['source_data']['gencode_data']['version'] = gencode_retrieve.get_release()
    logging.info(f"Saving processed data: {configuration['source_data']['gencode_data']['processed_file']}.")

    # Fetching Ensembl version and genome build:
    logging.info('Fetching Ensembl release...')
    ensembl_release_url = configuration['source_data']['ensembl_data']['version_url']
    ensembl_release = fetch_ensembl_version(ensembl_release_url)
    configuration['source_data']['ensembl_data']['release'] = ensembl_release
    logging.info(f'Current Ensembl release: {ensembl_release}')

    # Fetching the human genome:
    logging.info('Fetching the human genome sequence...')
    genome_retrieve = FetchGenome(configuration['source_data']['ensembl_data'])
    genome_retrieve.retrieve_data()
    genome_retrieve.parse_genome(chunk_size, tolerance, data_dir)

    # Save config file:
    logging.info(f"Saving updated configuration as {config_file.replace('json', 'updated.json')}")
    with open(config_file.replace('json', 'updated.json'), 'w') as f:
        json.dump(configuration, f, indent=4)


if __name__ == '__main__':
    main()
