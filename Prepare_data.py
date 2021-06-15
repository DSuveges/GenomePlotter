import os
import json
import argparse
import logging
import sys

# Importing custom functions:
from functions.input_parsers import (
    Fetch_gwas,
    Fetch_genome,
    Fetch_gencode,
    Fetch_cytobands,
    Fetch_ensembl_version
)

def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='This script fetches and parses input data for the genome plotter project'
    )
    parser.add_argument(
        '-d', '--dataDir', help='Folder into which the input data and the temporary files will be saved',
        required=True, type=str
    )
    parser.add_argument(
        '-c', '--config', help='JSON file with configuration data',
        required=True, type=str
    )
    parser.add_argument(
        '-l', '--logfile', help='Name of the logfile',
        required=False, type=str
    )
    parser.add_argument(
        '-s', '--chunk_size', help='Chunk size to pool genomic sequence',
        required=True, type=int
    )
    parser.add_argument(
        '-t', '--tolerance', help='Fraction of a chunk that cannot be N.',
        required=True, type=float
    )

    # Parse parameters:
    args = parser.parse_args()
    data_dir = args.dataDir
    data_dir = os.path.abspath(data_dir)
    config_file = args.config
    chunk_size = args.chunk_size
    tolerance = args.tolerance

    # Initialize logger:
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.logfile != '':
        handlers.append(logging.FileHandler(filename=args.logfile))

    # Initialize logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

    logging.info('Fetching data for Genome Plotter started....')

    # Checking if output dir exists:
    if not os.path.isdir(data_dir):
        raise ValueError(f'The provided folder ({data_dir}) does not exist')

    logging.info(f'Pre-processed data is saved to {data_dir}')

    # Reading configuration file:
    if not os.path.isfile(config_file):
        raise ValueError(f'The provided config file ({config_file}) does not exist.')

    with open(config_file) as f:
        try:
            configuration = json.load(f)
        except json.decoder.JSONDecodeError:
            raise ValueError(f'The provided config file ({config_file}) is not a valid JSON file.')

    logging.info(f'Configuration is read from {config_file}')

    # Update data folder:
    configuration['basic_parameters']['data_folder'] = data_dir
    configuration['basic_parameters']['chunk_size'] = chunk_size
    configuration['basic_parameters']['missing_tolerance'] = tolerance

    # Fetching GWAS Catalog data:
    logging.info('Fetching GWAS data...')
    gwas_retrieve = Fetch_gwas(configuration['source_data']['gwas_data'])
    gwas_retrieve.retrieve_data()
    gwas_retrieve.process_gwas_data()
    gwas_retrieve.save_gwas_data(data_dir)
    configuration['source_data']['gwas_data']['release_date'] = gwas_retrieve.get_release_date()

    # Fetching cytological bands:
    cytoband_url = configuration['source_data']['cytoband_data']['url']
    cytoband_output_file = f"{data_dir}/{configuration['source_data']['cytoband_data']['processed_file']}"
    logging.info('Fetching cytoband information.')
    cytoband_retrieve = Fetch_cytobands(cytoband_url)
    cytoband_retrieve.save_cytoband_data(cytoband_output_file)
    configuration['source_data']['cytoband_data']['genome_build'] = cytoband_retrieve.get_assembly_build()

    # Fetching GENCODE data:
    logging.info('Fetching GENCODE data.')
    gencode_retrieve = Fetch_gencode(configuration['source_data']['gencode_data'])
    gencode_retrieve.retrieve_data()
    gencode_retrieve.process_gencode_data()
    gencode_retrieve.save_gencode_data(data_dir)
    configuration['source_data']['gencode_data']['release_date'] = gencode_retrieve.get_release_date()
    configuration['source_data']['gencode_data']['version'] = gencode_retrieve.get_release()
    logging.info(f"Saving processed data: {configuration['source_data']['gencode_data']['processed_file']}.")

    # Fetching Ensembl version and genome build:
    logging.info('Fetching Ensembl release...')
    ensembl_release_url = configuration['source_data']['ensembl_data']['version_url']
    ensembl_release = Fetch_ensembl_version(ensembl_release_url)
    configuration['source_data']['ensembl_data']['release'] = ensembl_release
    logging.info(f'Current Ensembl release: {ensembl_release}')

    # Fetching the human genome:
    logging.info('Fetching the human genome sequence...')
    genome_retrieve = Fetch_genome(configuration['source_data']['ensembl_data'])
    genome_retrieve.retrieve_data()
    genome_retrieve.parse_genome(chunk_size, tolerance, data_dir)

    # Save config file:
    logging.info(f"Saving updated configuration as {config_file.replace('json', 'updated.json')}")
    with open(config_file.replace('json', 'updated.json'), 'w') as f:
        json.dump(configuration, f, indent=4)


if __name__ == '__main__':
    main()
