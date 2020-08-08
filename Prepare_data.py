import re
import os
import json
import argparse
import requests
import logging 
from datetime import datetime
from dateutil import parser
import pandas as pd
import ftplib
import sys

# Importing custom functions:
from functions.input_parsers import Fetch_gwas, Fetch_gencode, get_cytoband_data


# Save parameters
def save_parameters(d,config_file):
    return 1


# get Ensembl data:
def get_ensembl_data():
    return 1


# get ensembl version
def get_ensembl_version():
    # http://rest.ensembl.org/info/data/?content-type=application/json
    return 1


# Main 
def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='This script fetches and parses input data for the genome plotter project')
    parser.add_argument('-d', '--dataDir', help='Folder into which the input data and the temporary files will be saved', required=True, type=str)
    parser.add_argument('-c', '--config', help='JSON file with configuration data', required=True, type=str)    
    parser.add_argument('-l', '--logfile', help='Name of the logfile', required=False, type=str)

    # Parse parameters:
    args = parser.parse_args()
    data_dir = args.dataDir
    data_dir = os.path.abspath(data_dir)
    config_file = args.config

    # Initialize logger:
    handlers = [ logging.StreamHandler(sys.stdout) ]
    if args.logfile != '':
        handlers.append( logging.FileHandler(filename=args.logfile) )

    # Initialize logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

    logging.info('Fetching data for Genome Plotter started....')

    # Reading configuration file:
    with open(config_file) as f:
        configuration = json.load(f)

    # Update data folder:
    configuration['basic_parameters']['data_folder'] = data_dir

    # Fetching GWAS Catalog data:
    # gwas_retrieve = Fetch_gwas(configuration['source_data']['gwas_data'])
    # gwas_retrieve.retrieve_data()
    # gwas_retrieve.process_gwas_data()
    # gwas_retrieve.save_gwas_data(data_dir)
    # configuration['source_data']['gwas_data']['release_date'] = gwas_retrieve.get_release_date()

    # Fetching cytological bands:
    cytoband_url = configuration['source_data']['cytoband_data']['url']
    cytoband_output_file = f"{data_dir}/{configuration['source_data']['cytoband_data']['processed_file']}"
    logging.info('Fetching cytoband information.')
    get_cytoband_data(cytoband_url, cytoband_output_file)

    # Fetching GENCODE data:
    logging.info('Fetching GENCODE data.')
    gencode_retrieve = Fetch_gencode(configuration['source_data']['gencode_data'])
    gencode_retrieve.retrieve_data()
    gencode_retrieve.process_gencode_data()
    gencode_retrieve.save_gencode_data(data_dir)
    configuration['source_data']['gencode_data']['release_date'] = gencode_retrieve.get_release_date()
    configuration['source_data']['gencode_data']['version'] = gencode_retrieve.get_release()
    logging.info(f"Saving processed data: {configuration['source_data']['gencode_data']['processed_file']}.")




    # Save config file:
    logging.info('Saving updated configuration.')
    with open(config_file.replace('.json','_updated.json'), 'w') as f:
        json.dump(configuration, f, indent=4)



if __name__ == '__main__':
    main()

