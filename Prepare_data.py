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
from functions.input_parsers import Fetch_gwas
from functions.input_parsers import get_cytoband_data


# Save parameters
def save_parameters(d,config_file):
    return 1


# get Ensembl data:
def get_ensembl_data():
    return 1


# get GENCODE data:
def get_gencode_data():
    return 1


# get ensembl version
def get_ensembl_version():
    # http://rest.ensembl.org/info/data/?content-type=application/json
    return 1


# parse genome:
def parse_genome():
    return 1

       
def get_gwas_bed():
    
    # Open ftp connection and fetch GWAS data:
    logging.info('Fetching GWAS data from ftp...')
    gwas_obj = fetch_gwas_data(gwas_ftp_host)
    release_date = gwas_obj.get_release_date()
    gwas_df = gwas_obj.fetch_associations()
    logging.info('Done. Release date: {}, number of associations: {}'.format(release_date, len(gwas_df)))
    gwas_obj.close_connection()

    # Set proper types:
    logging.info('Processing GWAS data...')
    gwas_df.CHR_ID = gwas_df.CHR_ID.astype(str)

    # Filtering columns:
    filt = gwas_df.loc[(~gwas_df.CHR_ID.isna()) & 
                  (~gwas_df.CHR_POS.isna()),['CHR_ID', 'CHR_POS', 'SNPS']]

    filt = filt.loc[(~filt.CHR_ID.str.contains('x')) &
                  (~filt.CHR_ID.str.contains(';'))&
                   (filt.SNPS.str.contains('rs', case=False))]

    # Set proper types again:
    filt.CHR_POS = filt.CHR_POS.astype(int)

    # rename columns:
    filt['end'] = filt.CHR_POS
    filt.rename(columns={
        'CHR_POS':'start',
        'CHR_ID': '#chr',
        'SNPS':'rsID'
    }, inplace=True)

    # Order columns:
    filt = filt[['#chr','start','end','rsID']]
    
    # Save dataframe:
    logging.info('Saving file ({})'.format(gwas_output_filename))
    filt.to_csv(gwas_output_filename, sep='\t', compression='infer')
    
    # Return release date:
    return release_date

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
    gwas_retrieve = Fetch_gwas(configuration['source_data']['gwas_data'])
    gwas_retrieve.retrieve_data()
    gwas_retrieve.process_gwas_data()
    gwas_retrieve.save_gwas_data(data_dir)

    # Fetching cytological bands:
    cytoband_url = configuration['source_data']['cytoband_data']['url']
    cytoband_output_file = configuration['source_data']['cytoband_data']['processed_file']
    logging.info('Fetching cytoband information.')
    get_cytoband_data(cytoband_url, cytoband_output_file)


    # Save config file:
    logging.info('Saving updated configuration.')
    with open(config_file, 'w') as f:
        json.dump(configuration, f, indent=4)



if __name__ == '__main__':
    main()

