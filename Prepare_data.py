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

# A function to retrieve data from the web:
def fetch_json(URL):
    """Fetching JSON Data from web

    Args:
        URL (str): URL that will be fetched

    Returns:
        dict: the loaded JSON object
    """

    response = requests.get(URL)

    try:
        parsed_json = response.json
    except:
        raise('[Error] JSON from URL (%s) could not be fetched.'.format(URL))

    return parsed_json


# Save parameters
def save_parameters(d,config_file):
    return 1


# get Ensembl data:
def get_ensembl_data():
    return 1


# get GENCODE data:
def get_gencode_data():
    return 1


# get cytoband data:
def get_cytoband_data(url, outfile):
    response = requests.get(url)
    data = response.json()

    logging.info('Cytobands successfully fetched. Parsing.')

    bands = []
    for region in data['top_level_region']:
        if 'bands' in region:
            bands += region['bands']
            
            
    df = pd.DataFrame(bands)
    df.rename(columns={
        'id': 'name',
        'seq_region_name': 'chr',
        'stain': 'type'
    }, inplace=True)

    logging.info(f'Number of bands in the genome: {len(df)}.')
    logging.info(f'Saving cytoband file: {outfile}.')

    df = df[['chr', 'start', 'end', 'name', 'type']]
    df.sort_values(by = ['chr', 'start'], inplace=True)
    df.to_csv(outfile, sep='\t', index=False, compression='gzip')
    logging.info(f'Outputfile successfully saved.')

# get ensembl version
def get_ensembl_version():
    # http://rest.ensembl.org/info/data/?content-type=application/json
    return 1


# parse genome:
def parse_genome():
    return 1


class fetch_gwas_data(object):
    """
    This class is to retrieve the association table from the most recent GWAS Catalog release.
    Expects the ftp host address.

    It also returns the release date.
    """
    def __init__(self, url):
        self.FTP_HOST = url
        self.release_folder = '/pub/databases/gwas/releases/latest/'
        
        # Initialize connection and go to folder:
        self.ftp = ftplib.FTP(self.FTP_HOST,'anonymous','')
        self.ftp.cwd(self.release_folder)
        
    def get_release_date(self):
        
        # Get list of files and the date of modification:
        files = []
        self.ftp.dir(files.append)
        
        # Parse creation date of a file:
        release_date_string = ' '.join(files[1].split()[5:8])
        release_date = parser.parse(release_date_string)
        
        # Store date:
        self.release_date = release_date.strftime('%Y-%m-%d')
        return self.release_date
    
    def fetch_associations(self):
        return(pd.read_csv('ftp://{}/{}/gwas-catalog-associations_ontology-annotated.tsv '.format(self.FTP_HOST, self.release_folder),
                          sep='\t', dtype=str))
    
    def close_connection(self):
        self.ftp.close()

        
def get_gwas_bed(gwas_ftp_host, gwas_output_filename):
    
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
    gwas_host = configuration['source_data']['gwas_data']['host']
    gwas_output_filename = f"{data_dir}/{configuration['source_data']['gwas_data']['processed_file']}"
    logging.info(f'Fetching GWAS Catalog data from {gwas_host}. Data will be saved as {gwas_output_filename}')
    configuration['source_data']['gwas_data']['release_date'] = get_gwas_bed(gwas_host, gwas_output_filename)

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

