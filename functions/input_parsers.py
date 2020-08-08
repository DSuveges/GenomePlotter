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
import logging


# Load custom packages:
from .Fetch_from_ftp import Fetch_from_ftp


class Fetch_gwas(Fetch_from_ftp):
    def __init__(self,gwas_parameters):
        """
        Based on the gwas source related parameters, fetches, parses and saves data
        """
        
        # Extract parameters:
        self.host = gwas_parameters['host']
        self.path = gwas_parameters['path']
        self.source_file = gwas_parameters['source_file']
        self.processed_file = gwas_parameters['processed_file']
        self.release_date = None
        
        # Initialize host:
        Fetch_from_ftp.__init__(self, self.host)
        
        
    def retrieve_data(self):
        """
        Retrieving release data and association table.
        Then the connection is closed.
        """
        
        # Get release date
        self.release_date = self.fetch_last_update_date(self.path)
    
        # Parse data
        self.fetch_tsv(self.path, self.source_file)
        
        # Close connection:
        self.close_connection()
        
        
    # Processing gwas data:
    def process_gwas_data(self):
        
        gwas_df = self.tsv_data
        gwas_df.CHR_ID = gwas_df.CHR_ID.astype(str)

        # Filtering columns:
        filt = gwas_df.loc[(~gwas_df.CHR_ID.isna()) & 
              (~gwas_df.CHR_POS.isna()),['CHR_ID', 'CHR_POS', 'SNPS']]

        filt = filt.loc[
              (~filt.CHR_ID.str.contains('x')) &
              (~filt.CHR_ID.str.contains(';')) &
              (filt.SNPS.str.contains('rs', case=False))
        ]

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
        self.gwas_df = filt

        
    # Save data
    def save_gwas_data(self,data_dir):
        gwas_output_filename = f'{data_dir}/{self.processed_file}'
        self.gwas_df.to_csv(gwas_output_filename, sep='\t', compression='infer', index=False)
        
    
    # Extract release date:
    def get_release_date(self):
        return self.release_date


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
