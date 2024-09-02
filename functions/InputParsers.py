from __future__ import annotations

import logging
import re

import pandas as pd
import requests

# Load custom packages:
from .FetchFromFtp import FetchFromFtp

logger = logging.getLogger(__name__)

# get ensembl version
def fetch_ensembl_version(url):
    response = requests.get(url)
    data = response.json()
    return data['releases'][0]


# Fetch GWAS Catalog data and parse:
class FetchGwas(FetchFromFtp):

    """Function to fetch GWAS data from ftp and format as bed"""

    def __init__(self, gwas_parameters):
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

        logger.info(f'Successfully fetched {len(self.tsv_data)} GWAS associations.')

        # Close connection:
        self.close_connection()

    # Processing gwas data:
    def process_gwas_data(self):

        gwas_df = self.tsv_data
        gwas_df.CHR_ID = gwas_df.CHR_ID.astype(str)

        # Filtering columns:
        filt = gwas_df.loc[
            (~gwas_df.CHR_ID.isna())
            & (~gwas_df.CHR_POS.isna()), ['CHR_ID', 'CHR_POS', 'SNPS']
        ]

        filt = filt.loc[
            (~filt.CHR_ID.str.contains('x'))
            & (~filt.CHR_ID.str.contains(';'))
            & (filt.SNPS.str.contains('rs', case=False))
        ]

        logger.info(f'Number of filtered associations:  {len(filt)}.')
        logger.info('Formatting data...')

        # Set proper types again:
        filt.CHR_POS = filt.CHR_POS.astype(int)

        # rename columns:
        filt['end'] = filt.CHR_POS
        filt.rename(
            columns={
                'CHR_POS': 'start',
                'CHR_ID': '#chr',
                'SNPS': 'rsID'
            }, inplace=True
        )

        # Order columns and getting rid of duplicates:
        filt = filt[['#chr', 'start', 'end', 'rsID']].drop_duplicates()

        # Save dataframe:
        self.gwas_df = filt

    # Save data
    def save_gwas_data(self, data_dir):
        gwas_output_filename = f'{data_dir}/{self.processed_file}'
        self.gwas_df.to_csv(gwas_output_filename, sep='\t', compression='infer', index=False)

    # Extract release date:
    def get_release_date(self):
        return self.release_date


# get cytoband data:
class FetchCytobands(object):

    """Function to retrieve cytogenic bands from Ensembl"""

    def __init__(self, url):

        response = requests.get(url)
        data = response.json()

        logger.info('Cytobands successfully fetched. Parsing.')

        # Saving assembly:
        self.assembly = data['default_coord_system_version']
        logger.info(f"Current genome assembly: {self.assembly}")

        bands = []
        for region in data['top_level_region']:
            if 'bands' in region:
                bands += region['bands']

        df = pd.DataFrame(bands)
        df.rename(
            columns={
                'id': 'name',
                'seq_region_name': 'chr',
                'stain': 'type'
            }, inplace=True
        )

        logger.info(f'Number of bands in the genome: {len(df)}.')

        df = df[['chr', 'start', 'end', 'name', 'type']]
        self.cytobands = df.sort_values(by=['chr', 'start'])

    def save_cytoband_data(self, outfile):
        logger.info(f'Saving cytoband file: {outfile}.')
        self.cytobands.to_csv(outfile, sep='\t', index=False, compression='gzip')

    def get_assembly_build(self):
        return self.assembly


# Fetch and process Gencode data:
class FetchGenome(FetchFromFtp):

    """Function to retrieve genome sequence from Ensembl"""

    def __init__(self, ensembl_parameters):

        # These are the parameters required to fetch the data:
        self.host = ensembl_parameters['host']
        self.path = ensembl_parameters['path'].format(ensembl_parameters['release'])
        self.source_file = ensembl_parameters['source_file']
        self.parsed_file = ensembl_parameters['processed_file']

        # Initialize host:
        FetchFromFtp.__init__(self, self.host)

    def retrieve_data(self):

        ftp = FetchFromFtp(self.host)
        self.resp = ftp.fetch_file(self.path, self.source_file)

        logger.info('Sequence data successfully fetched. Parsing...')

    def parse_genome(self, chunk_size, threshold, data_folder):

        self.chunk_size = chunk_size
        self.threshold = threshold
        self.data_folder = data_folder

        # This variable contains sequence for one chromosome:
        chrom_data = ''
        chrom_name = None

        for line in self.resp:
            line = line.decode("utf-8")

            # Check if it is a header:
            if re.match('>', line):

                # If there's data in the buffer, save it:
                if chrom_data:
                    if len(chrom_name) < 3:
                        logger.info(f'Parsing chromosome {chrom_name} is done.')
                        self.process_chromosome(chrom_data, chrom_name)
                    else:
                        logger.info(f'Chromosome {chrom_name} is skipped.')

                    chrom_data = ''

                x = re.match(r'>(\S+) ', line)
                chrom_name = x.group(1)

            chrom_data += line.strip()

        # The last chunk is passed:
        logger.info(f'Parsing chromosome {chrom_name} is done.')
        self.process_chromosome(chrom_data, chrom_name)

    def process_chromosome(self, chrom_data, chr_name):

        raw_data = []
        file_name = f"{self.data_folder}/{self.parsed_file.format(chr_name)}"
        chunk_size = self.chunk_size
        threshold = self.threshold

        for i in range(0, len(chrom_data), chunk_size):
            # Removing Ns:
            chunk = chrom_data[i: i + chunk_size].replace('N', '')

            # Skipping chunks where the Ns are above the threshold:
            if len(chunk) < chunk_size * threshold:
                gc_content = None
            else:
                # Calculate GC content:
                gc_content = (chunk.count('G') + chunk.count('C')) / len(chunk)

            raw_data.append({
                'chr': chr_name,
                'start': i,
                'end': i + chunk_size,
                'GC_ratio': gc_content
            })

        # Save data:
        df = pd.DataFrame(raw_data)
        df.to_csv(file_name, sep='\t', compression='infer', index=False, na_rep='NA')
