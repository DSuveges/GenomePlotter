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


# Fetch and process Gencode data:
class Fetch_gencode(Fetch_from_ftp):
    def __init__(self,gencode_parameters):
        """
        Based on the gwas source related parameters, fetches, parses and saves data
        """
        
        # Extract parameters:
        self.host = gencode_parameters['host']
        self.path = gencode_parameters['path']
        self.source_file = gencode_parameters['source_file']
        self.processed_file = gencode_parameters['processed_file']
        self.release_date = None
        self.release = None
        
        # Initialize host:
        Fetch_from_ftp.__init__(self, self.host)
        
        
    def retrieve_data(self):
        """
        Retrieving release data and association table.
        Then the connection is closed.
        """
        
        # Get last release:
        file_list = self.fetch_file_list(self.path)
        releases = [int(x.replace('release_','')) for x in file_list if re.match('release_\d+$', x)]
        self.release = max(releases)

        # Get date of last release: 
        path_to_last_release = '{}/release_{}/'.format(self.path, self.release)
        self.release_date = self.fetch_last_update_date(path_to_last_release)

        # Fetch data:
        self.fetch_tsv(path_to_last_release, self.source_file.format(self.release), skiprows=5, header=None)

        # Parse data:
        columns = {0: 'chr', 2: 'type', 3: 'start', 4: 'end', 8: 'annotation'}
        self.gencode_raw = self.tsv_data.rename(columns=columns)[columns.values()]

        # Close connection:
        self.close_connection()

        logging.info('GENCODE data retrieved.')
        logging.info(f'Number of rows in the table: {len(self.gencode_raw)}')
        logging.info(f'Number of genes in {self.release} release: {len(self.gencode_raw.loc[self.gencode_raw.type == "gene"])}')
        
    # Processing gwas data:
    def process_gencode_data(self):

        # Parsing gtf annotation:
        logging.info('Parsing GTF annotation.')
        parsed_annotation = self.gencode_raw.annotation.apply(
            lambda annotation: { 
                x.strip().split(' ', 1)[0] : x.strip().split(' ', 1)[1].replace('"','') for x in annotation.split(';') if x != ''
            }
        )

        # Merging annotation with coordinates:
        gencode_df_updated = self.gencode_raw.merge(pd.DataFrame(parsed_annotation.tolist()), left_index=True, right_index=True)

        # Drop unparsed annotation column:
        gencode_df_updated.drop(['annotation'], axis=1, inplace=True)

        # Filtering for protein coding genes:
        gencode_df_updated = gencode_df_updated.loc[gencode_df_updated.gene_type == 'protein_coding']
        logging.info(f"Number of protein coding genes: {len(gencode_df_updated.loc[(gencode_df_updated.gene_type == 'protein_coding')&(gencode_df_updated.type == 'gene')])}")

        # Updating types:
        gencode_df_updated['start'] = gencode_df_updated.start.astype(int)
        gencode_df_updated['end'] = gencode_df_updated.end.astype(int)
        
        # Initialize empty dataframe:
        processed = pd.DataFrame(columns=['chr','start','end','type','gene_id','gene_name','transcript_id'])

        logging.info("Generate exon/intron annotations for the canonical transcripts for each gene... (it will take a while.)")
        for gene_id, features in gencode_df_updated.groupby(['gene_id']):

            # Adding length to all features:
            features['length'] = features[['start', 'end']].apply(lambda row: row['end'] - row['start'], axis=1)
            
            # Selecting protein coding transcripts:
            transcripts = features.loc[(features.type == 'transcript') & 
                                       (features.transcript_type == 'protein_coding')]
            
            # If no protein coding transcript is found, we skip gene:
            if len(transcripts) == 0:
                continue

            # Adding CDS length to all transcripts:
            transcripts['CDS_length'] = transcripts.transcript_id.apply(lambda t_id: features.loc[(features.type == 'CDS')&(features.transcript_id == t_id)].length.sum())
            
            # Get canonical transcript and properties:
            canonical_transcript_id = self.get_canonical_transcript(transcripts)
            [start, end] = transcripts.loc[transcripts.transcript_id == canonical_transcript_id, ['start', 'end']].iloc[0].tolist()
            
            # Generate exon-intron splice:
            gene_df = self.generate_exon_intron_structure(gene_id, canonical_transcript_id, start, end, features.loc[(features.transcript_id == canonical_transcript_id) &
                                                                             (features.type == 'exon')])
            # appending to existing data:
            processed = processed.append(gene_df, ignore_index=True)
            

        # Saving data:
        self.processed = processed
        
        
    # Save data
    def save_gencode_data(self,data_dir):
        gwas_output_filename = f'{data_dir}/{self.processed_file}'
        self.processed.to_csv(gwas_output_filename, sep='\t', compression='infer', index=False)
        
    
    # Extract release date:
    def get_release_date(self):
        return self.release_date


    # Extract release date:
    def get_release(self):
        return self.release

    @staticmethod
    def get_canonical_transcript(transcripts):
        ##
        ## Selecting canonical transcript following Ensembl guidelines:
        ## http://www.ensembl.org/Help/Glossary
        ##
        
        # First choice: selecting the longest ccds transcript:
        if transcripts.ccdsid.notna().any():
            longest_CDS = transcripts.loc[~transcripts.ccdsid.isna()].CDS_length.max()
            canonical_transcript_id = transcripts.loc[ longest_CDS == transcripts.CDS_length].transcript_id.tolist()[0]

        # Secong choice: selecting the longest havana transcript:
        elif transcripts.havana_transcript.notna().any():
            longest_CDS = transcripts.loc[~transcripts.havana_transcript.isna()].CDS_length.max()
            canonical_transcript_id = transcripts.loc[ longest_CDS == transcripts.CDS_length].transcript_id.tolist()[0]

        # Third choice: selecting the longest protein coding transcript:
        else:
            longest_CDS = transcripts.CDS_length.max()
            canonical_transcript_id = transcripts.loc[ longest_CDS == transcripts.CDS_length].transcript_id.tolist()[0]
        
        return canonical_transcript_id

    @staticmethod
    def generate_exon_intron_structure(gene_id, transcript_id, start, end, exons):
        
        # Looping through all exons generate coordinates of the corresponding intron(s)
        new_features = []
        for i, values in exons[['start', 'end']].sort_values(by='start', axis=0, ascending=True).iterrows():
            es = values['start']
            ee = values['end']

            # Adding intron:
            if es > start:
                new_features.append({'start': start,'end': es,'type': 'intron'})

            # Adding exon to feature:
            new_features.append({'start': es,'end': ee,'type': 'exon'})

            # updating start position:
            start = ee

        # Adding final exon if exists:
        if start < end:
            new_features.append({'start': start,'end': transript_end,'type': 'intron'})   

        # Generate dataframe + add extra features:
        df = pd.DataFrame(new_features)
        df['gene_id'] = gene_id
        df['transcript_id'] = transcript_id
        df['chr'] = exons.chr.tolist()[0].replace('chr','')
        df['gene_name'] = exons.gene_name.tolist()[0]
        
        return df

