import numpy as np
import pickle
import logging

import pybedtools

from .ColorFunctions import color_darkener, color_picker

class DataIntegrator(object):

    """
    This function assigns color for each chunk in the chromosome
    based on GC content + GENCODE annotation + darkened fraction
    """

    __required_columns = ['chr', 'start', 'end']

    def __init__(self, genome_df):

        self.__genome__ = genome_df.copy()
        self.chromosome_name = genome_df.iloc[0]['chr']

        logging.info(f'Integrating data on chromosome: {self.chromosome_name}')
        logging.info(f'Number of chunks on chromosome {self.chromosome_name}: {len(self.__genome__)}')

        # Testing columns:
        for col in self.__required_columns:
            if col not in genome_df.columns:
                raise ValueError(f'Manadatory colum: {col} is not found in the provided dataframe.')

    def add_xy_coordinates(self, width):
        self.__width__ = width

        self.__genome__ = (
            self.__genome__
            .assign(
                # Get x position of the chunk:
                x=self.__genome__.index.astype(int) % width,
                # Set y position of the chunk:
                y=self.__genome__.index.astype(int) / width
            )
            # Set proper type:
            .astype({'y': 'int32'})
        )
        logging.info(f'Number of chunks in one row: {width}')
        logging.info(f'Number of rows: {self.__genome__.y.max()}')

    def add_genes(self, GENCODE_df):
        logging.info(f'Number of gencode features: {len(GENCODE_df)}')

        # Filtering GENCODE data:
        GENCODE_df = GENCODE_df.loc[GENCODE_df.chr == self.chromosome_name]

        logging.info(f'Number of gencode features on chromosome {self.chromosome_name}: {len(GENCODE_df)}')

        # Creating bedtools objects:
        gencode_bed = pybedtools.bedtool.BedTool.from_dataframe(GENCODE_df.rename(columns={'chr': 'chrom'}))
        chrom_bed = pybedtools.bedtool.BedTool.from_dataframe(self.__genome__.rename(columns={'chr': 'chrom'}))

        # Run intersectbed and extract result as dataframe:
        try:
            GencodeIntersect = chrom_bed.intersect(gencode_bed, wa=True, wb=True)
            intersect_df = GencodeIntersect.to_dataframe(
                header=None,
                names=[
                    'chr', 'start', 'end', 'GC_ratio', 'x', 'y', 'chr_2',
                    'start_2', 'end_2', 'gene_id', 'gene_name',
                    'transcript_id', 'type'
                ]
            )
        except Exception as e:
            print('Gencode data:')
            print(gencode_bed.head())

            print('Chromosome data:')
            print(chrom_bed.head())

            raise e

        # Parse out results:
        GENCODE_chunks = intersect_df.groupby('start').apply(
            lambda x: 'exon' if 'exon' in x.type.unique() else 'gene'
        )
        GENCODE_chunks.name = "GENCODE"

        # Updating index:
        genome_df = self.__genome__.merge(GENCODE_chunks, left_on='start', right_index=True, how='left')
        genome_df.GENCODE.fillna('intergenic', inplace=True)

        # Adding annotation to df:
        self.__genome__ = genome_df

    def get_data(self):
        return(self.__genome__.copy())

    def add_centromere(self, cytoband_df):

        logging.info(self.__genome__.GENCODE.value_counts())
        chromosome = self.__genome__.chr[1]
        centromer_loc = cytoband_df.loc[
            (cytoband_df.chr == str(chromosome))
            & (cytoband_df.type == 'acen'), ['start', 'end']
        ]
        centromer_loc = (int(centromer_loc.start.min()),
                         int(centromer_loc.end.max()))

        # If GENCODE column is missing, let's initialize:
        if 'GENCODE' not in self.__genome__.columns:
            self.__genome__['GENCODE'] = None

        # Assigning centromere:
        self.__genome__.loc[
            (self.__genome__.end > centromer_loc[0])
            & (self.__genome__.start < centromer_loc[1]), 'GENCODE'
        ] = 'centromere'

    def assign_hetero(self):
        self.__genome__.loc[self.__genome__.GC_ratio.isnull(), 'GENCODE'] = 'heterochromatin'

    def add_colors(self, colors, darkStart, darkMax, dummy=False):
        """
        Colors are also assigned to dummy: only color for the dummy + color for the centromere
        """

        if dummy:
            self.__genome__['color'] = self.__genome__.GENCODE.apply(
                lambda x: colors['centromere'][0] if x == 'centromere' else colors['dummy']
            )

        else:
            self.__genome__['color'] = self.__genome__.apply(color_picker, axis=1, args=(colors))

            self.__genome__['color'] = self.__genome__.apply(
                color_darkener, axis=1, args=(self.__width__, darkStart, darkMax)
            )

    def save_pkl(self, fileName):
        pickle.dump(self.__genome__, open(fileName, "wb"))
