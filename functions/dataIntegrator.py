import pandas as pd 
import numpy as np
import pickle
import pybedtools 
from .color_fun import *

class dataIntegrator(object):

    def __init__(self, genome_df):
        # Needs to be tested if that's a dataframe and all the columns are added!
        self.__genome__ = genome_df.copy()

    def add_xy_coordinates(self, width):
        self.__width__ = width
        self.__genome__['x'] = self.__genome__.index % width
        self.__genome__['y'] = [int(x / width) for x in self.__genome__.index ]

    def add_genes(self, GENCODE_df):
        # Creating bedtools objects:
        gencode_bed = pybedtools.bedtool.BedTool.from_dataframe(GENCODE_df)
        chrom_bed = pybedtools.bedtool.BedTool.from_dataframe(self.__genome__[['chr','start', 'end']])

        # Run intersectbed and extract result as dataframe:
        GencodeIntersect = chrom_bed.intersect(gencode_bed, wa = True, wb = True)
        intersect_df = pybedtools.bedtool.BedTool.to_dataframe(GencodeIntersect)

        # Parse out results:
        GENCODE_chunks = intersect_df.groupby('start').apply(lambda x: 'exon' if 'exon' in x.thickStart.unique() else 'gene' )
        GENCODE_chunks.names = "GENCODE"

        # Updating index:
        GENCODE_chunks.index = self.__genome__.loc[self.__genome__.start.isin(GENCODE_chunks.index)].index

        # Adding annotation to df:
        self.__genome__['GENCODE'] = 'intergenic' # By default, all GENCODE values are intergenic
        self.__genome__.GENCODE.update(GENCODE_chunks) # This value will be overwritten if overlaps with exon or gene

    def get_data(self):
        return(self.__genome__.copy())

    def add_centromere(self, cytoband_df):
        chromosome = self.__genome__.chr[1]
        centromer_loc = cytoband_df.loc[(cytoband_df.chr == str(chromosome)) & (cytoband_df.type == 'acen'),['start', 'end']]
        centromer_loc = (int(centromer_loc.start.min()),
                         int(centromer_loc.end.max()))

        # Assigning centromere:
        self.__genome__.loc[(self.__genome__.end > centromer_loc[0]) & (self.__genome__.start < centromer_loc[1]) , 'GENCODE'] = 'centromere'

    def assign_hetero(self):
        self.__genome__.loc[self.__genome__.GC_ratio.isnull(), 'GENCODE'] = 'heterochromatin'

    def add_colors(self, colors, darkStart, darkMax):
        self.__genome__['color'] = self.__genome__.apply(lambda x: colors[x['GENCODE']][int(x[3]*20)] if not np.isnan(x[3]) else colors['heterochromatin'][0], axis = 1)
        self.__genome__['color'] = self.__genome__.apply(color_darkener, axis = 1, args=(self.__width__, darkStart, darkMax))

    def save_pkl(self, fileName):
        pickle.dump( self.__genome__, open( fileName, "wb" ) )
