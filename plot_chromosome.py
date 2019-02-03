#!/usr/bin/env python3

# Importing standard libraries:
import gzip # For reading data
import pandas as pd # data handling
import colorsys # Generate color gradient.
import sys # Read command line arguments:
import argparse # Now command line arugments are properly set.
import pickle

# Importing custom functions:
from functions.helper_fun import *
from functions.color_fun import *
from functions.dataIntegrator import dataIntegrator
from functions.chromosome_plotter import chromosome_plotter
from functions.readConfig import readConfig

# reading parameters:
config = readConfig()

if __name__ == '__main__':

    # Processing command line parameters:
    parser = argparse.ArgumentParser(description='Script to plot genome chunks colored based on GC content and gene annotation. See github: https://github.com/DSuveges/GenomePlotter')
    parser.add_argument('-c', '--chromosome', type = str, help='Selected chromosome to process', required = True)
    parser.add_argument('-w', '--width', type = int, help='Number of chunks in one row.', default = 200)
    parser.add_argument('-p', '--pixel', type = int, help='The size of a plotted chunk in pixels (default: 3).', default = 9)
    parser.add_argument('-s', '--darkStart', type = float, help='Fraction of the width from where the colors start getting darker (default: 0.75).', default = 0.75)
    parser.add_argument('-m', '--darkMax', type = float, help='How dark a pixel can get at the right end of the plot (default: 0.15).', default = 0.15)
    parser.add_argument('-f', '--folder', type = str, help = 'The working directory (default is the current working directory)', default='.')
    parser.add_argument('-t', '--test', type = int, help = 'The number of chunks to be read (by default the whole chromosome is processed.)', default = 0)
    parser.add_argument('-d', '--dummy', type = bool, help = 'If instead of the chunks, a dummy is drawn with identical dimensions', default = False)
    parser.add_argument('--config', type = str, help = 'Specifying json file containing custom configuration', default = 'config.json')

    # Extracting submitted options:
    args = parser.parse_args()
    chromosome = args.chromosome
    width = args.width
    pixel = args.pixel
    darkStart = args.darkStart
    darkMax = args.darkMax
    dummy = args.dummy
    configFile = args.config

    ##
    ## Checking input parameters.... will be implemented later. At this point we want something functional...
    ##

    # Reading custom configuration:
    if configFile:
        config.updateConfig(configFile)

    # Parsing config files:
    colorScheme = config.getColorScheme()
    inputFiles = config.getFiles()

    # Processing colors:
    colors_GENCODE = {
        'centromere'      : linear_gradient(colorScheme['centromere'], n=20),
        'heterochromatin' : linear_gradient(colorScheme['heterochromatin'], finish_hex=colorScheme['heterochromatin'], n=20), # Monochrome, no gradient!
        'intergenic'      : linear_gradient(colorScheme['intergenic'], n=20),
        'exon'            : linear_gradient(colorScheme['exon'], n=20),
        'gene'            : linear_gradient(colorScheme['gene'], n=20)
    } 

    # 
    cytoband_file = inputFiles['cytoband_file']
    Chromosome_file = (inputFiles['Chromosome_file_pattern'] % str(chromosome))
    GWAS_file_loc = inputFiles['GWAS_file_loc']
    GENCODE_file_loc = inputFiles['GENCODE_file_loc']

    # Reading datafile:
    print("[Info] Reading input files...")
    chr_df = pd.read_csv(Chromosome_file, compression='gzip', sep='\t', quotechar='"', header=0, dtype={'chr': str, 'start' : int, 'end' : int, 'GC_ratio': float})
    GENCODE_df = pd.read_csv(GENCODE_file_loc, compression='gzip', sep='\t', header = 0, dtype = {'chr' : str, 'start' : int, 'end': int, 'type' : str })
    cyb_df = pd.read_csv(cytoband_file, compression='gzip', sep='\t', header = 0, dtype = {'chr' : str, 'start' : int, 'end': int, 'name' : str, 'type' : str })

    ##
    ## Integrating cytoband, sequence and gene data:
    ##
    print("[Info] Integrating data...")
    integrator = dataIntegrator(chr_df)
    integrator.add_xy_coordinates(width)
    integrator.add_genes(GENCODE_df)
    integrator.add_centromere(cyb_df)
    integrator.assign_hetero()
    integrator.add_colors(colors_GENCODE, darkStart, darkMax)
    integrator.save_pkl("data/IntegratedData.chr%s.pkl" % chromosome)
    integratedData = integrator.get_data()

    ##
    ## Generate chromosome plot
    ##
    print("[Info] Generating plots...")
    x = chromosome_plotter(integratedData, pixel=pixel)# 
    x.draw_dummy()
    x.add_centromere(cyb_df)
    x.wrap_svg("test.svg")
    x.save_png("data/chr%s.dummy.png" % chromosome)
    x.pickle_data("data/chr%s.dummy.pkl" % chromosome)

    x.draw_chromosome()
    x.add_centromere(cyb_df)
    x.wrap_svg("test.svg")
    x.save_png("data/chr%s.png" % chromosome)
    x.pickle_data("data/chr%s.pkl" % chromosome)