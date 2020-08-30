#!/usr/bin/env python3

# Importing standard libraries:
import gzip
import pandas as pd 
import colorsys 
import sys 
import argparse 
import pickle
import os

# Importing custom functions:
from functions.helper_fun import *
from functions.color_fun import *
# from functions.dataIntegrator import dataIntegrator
from functions.chromosome_plotter import chromosome_plotter
from functions.ConfigManager import ConfigManager

if __name__ == '__main__':

    # Processing command line parameters:
    parser = argparse.ArgumentParser(description='Script to plot genome chunks colored based on GC content and gene annotation. See github: https://github.com/DSuveges/GenomePlotter')
    parser.add_argument('-c', '--chromosome', type = str, help='Selected chromosome to process', required = True)
    parser.add_argument('-w', '--width', type = int, help='Number of chunks in one row.', default = 200)
    parser.add_argument('-p', '--pixel', type = int, help='The size of a plotted chunk in pixels (default: 3).', default = 9)
    parser.add_argument('-s', '--darkStart', type = float, help='Fraction of the width from where the colors start getting darker (default: 0.75).', default = 0.75)
    parser.add_argument('-m', '--darkMax', type = float, help='How dark a pixel can get at the right end of the plot (default: 0.15).', default = 0.15)
    parser.add_argument('-f', '--folder', type = str, help = 'Folder into which the plots are saved.', required=True)
    parser.add_argument('-t', '--test', type = int, help = 'The number of chunks to be read (by default the whole chromosome is processed.)', default = 0)
    parser.add_argument('-d', '--dummy', help = 'If instead of the chunks, a dummy is drawn with identical dimensions', action='store_true')
    parser.add_argument('--config', type = str, help = 'Specifying json file containing custom configuration', required=True)

    # python plot_chromosome.py -c 22 -w 200 -p 9 -s 0.75 -m 0.15 -f plots/ --config config.json

    # Extracting submitted options:
    args = parser.parse_args()
    chromosome = args.chromosome
    width = args.width
    pixel = args.pixel
    darkStart = args.darkStart
    darkMax = args.darkMax
    dummy = args.dummy
    config_file = args.config
    plot_folder = os.path.abspath(args.folder)

    # initialize config manager:
    cm = ConfigManager(config_file)

    # Set new configuration:
    cm.set_width(width)
    cm.set_pixel_size(pixel)
    cm.set_dark_start(darkStart)
    cm.set_dark_max(darkMax)
    cm.set_plot_folder(plot_folder)

    # Parsing config files:
    colorScheme = cm.get_color_scheme()
    data_folder = cm.get_data_folder()

    # Processing colors:
    colors_GENCODE = {
        'centromere'      : linear_gradient(colorScheme['centromere'], n=20),
        'heterochromatin' : linear_gradient(colorScheme['heterochromatin'], finish_hex=colorScheme['heterochromatin'], n=20), # Monochrome, no gradient!
        'intergenic'      : linear_gradient(colorScheme['intergenic'], n=20),
        'exon'            : linear_gradient(colorScheme['exon'], n=20),
        'gene'            : linear_gradient(colorScheme['gene'], n=20)
    } 

    # The config manager returns the path to all required files and also does a check on the existence:
    cytoband_file = cm.get_cytoband_file()
    chromosome_file = cm.get_chromosome_file(chromosome)
    gwas_file = cm.get_gwas_file()
    gencode_file = cm.get_gencode_file()

    # Updating config file:
    cm.save_config('pocok.json')

    # Reading datafile:
    print("[Info] Reading input files...")
    chr_df = pd.read_csv(chromosome_file, compression='gzip', sep='\t', quotechar='"', header=0, dtype={'chr': str, 'start' : int, 'end' : int, 'GC_ratio': float})
    GENCODE_df = pd.read_csv(gencode_file, compression='gzip', sep='\t', header = 0, dtype = {'chr' : str, 'start' : int, 'end': int, 'type' : str })
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

    # ##
    # ## Generate chromosome plot
    # ##
    # print("[Info] Generating plots...")
    # x = chromosome_plotter(integratedData, pixel=pixel)# 
    # x.draw_dummy()
    # x.add_centromere(cyb_df)
    # # x.wrap_svg("test.svg")
    # # x.save_png("data/chr%s.dummy.png" % chromosome)
    # x.pickle_data("data/chr%s.dummy.pkl" % chromosome)

    # # If only dummy is requested, exiting:
    # if dummy:
    #     quit()

    # x.draw_chromosome()
    # x.add_centromere(cyb_df)
    # # x.wrap_svg("test.svg")
    # # x.save_png("data/chr%s.png" % chromosome)
    # x.pickle_data("data/chr%s.pkl" % chromosome)


