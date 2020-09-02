#!/usr/bin/env python3


# Importing standard libraries:
import gzip
import pandas as pd 
import colorsys 
import sys 
import argparse 
import pickle
import os
import logging 


# Importing custom functions:
from functions.helper_fun import *
from functions.color_fun import *
from functions.dataIntegrator import dataIntegrator
from functions.chromosome_plotter import chromosome_plotter
from functions.ConfigManager import ConfigManager
from functions.svg_handler import svg_handler
from functions.gwas_annotator import gwas_annotator

def gwas_annotation_wrapper(config_manager, chromosome):
    # Extract config values:
    color_scheme = config_manager.get_color_scheme()
    gwas_color = color_scheme['gwas_point']
    gwas_file = config_manager.get_gwas_file()
    pixel = config_manager.get_pixel()
    chunk_size = config_manager.get_chunk_size()
    width = config_manager.get_width()
    gwas_file = config_manager.get_gwas_file()

    logging.info(f'Generating GWAS annotation from file: {gwas_file}.')

    gwasAnnot = gwas_annotator(chromosome=chromosome, gwasColor=gwas_color, pixel=pixel, 
                               chunkSize=chunk_size, gwasFile=gwas_file, width=width)
    return gwasAnnot.generateGWAS()

def integrator_wrapper(config_manager, is_dummy=False):
    """
    Opens file, integrates data, returns with integrated dataframe
    """

    # Extracting colors:
    colorScheme = config_manager.get_color_scheme()

    # Processing colors:
    colors_GENCODE = {
        'centromere'      : linear_gradient(colorScheme['centromere'], n=20),
        'heterochromatin' : linear_gradient(colorScheme['heterochromatin'], finish_hex=colorScheme['heterochromatin'], n=20), # Monochrome, no gradient!
        'intergenic'      : linear_gradient(colorScheme['intergenic'], n=20),
        'exon'            : linear_gradient(colorScheme['exon'], n=20),
        'gene'            : linear_gradient(colorScheme['gene'], n=20),
        'dummy'           : colorScheme['dummy']
    } 

    # Extracting file locations from config:
    cytoband_file = config_manager.get_cytoband_file()
    chromosome_file = config_manager.get_chromosome_file(chromosome)
    gencode_file = config_manager.get_gencode_file()

    # Updating config file:
    logging.info('Updating config file.')

    # Reading datafile:
    logging.info('Reading input files.')
    chr_df = pd.read_csv(chromosome_file, compression='gzip', sep='\t', quotechar='"', header=0, dtype={'chr': str, 'start' : int, 'end' : int, 'GC_ratio': float})
    GENCODE_df = pd.read_csv(gencode_file, compression='gzip', sep='\t', header = 0, dtype = {'chr' : str, 'start' : int, 'end': int, 'type' : str })
    cyb_df = pd.read_csv(cytoband_file, compression='gzip', sep='\t', header = 0, dtype = {'chr' : str, 'start' : int, 'end': int, 'name' : str, 'type' : str })
    logging.info(f'Number of genome chunks: {len(chr_df)}')
    logging.info(f'Number of GENCODE annotations in the genome: {len(GENCODE_df)}')
    logging.info(f'Number of cytological bands in the genome: {len(cyb_df)}')

    ##
    ## Integrating cytoband, sequence and gene data:
    ##
    logging.info(f'Integrating data...')

    # Initialize data integrator:
    integrator = dataIntegrator(chr_df)

    # Convert genomic coordinates to plot coordinates:
    integrator.add_xy_coordinates(width)

    # Adding GENCODE annotation to genomic data:
    integrator.add_genes(GENCODE_df)

    # Adding cytological band information to the data:
    integrator.add_centromere(cyb_df)

    # Assigning heterocromatic regions:
    integrator.assign_hetero()

    # Assigning colors to individual regions:
    logging.info(f'Calculating colors for each chunk.')
    integrator.add_colors(colors_GENCODE, darkStart, darkMax, is_dummy)

    # Extract integrated data:
    integratedData = integrator.get_data()

    # Save data for diagnostic purposes:
    integratedData.to_csv('cica.tsv.gz', sep='\t', index=False, compression='infer')

    return integratedData



if __name__ == '__main__':

    # Processing command line parameters:
    parser = argparse.ArgumentParser(description='Script to plot genome chunks colored based on GC content and gene annotation. See github: https://github.com/DSuveges/GenomePlotter')
    parser.add_argument('-c', '--chromosome', type=str, help='Selected chromosome to process', required = True)
    parser.add_argument('-w', '--width', type=int, help='Number of chunks in one row.', default = 200)
    parser.add_argument('-p', '--pixel', type=int, help='The size of a plotted chunk in pixels (default: 3).', default = 9)
    parser.add_argument('-s', '--darkStart', type=float, help='Fraction of the width from where the colors start getting darker (default: 0.75).', default = 0.75)
    parser.add_argument('-m', '--darkMax', type=float, help='How dark a pixel can get at the right end of the plot (default: 0.15).', default = 0.15)
    parser.add_argument('-f', '--folder', type=str, help = 'Folder into which the plots are saved.', required=True)
    parser.add_argument('-t', '--test', type=int, help='The number of chunks to be read (by default the whole chromosome is processed.)', default = 0)
    parser.add_argument('--dummy', help='If instead of the chunks, a dummy is drawn with identical dimensions', action='store_true')
    parser.add_argument('--config', type=str, help='Specifying json file containing custom configuration', required=True)
    parser.add_argument('-l', '--logFile', type=str, help='File into which the logs are generated.', required=False, default='plot_chromosome.log')

    # python plot_chromosome.py -c 22 -w 200 -p 9 -s 0.75 -m 0.15 -f plots/ --config config.json -l logfile.log

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

    # Initialize logger:
    handlers = [ logging.StreamHandler(sys.stdout) ]
    if args.logFile != '':
        handlers.append( logging.FileHandler(filename=args.logFile) )

    # Initialize logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

    logging.info(f'Generating plot for chromosome: {chromosome}')
    logging.info('Processing parameters.')

    # Output file name:
    output_filename = f'{plot_folder}/chr{chromosome}_dummy.png' if dummy else f'{plot_folder}/chr{chromosome}.png'

    # initialize config manager:
    config_manager = ConfigManager(config_file)

    # Set new configuration:
    config_manager.set_width(width)
    config_manager.set_pixel_size(pixel)
    config_manager.set_dark_start(darkStart)
    config_manager.set_dark_max(darkMax)
    config_manager.set_plot_folder(plot_folder)

    # Updating config file:
    config_manager.save_config('pocok.json')

    # Integrating data:
    integratedData = integrator_wrapper(config_manager, dummy)

    ##
    ## Generate chromosome plot
    ##
    logging.info(f'Initializing plot.')
    x = chromosome_plotter(integratedData, pixel=pixel)# 

    if dummy:
        logging.info(f'Generating dummy plot.')
        x.draw_dummy()

    else:
        logging.info(f'Generating plot.')
        x.draw_chromosome()

    # Extract data after plotting:
    plot_width = x.get_plot_with()
    plot_height = x.get_plot_height()
    plot_data = x.return_svg()

    # Initialize svg wrapper object with the returned data:
    chromosomeSvgObject = svg_handler(plot_data, plot_width, plot_height)

    # Generate gwas annitation:
    gwas_annotation = gwas_annotation_wrapper(config_manager, chromosome)
    chromosomeSvgObject.appendSvg(gwas_annotation)

    # Save file:
    logging.info(f'Saving image: {output_filename}')
    chromosomeSvgObject.savePng(output_filename)

    logging.info(f'All done.')

