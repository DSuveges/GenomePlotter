import sys # Read command line arguments:
import argparse # Now command line arugments are properly set.
import pickle

# Custom functions:
from functions.readConfig import readConfig
from functions.svg_handler import svg_handler
from functions.gwas_annotator import gwas_annotator
from functions.cytoband_annotator import cytoband_annotator

# Annotator modules: all annotation types are added using a dedicated module eg. GWAS 

# Main code block:
if __name__ == '__main__':
    """
    This script adds annotation to the chromosome svg file.
    """
    
    # Parsing command line arguments:
    parser = argparse.ArgumentParser(description = "This script adds annotation to a chromosome plot.")
    parser.add_argument('--chromosome', type=str, help='Name of the Chromosome.', required = True)
    parser.add_argument('--dummy', default=False, help='Dummy is used or not.', action='store_true')
    args = parser.parse_args()

    # Parsing input file name:
    chromosome = args.chromosome
    fileName = 'data/chr{}.dummy.pkl'.format(chromosome) if args.dummy else 'data/chr{}.pkl'.format(chromosome)

    # Reading data:
    with open(fileName, 'rb') as f:
        data = pickle.load(f)

    ## extractig data:
    chunkSize = data['chunk_size']
    svg_data = data['data']
    screenWidth = data['width']
    screenHeight = data['height']
    pixel = data['pixel']
    width = screenWidth / pixel

    ## Reading configuration:
    config = readConfig()
    filesConfig = config.getFiles()
    colorConfig = config.getColorScheme()

    ###
    ### Adding GWAS hits:
    ###
    gwasColor = colorConfig['gwas_point']
    gwasFile = filesConfig['GWAS_file_loc']
    gwasAnnot = gwas_annotator(chromosome=chromosome, gwasColor=gwasColor, pixel=pixel, 
                               chunkSize=chunkSize, gwasFile=gwasFile, width=width)
    points = gwasAnnot.generateGWAS()
    
    # Adding points to the svg image:
    # Initialize the first object
    chromosomeSvgObject = svg_handler(svg_data, screenWidth, screenHeight)

    # Adding Annotation to the chromosome:
    chromosomeSvgObject.appendSvg(points)

    # Saving png:
    chromosomeSvgObject.savePng("plots/chr{}_GWAS.png".format(chromosome))

    ###
    ### Generate cytobands:
    ###
    # Reading cytoband file:
    bandFile = filesConfig['cytoband_file']
    cytbandColors = colorConfig['cytobandColors']
    cyb = cytoband_annotator(bandFile=bandFile,chromosome=chromosome,cytbandColors=cytbandColors,chunkSize=chunkSize, 
                             pixel=pixel, width=width)
    cyb.generate_bands()
    (cyb_width, cyb_height) = cyb.get_dimensions()
    cyb_svg = svg_handler(cyb.return_svg(), cyb_width, cyb_height)


    ##
    ## Merging svgs
    ##
    chromosomeSvgObject.group(translate = (cyb_width,0))
    chromosomeSvgObject.mergeSvg(cyb_svg)
    chromosomeSvgObject.savePng("plots/chr{}_GWAS_cytoband.png".format(chromosome))

