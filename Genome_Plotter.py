#!/Users/ds26/anaconda/bin/python

# v.2.1 Last modified: 2018.01.10
    # Adding command line parameter support.
    # All applied parameters can be passed as command line arguments.
    # Flexible margins around plot.
    # Adding support for cytobands
    # Adding support for custom associations and their plot on figure.

# v.2.0 Last modified: 2018.01.06
    # Helper functions are moved to separate files.
    # Use genomic coordinates to find x,y coordinates on the plot.
    # The x and y coordinates can be fixed.

# TODO:
    # Some of the code can be written in oop. Especially the plotting.
    # Implement the addition of the custom markers.

# Importing standard libraries:
import gzip # For reading data
import pybedtools # For finding overlap between our chunks and genes and GWAS signals
import pandas as pd # data handling
import numpy as np # Working with large arrays
import colorsys # Generate color gradient.
import cairosvg # Converting svg to image
import os.path # checking if the datafies are already there.
import sys #
import argparse # Now command line arugments are properly set.

# Importing custom functions:
from functions.helper_fun import *
from functions.color_fun import *
from functions.SVG_plot import SVG_plot

# Defining colors for different structures:
colors_GENCODE = {
    'centromere'      : linear_gradient('#9393FF', n=20),
    'heterochromatin' : linear_gradient('#F9D2C2', finish_hex='#ffc6af', n=20), # Monochrome, no gradient!
    'intergenic'      : linear_gradient('#A3E0D1', n=20),
    'exon'            : linear_gradient('#FFD326', n=20),
    'gene'            : linear_gradient('#6CB8CC', n=20)
}

# Adding command line parameters:
parser = argparse.ArgumentParser(description='Script to plot chromosomes with elaborate annotations. See github: https://github.com/DSuveges/GenomePlotter')
parser.add_argument('-c', '--chromosome', type = str, help='Selected chromosome to process', required = True)
parser.add_argument('-d', '--dimension', type = int, help='Fixed dimension (height of width) of the plot (200 chunks by default).', default = 200)
parser.add_argument('-a', '--axis', type = int, help='The fixed axis of the plot (1 - width, 2 - height, 1 by default)', default = 1)
parser.add_argument('-p', '--pixel', type = int, help='The size of a plotted chunk in pixel (default: 3).', default = 3)
parser.add_argument('-s', '--darkStart', type = float, help='Fraction of the width from where the colors start getting darker (default:  0.75).', default = 0.75)
parser.add_argument('-m', '--darkMax', type = float, help='How dark a pixel can get at the right end of the plot (default: 0.15).', default = 0.15)
parser.add_argument('-f', '--folder', type = str, help = 'The working directory (default is the current working directory)', default='.')
parser.add_argument('-t', '--test', type = int, help = 'The number of chunks to be read (by default the whole chromosome is processed.)', default = 0)
parser.add_argument('--write-pickle', type = int, help = 'The resulting svg object is saved into a serialized file for further analysis.', default = 0)

# Extracting submitted options:
args = parser.parse_args()
chromosome = args.chromosome
dimension = args.dimension
axis = args.axis
pixel = args.pixel
darkStart = args.darkStart
darkMax = args.darkMax

# Understanding the provided axis:
fixed_dim = "Width"
if axis == 2:
    fixed_dim = "Height"

# Print report:
print("[Info %s] Processing chromosome: %s." % (get_now(), chromosome))
print("[Info %s] Fixed dimension: %s, length: %s chunks." % (get_now(), fixed_dim, dimension))

try:
    workingDir = args.folder
except:
    workingDir = os.getcwd()
print("[Info %s] Working directory: %s." % (get_now(), workingDir))

# Checking input files:
cytoband_file = workingDir + '/data/cytoBand.GRCh38.bed.bgz'
Chromosome_file_loc = workingDir + "/data/Processed_chr%s.bed.gz"
GWAS_file_loc = workingDir + "/data/processed_GWAS.bed.gz"
GENCODE_file_loc = workingDir + "/data/GENCODE.merged.bed.gz"
outputDir = workingDir + "/plots"

if not os.path.isfile(GWAS_file_loc):
    print ("GWAS file is missing. Run `prepare_data.sh` first!")
if not os.path.isfile(GENCODE_file_loc):
    print ("GENCODE file is missing. Run `prepare_data.sh` first!")
if not os.path.isfile(Chromosome_file_loc % 11):
    print ("Processed chromosome files are missing! Run `prepare_data.sh` first!")

# Creating plot folder if does not exist:
if not os.path.isdir(outputDir):
    print ("[Info %s] Creating output folder: %s." %( get_now(), outputDir))
    os.mkdir(outputDir)
else:
    print ("[Info %s] Output folder (%s) already exists." %( get_now(), outputDir))

# Reading datafile:
dataFile = Chromosome_file_loc % chromosome
print("[Info %s] reading file %s... " % (get_now(), dataFile))
chr_dataf = pd.read_csv(dataFile, compression='gzip', sep='\t', quotechar='"')

# If we want to do test, we can restrict the dataframe to n lines:
if args.test:
    chr_dataf = chr_dataf.ix[0:args.test]
    print("[Info %s] Test mode is on. The first %s rows will be kept." %(get_now(), args.test))

# the start position will be the new index:
chr_dataf = chr_dataf.set_index('start', drop = False)

chunk_size = chr_dataf.head(1).apply(lambda x: x['end'] - x['start'], axis = 1).tolist()[0]
min_pos  = chr_dataf.start.min()
print("[Info %s] Number of chunks read: %s" %(get_now(), chr_dataf.shape[0]))

# Get the width (number of chunks plotted in one row)
width = dimension # We are dealing with fixed column numbers
if axis == 2:  width = int(chr_dataf.shape[0]/dimension) # We are dealing with fixed number of rows.

# Adding plot coordinates to the dataframe:
print("[Info %s] Chunk size: %s bp, plot width: %s" %(get_now(), chunk_size, width))
print("[Info %s] Calculating plot coordinates for each chunk." %(get_now()))
chr_dataf = chr_dataf.apply(generate_xy, args = (min_pos, chunk_size, width), axis =1)
print(chr_dataf.head())

##
## Get intersecting GENCODE features:
##

# Reading GENCODE file:
print("[Info %s] Reading GENCODE file (%s)..." % (get_now(), GENCODE_file_loc))
GENCODE_bed = pybedtools.BedTool(GENCODE_file_loc)

# Run intersectbed:
chr_data_bed = pybedtools.BedTool(dataFile)

print("[Info %s] Selecting intersecting features." %(get_now()))
GencodeIntersect = chr_data_bed.intersect(GENCODE_bed, wa = True, wb = True)
GC_INT = pybedtools.bedtool.BedTool.to_dataframe(GencodeIntersect)

# Assign features to each chunk:
print("[Info %s] Assign intersecting features to each chunk." % (get_now()))
GENCODE_chunks = GC_INT.groupby('start').apply(lambda x: 'exon' if 'exon' in x.thickEnd.unique() else 'gene' )
GENCODE_chunks.names = "GENCODE"

chr_dataf['GENCODE'] = 'intergenic' # By default, all GENCODE values are intergenic
chr_dataf.GENCODE.update(GENCODE_chunks) # This value will be overwritten if overlaps with exon or gene
print(chr_dataf.head())

# Reading cytoband file:
print("[Info %s] Opening and processing cytoband file..." %(get_now()))
cyb_df = pd.read_csv(cytoband_file, compression='gzip', sep='\t')
cyb_df = cyb_df.loc[cyb_df.chr == chromosome] # Selecting only the relevant rows

# Centromeres actually won't be used.... these rows will be deleted:
centromer_loc = cyb_df.loc[(cyb_df.chr == chromosome) & (cyb_df.type == 'acen'),['start', 'end']]
centromer_loc = (int(centromer_loc.start.min()),
                 int(centromer_loc.end.max()))

# Calculating proper y coordinate for both the start and end position of each cytoband (x coordinate won't be used:):
cyb_df = cyb_df.apply(generate_xy, axis = 1, args = (min_pos, chunk_size, width), y = 'y1')
cyb_df = cyb_df.apply(generate_xy, axis = 1, args = (min_pos, chunk_size, width), position_column = 'end', y = 'y2')

print(cyb_df.head())

# Assigning centromere:
chr_dataf.loc[(chr_dataf.end > centromer_loc[0]) & (chr_dataf.start < centromer_loc[1]) , 'GENCODE'] = 'centromere'

# Assigning heterocromatin where no GC ratio is available:
chr_dataf.loc[chr_dataf.GC_ratio.isnull(), 'GENCODE'] = 'heterochromatin'

##
## Get colors based on GENCODE feature:
##
print("[Info %s] Assigning colors to each chunk." % (get_now()))
chr_dataf['color'] = chr_dataf.apply(lambda x: colors_GENCODE[x['GENCODE']][int(x[3]*20)] if not np.isnan(x[3]) else colors_GENCODE['heterochromatin'][0], axis = 1)

# Get the colors darker based on column number:
print("[Info %s] Apply darkness filter." % (get_now()))
chr_dataf['color'] = chr_dataf.apply(color_darkener, axis = 1, args=(width, darkStart, darkMax))
print(chr_dataf.head())

##
## Overlapping GWAS signals:
##
print("[Info %s] Opening GWAS file.." %(get_now()))
GWAS_file = pybedtools.BedTool(GWAS_file_loc)
full_GWAS = pybedtools.bedtool.BedTool.to_dataframe(GWAS_file)
GWAS_df = full_GWAS[full_GWAS.chrom == str(chromosome)]
print("[Info %s] Calculating plot coordinates for overlapping GWAS sigals." %(get_now()))
GWAS_df = GWAS_df.apply(generate_xy, args = (min_pos, chunk_size, width), axis =1)

##
## Let's draw the plot:
##
# Adding some custom annotation:

outputFileName = outputDir +("/chr%s.w.%s.c.%s" %(chromosome, int(width), int(chunk_size)))
plot = SVG_plot(chr_dataf.x.max(), chr_dataf.y.max(), pixel, margins = [250, 100, 100, 100])

print("[Info %s] Drawing chunks (%s of them)... might take a while to complete." % (get_now(), chr_dataf.shape[0]))
chr_dataf.apply(plot.draw_chunk, axis = 1)

print("[Info %s] Adding GWAS hits." %(get_now()))
GWAS_df.apply(plot.draw_GWAS, axis = 1 )

print("[Info %s] Marking centromere." %(get_now()))
plot.mark_centromere(int(centromer_loc[0]/chunk_size), int(centromer_loc[1]/chunk_size))

GWAS_df.sample(frac = 0.20).apply(plot.add_assoc, axis = 1 )

print("[Info %s] Adding cytoband ruler." %(get_now()))
cyb_df.loc[cyb_df.chr == chromosome].apply(plot.draw_cytoband, axis = 1)

print("[Info %s] Saving svg file.." % (get_now()))
plot.save_svg(outputFileName + ".svg")

print("[Info %s] Saving png file.." % (get_now()))
plot.save_png(outputFileName + ".png")
