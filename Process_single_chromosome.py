#!/usr/bin/env python

##
## This script helps prototyping 
## 

# Importing standard libraries:
import gzip # For reading data
import pandas as pd # data handling
import colorsys # Generate color gradient.
import sys # Read command line arguments:
#import argparse # Now command line arugments are properly set.
import pickle

# Importing custom functions:
from functions.helper_fun import *
from functions.color_fun import *
from functions.dataIntegrator import dataIntegrator
from functions.chromosome_plotter import chromosome_plotter

##
## Change this dictionary to change colors of the plot:
##
colors_GENCODE = {
    'centromere'      : linear_gradient('#9393FF', n=20),
    'heterochromatin' : linear_gradient('#F9D2C2', finish_hex='#ffc6af', n=20), # Monochrome, no gradient!
    'intergenic'      : linear_gradient('#A3E0D1', n=20),
    'exon'            : linear_gradient('#FFD326', n=20),
    'gene'            : linear_gradient('#6CB8CC', n=20)
} 

workingDir = '/Users/dsuveges/Project/GenomePlotter'

# Input files are quite hardcoded:
cytoband_file = workingDir + '/data/cytoBand.GRCh38.bed.bgz'
Chromosome_file_loc = workingDir + "/data/Processed_chr%s.bed.gz"
GWAS_file_loc = workingDir + "/data/processed_GWAS.bed.gz"
GENCODE_file_loc = workingDir + "/data/GENCODE.merged.bed.gz"

# Chromosome is read from the command line:
chromosome = sys.argv[1]

# Other parameters are hardcoded, no need to change:
width = 200
darkStart = 0.75
darkMax = 0.15
pixel = 9

# Reading datafile:
print("[Info] Reading input files...")
chr_df = pd.read_csv(Chromosome_file_loc % chromosome, compression='gzip', sep='\t', quotechar='"')
GENCODE_df = pd.read_csv(GENCODE_file_loc, compression='gzip', sep='\t')
cyb_df = pd.read_csv(cytoband_file, compression='gzip', sep='\t')

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

# pickleFileName = "data/IntegratedData.chr%s.pkl" % chromosome
# Reading pickle:
# pickle_in = open(pickleFileName,"rb")
# integratedData = pickle.load(pickle_in)

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

