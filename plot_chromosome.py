#!/usr/bin/env python3

import argparse # Now command line arugments are properly set.



if __name__ == '__main__':

    # Processing command line parameters:
    parser = argparse.ArgumentParser(description='Script to plot genome chunks colored based on GC content and gene annotation. See github: https://github.com/DSuveges/GenomePlotter')
    parser.add_argument('-c', '--chromosome', type = str, help='Selected chromosome to process', required = True)
    parser.add_argument('-w', '--width', type = int, help='Number of chunks in one row.', default = 200)
    parser.add_argument('-p', '--pixel', type = int, help='The size of a plotted chunk in pixels (default: 3).', default = 3)
    parser.add_argument('-s', '--darkStart', type = float, help='Fraction of the width from where the colors start getting darker (default: 0.75).', default = 0.75)
    parser.add_argument('-m', '--darkMax', type = float, help='How dark a pixel can get at the right end of the plot (default: 0.15).', default = 0.15)
    parser.add_argument('-f', '--folder', type = str, help = 'The working directory (default is the current working directory)', default='.')
    parser.add_argument('-t', '--test', type = int, help = 'The number of chunks to be read (by default the whole chromosome is processed.)', default = 0)
    parser.add_argument('-d', '--dummy', type = bool, help = 'If instead of the chunks, a dummy is drawn with identical dimensions', default = False)

    # Extracting submitted options:
    args = parser.parse_args()
    chromosome = args.chromosome
    width = args.width
    pixel = args.pixel
    darkStart = args.darkStart
    darkMax = args.darkMax
    dummy = args.dummy

    ##
    ## Checking input parameters.... will be implemented later. At this point we want something functional...
    ##

    # 