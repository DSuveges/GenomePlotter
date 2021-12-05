import sys

import pandas as pd
import logging

from functions.CustomGenePlotter import CustomGeneIntegrator, GenerateArrowPlot, get_translation_start_chunk
from functions.ConfigManager import ConfigManager
from functions.svg_handler import svg_handler
from functions.ColorFunctions import ColorPicker
from functions.ChromosomePlotter import ChromosomePlotter


def main(gene_name):
    # Initialize logger:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Config manager:
    cm = ConfigManager('config.updated.json')

    # Color picker:
    cp = ColorPicker(cm.get_chromosome_colors(), 0.9, 0.9, 20, 4000)

    # Get custom gene integrator:
    cgi = CustomGeneIntegrator(gene_name, cm)
    cgi.integrate(cp)

    # Get chromosome plotter:
    chrp = ChromosomePlotter(cgi.get_integrated_data(), pixel=18)
    logging.info(cgi.get_integrated_data())
    chrp.draw_chromosome()

    # Get arrow plot:
    arrow_width = 50
    gap = GenerateArrowPlot(cm, arrow_width=arrow_width)

    svg_string, gene_width = gap.generate_arrow_polt(gene_name)

    # Generate svg object:
    chunks_svg = svg_handler(chrp.return_svg(), cgi.get_data_width() * 18, 20, background='white')

    arrow_svg = svg_handler(svg_string, gene_width, 100, background='white')
    arrow_svg.group(translate=(10_000 / 450 * 18, 100))

    # get translation start position and chunk:
    translation_start = gap.get_translation_start()
    translation_start_chunk = get_translation_start_chunk(translation_start, cgi.get_integrated_data())

    # Get translation start chunk:
    logging.info(translation_start_chunk)

    # Merging chunks with the arrow:
    chunks_svg.mergeSvg(arrow_svg)
    chunks_svg.group(translate=(10, 10))
    chunks_svg.savePng(f'gene_{gene_name}.png')


if __name__ == '__main__':

    gene_name = sys.argv[1]

    main(gene_name)
