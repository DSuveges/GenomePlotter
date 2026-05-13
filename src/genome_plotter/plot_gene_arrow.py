"""Plotting script for gene arrow diagrams."""

from __future__ import annotations

import argparse
import os

import pandas as pd
import yaml
from loguru import logger

from genome_plotter import LOG_FORMAT
from genome_plotter.functions.ColorFunctions import ColorPicker
from genome_plotter.functions.ConfigManager import Config
from genome_plotter.functions.CustomGenePlotter import (
    CustomGeneIntegrator,
    GenerateArrowPlot,
)
from genome_plotter.functions.svg_handler import SvgHandler


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate an arrow plot (gene structure diagram) for a given gene. "
        "See github: https://github.com/DSuveges/GenomePlotter"
    )
    parser.add_argument(
        "-g",
        "--gene",
        help="Gene name (e.g. BRCA1)",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-d",
        "--data_folder",
        help="Path to the data folder containing processed files",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Path to config JSON file",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file basename (defaults to {gene_name}_arrow)",
        required=False,
        type=str,
        default=None,
    )
    parser.add_argument(
        "-s",
        "--scale-factor",
        help="Scaling factor applied to all plot dimensions (default: 30)",
        type=int,
        default=30,
    )

    return parser.parse_args()


def main() -> None:
    """Entry point for the plot-gene-arrow CLI command."""
    args = parse_arguments()

    gene_name = args.gene
    data_folder = os.path.abspath(args.data_folder)
    config_file = args.config
    output_basename = args.output if args.output else f"{gene_name}_arrow"
    scale_factor = args.scale_factor

    # Initialise logger:
    logger.add("genome_plotter.log", level="DEBUG", format=LOG_FORMAT)

    logger.info(f"Generating arrow plot for gene: {gene_name}")

    # Initialise configuration:
    with open(config_file) as f:
        try:
            config_manager = Config.model_validate(yaml.safe_load(f))
        except yaml.YAMLError as e:
            raise ValueError(
                f"The provided config file ({config_file}) is not a valid YAML file."
            ) from e

    # Set the data folder on the config:
    config_manager.basic_parameters.data_folder = data_folder

    # Generate arrow plot:
    logger.info("Initializing arrow plot generator.")
    arrow_plotter = GenerateArrowPlot(
        pd.read_csv(
            config_manager.get_gencode_arrow_file(), sep="\t", compression="infer"
        ),
        config_manager,
        scale_factor=scale_factor,
    )

    logger.info(f"Generating arrow plot for gene: {gene_name}")
    svg_arrow, width = arrow_plotter.generate_arrow_polt(gene_name)

    # Integrate genome chunk data and generate the chunk layer:
    logger.info("Integrating genome chunk data.")
    color_picker = ColorPicker(
        config_manager.color_schema.chromosome_colors, count=30
    )
    integrator = CustomGeneIntegrator(gene_name, config_manager)
    integrator.integrate(color_picker)
    chunk_extension = 10
    svg_chunks = arrow_plotter.generate_chunk_svg(
        integrator.get_integrated_data(), extension=chunk_extension
    )

    # Layout: chunk layer on top, arrow below with a gap equal to one pixel_size:
    px = arrow_plotter.pixel_size
    gap = px
    arrow_y_offset = px + gap
    svg_arrow_shifted = f'<g transform="translate(0 {arrow_y_offset})">\n{svg_arrow}\n</g>\n'

    # The extension adds extra chunks left and right of the gene; expand canvas accordingly.
    # Left extension is handled by shifting translate_x so negative-x chunks land inside the canvas.
    extension_px = chunk_extension * px
    canvas_height = arrow_y_offset + px * 3
    svg_obj = SvgHandler(svg_chunks + svg_arrow_shifted, width + extension_px, canvas_height)
    svg_obj.group(translate=(px + extension_px, px))

    # Save as SVG and PNG:
    svg_filename = f"{output_basename}.svg"
    png_filename = f"{output_basename}.png"

    logger.info(f"Saving SVG: {svg_filename}")
    svg_obj.save_svg(svg_filename)

    logger.info(f"Saving PNG: {png_filename}")
    svg_obj.save_png(png_filename)

    logger.info("All done.")


if __name__ == "__main__":
    main()
