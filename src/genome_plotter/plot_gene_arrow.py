"""Plotting script for gene arrow diagrams."""

from __future__ import annotations

import argparse
import json
import logging.config
import os

import pandas as pd
import yaml

from genome_plotter.functions.ConfigManager import Config
from genome_plotter.functions.CustomGenePlotter import GenerateArrowPlot
from genome_plotter.functions.svg_handler import svg_handler

logger = logging.getLogger(__name__)


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
        "--data-folder",
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
        "-w",
        "--arrow-width",
        help="Arrow height in pixels (default: 10)",
        type=int,
        default=10,
    )

    return parser.parse_args()


def main() -> None:
    """Entry point for the plot-gene-arrow CLI command."""
    args = parse_arguments()

    gene_name = args.gene
    data_folder = os.path.abspath(args.data_folder)
    config_file = args.config
    arrow_width = args.arrow_width
    output_basename = args.output if args.output else f"{gene_name}_arrow"

    # Initialise logger:
    logger_config_path = os.path.join(os.path.dirname(__file__), "logger_config.yaml")
    with open(logger_config_path, "r") as stream:
        logger_config = yaml.safe_load(stream)

    logging.config.dictConfig(logger_config)
    logger = logging.getLogger(__name__)

    logger.info(f"Generating arrow plot for gene: {gene_name}")

    # Initialise configuration:
    with open(config_file) as f:
        try:
            config_manager = Config(**json.load(f))
        except json.decoder.JSONDecodeError:
            raise ValueError(
                f"The provided config file ({config_file}) is not a valid JSON file."
            )

    # Set the data folder on the config:
    config_manager.basic_parameters.data_folder = data_folder

    # Generate arrow plot:
    logger.info("Initializing arrow plot generator.")
    arrow_plotter = GenerateArrowPlot(
        pd.read_csv(
            config_manager.get_gencode_arrow_file(), sep="\t", compression="infer"
        ),
        config_manager,
        arrow_width,
    )

    logger.info(f"Generating arrow plot for gene: {gene_name}")
    svg_string, width = arrow_plotter.generate_arrow_polt(gene_name)

    # Wrap in svg_handler:
    svg_obj = svg_handler(svg_string, width, arrow_width)

    # Save as SVG and PNG:
    svg_filename = f"{output_basename}.svg"
    png_filename = f"{output_basename}.png"

    logger.info(f"Saving SVG: {svg_filename}")
    svg_obj.saveSvg(svg_filename)

    logger.info(f"Saving PNG: {png_filename}")
    svg_obj.savePng(png_filename)

    logger.info("All done.")


if __name__ == "__main__":
    main()
