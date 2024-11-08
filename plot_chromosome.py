"""Plotting script for chromosome data."""

from __future__ import annotations

import argparse
import json
import logging.config
import os
from dataclasses import asdict

import pandas as pd
import yaml

from functions.ChromosomePlotter import ChromosomePlotter

# Importing custom functions:
from functions.ColorFunctions import ColorPicker
from functions.ConfigManager import Config
from functions.CytobandAnnotator import CytobandAnnotator, get_centromere_position
from functions.DataIntegrator import DataIntegrator
from functions.GeneAnnotator import GeneAnnotator
from functions.GwasAnnotator import gwas_annotator
from functions.svg_handler import svg_handler


def genes_annotation_wrapper(
    config_manager: Config, chromosome: str, height: int, gene_filename: str
) -> GeneAnnotator:
    """Wrapper function to generate gene annotation.

    Args:
        config_manager (Config): Configuration manager object.
        chromosome (str): Chromosome to process.
        height (int): Height of the plot.
        gene_filename (str): File containing gene annotations.

    Returns:
        GeneAnnotator: Gene annotator object.
    """
    # Extractig config values:
    pixel = config_manager.plot_parameters.pixel_size
    chunk_size = config_manager.basic_parameters.chunk_size
    width = config_manager.plot_parameters.width
    centromere_position = get_centromere_position(
        config_manager.get_cytoband_file(), chromosome
    )

    # gene_file, centromerePosition, chromosome, chunk_size, pixel, height, width
    gene_annotator_object = GeneAnnotator(
        gene_filename, centromere_position, chromosome, chunk_size, pixel, height, width
    )
    return gene_annotator_object


def cytoband_annotation_wrapper(
    config_manager: Config, chromosome: str
) -> CytobandAnnotator:
    """Wrapper function to generate cytoband annotation.

    Args:
        config_manager (Config): Configuration manager object.
        chromosome (str): Chromosome to process.

    Returns:
        CytobandAnnotator: Cytoband annotator object.
    """
    # Extract config values:
    color_scheme = asdict(config_manager.color_schema)
    pixel = config_manager.plot_parameters.pixel_size
    chunk_size = config_manager.basic_parameters.chunk_size
    width = config_manager.plot_parameters.width
    cytoband_file = config_manager.get_cytoband_file()

    logger.info(f"Generating cytological band from file: {cytoband_file}.")
    # pixel, chromosome, bandFile, chunkSize, width, cytbandColors
    cytoband_annot = CytobandAnnotator(
        pixel,
        chromosome,
        cytoband_file,
        chunk_size,
        width,
        color_scheme["cytoband_colors"],
    )
    cytoband_annot.generate_bands()
    return cytoband_annot


def gwas_annotation_wrapper(config_manager: Config, chromosome: str) -> gwas_annotator:
    # Extract config values:
    color_scheme = asdict(config_manager.color_schema)
    gwas_color = color_scheme["gwas_point"]
    pixel = config_manager.plot_parameters.pixel_size
    chunk_size = config_manager.basic_parameters.chunk_size
    width = config_manager.plot_parameters.width
    gwas_file = config_manager.get_gwas_file()

    logger.info(f"Generating GWAS annotation from file: {gwas_file}.")

    gwasAnnot = gwas_annotator(
        chromosome=chromosome,
        gwas_color=gwas_color,
        pixel=pixel,
        chunk_size=chunk_size,
        gwas_file=gwas_file,
        width=width,
    )
    return gwasAnnot.generate_gwas()


def integrator_wrapper(
    config_manager: Config, dummy: bool, chromosome: str
) -> pd.DataFrame:
    """Integrate input data.

    Args:
        config_manager (Config): Configuration object.
        dummy (bool): Flag to indicate if dummy data should be generated.
        chromosome (str): Chromosome to process.

    Returns:
        pd.DataFrame: Integrated data.
    """

    # Extracting parameters from config:
    cytoband_file = config_manager.get_cytoband_file()
    chromosome_file = config_manager.get_chromosome_file(chromosome)
    gencode_file = config_manager.get_gencode_file()
    dark_start = config_manager.plot_parameters.dark_start
    dark_max = config_manager.plot_parameters.dark_max
    color_map = config_manager.color_schema.chromosome_colors

    width = config_manager.plot_parameters.width

    # Initialize color picker object:
    color_picker = ColorPicker(
        color_map, width=width, dark_threshold=dark_start, dark_max=dark_max, count=30
    )

    # Reading datafiles:
    logger.info("Reading input files.")
    chr_df = pd.read_csv(
        chromosome_file,
        compression="gzip",
        sep="\t",
        quotechar='"',
        header=0,
        dtype={"chr": str, "start": int, "end": int, "GC_ratio": float},
    )
    GENCODE_df = pd.read_csv(
        gencode_file,
        compression="gzip",
        sep="\t",
        header=0,
        dtype={"chr": str, "start": int, "end": int, "type": str},
    )
    cyb_df = pd.read_csv(
        cytoband_file,
        compression="gzip",
        sep="\t",
        header=0,
        dtype={"chr": str, "start": int, "end": int, "name": str, "type": str},
    )
    logger.info(f"Number of genome chunks: {len(chr_df):,}")
    logger.info(f"Number of GENCODE annotations in the genome: {len(GENCODE_df):,}")
    logger.info(f"Number of cytological bands in the genome: {len(cyb_df):,}")

    # Integrating cytoband, sequence and gene data:
    logger.info("Integrating data...")

    # Initialize data integrator:
    integrator = DataIntegrator(chr_df)

    # Convert genomic coordinates to plot coordinates:
    integrator.add_xy_coordinates(width)

    # Downstream processing depends on dummy status:
    if dummy:
        # Adding cytological band information to the data:
        integrator.add_centromere(cyb_df)

        # Adding dummy GENCODE annotation to genomic data:
        integrator.add_dummy()

    else:
        # Adding GENCODE annotation to genomic data:
        integrator.add_genes(GENCODE_df)

        # Adding cytological band information to the data:
        integrator.add_centromere(cyb_df)

        # Assigning heterocromatic regions:
        integrator.assign_hetero()

    # Assigning colors to individual regions:
    integrator.add_colors(color_picker)

    # Extract integrated data:
    integratedData = integrator.get_data()

    # Save data for diagnostic purposes:
    integratedData.to_csv("cica.tsv.gz", sep="\t", index=False, compression="infer")

    return integratedData


def parse_arguments() -> argparse.Namespace:
    # Processing command line parameters:
    parser = argparse.ArgumentParser(
        description="Script to plot genome chunks colored based on GC content and gene annotation. \
            See github: https://github.com/DSuveges/GenomePlotter"
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        help="Selected chromosome to process",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-w", "--width", help="Number of chunks in one row.", type=int, default=200
    )
    parser.add_argument(
        "-p",
        "--pixel",
        help="The size of a plotted chunk in pixels (default: 3).",
        type=int,
        default=9,
    )
    parser.add_argument(
        "-s",
        "--darkStart",
        help="Fraction of the width from where the colors start getting darker (default: 0.75).",
        type=float,
        default=0.75,
    )
    parser.add_argument(
        "-m",
        "--darkMax",
        help="How dark a pixel can get at the right end of the plot (default: 0.15).",
        type=float,
        default=0.15,
    )
    parser.add_argument(
        "-f",
        "--folder",
        help="Folder into which the plots are saved.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--textFile",
        help="Flag to indicate if svg file should also be saved.",
        action="store_true",
    )
    parser.add_argument(
        "-g",
        "--geneFile",
        help="A .bed file with genes to add to the chromosome.",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-t",
        "--test",
        help="The number of chunks to be read (by default the whole chromosome is processed.)",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--dummy",
        help="If instead of the chunks, a dummy is drawn with identical dimensions",
        action="store_true",
    )
    parser.add_argument(
        "--config",
        help="Specifying json file containing custom configuration",
        type=str,
        required=True,
    )

    return parser.parse_args()


if __name__ == "__main__":
    # Extracting submitted options:
    args = parse_arguments()

    chromosome = args.chromosome
    width = args.width
    pixel = args.pixel
    dark_start = args.darkStart
    dark_max = args.darkMax
    dummy = args.dummy
    config_file = args.config
    plot_folder = os.path.abspath(args.folder)

    # Initialise logger:
    with open("logger_config.yaml", "r") as stream:
        logger_config = yaml.safe_load(stream)

    logging.config.dictConfig(logger_config)
    logger = logging.getLogger(__name__)

    # Loading config:
    with open(args.config) as f:
        try:
            configuration = Config(**json.load(f))
        except json.decoder.JSONDecodeError:
            raise ValueError(
                f"The provided config file ({args.config}) is not a valid JSON file."
            )

    # Reporting parameters:
    logger.info(f"Generating plot for chromosome: {chromosome}")
    logger.info("Processing parameters.")
    logger.info(f"Number of chunks in one row: {width}")
    logger.info(f"Pixel size: {pixel}")
    logger.info(f"Dark start: {dark_start}, dark max: {dark_max}")
    logger.info(f"Plot is going to be saved into folder: {plot_folder}")
    if dummy:
        logger.info("Creating dummy without chromosome details.")

    # Output file name:
    output_filename = (
        f"{plot_folder}/chr{chromosome}_dummy.png"
        if dummy
        else f"{plot_folder}/chr{chromosome}.png"
    )

    # Initilise configuration:
    with open(config_file) as f:
        try:
            config_manager = Config(**json.load(f))
        except json.decoder.JSONDecodeError:
            raise ValueError(
                f"The provided config file ({config_file}) is not a valid JSON file."
            )

    # Set new configuration:
    config_manager.plot_parameters.width = width
    config_manager.plot_parameters.pixel_size = pixel
    config_manager.plot_parameters.dark_start = dark_start
    config_manager.plot_parameters.dark_max = dark_max
    config_manager.basic_parameters.plot_folder = plot_folder

    # Updating config file:
    logger.info(f"Updating config file: {config_file}")
    config_manager.save(config_file)

    # Integrating data:
    logger.info("Integrating data...")
    integratedData = integrator_wrapper(config_manager, dummy, chromosome)

    # Generate chromosome plot
    logger.info("Initializing plot.")
    x = ChromosomePlotter(integratedData, pixel=pixel)

    if dummy:
        logger.info("Generating dummy plot.")
        x.draw_dummy()
    else:
        logger.info("Generating plot.")
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

    # Generate cytoband annotation:
    cytoband_annotation_obj = cytoband_annotation_wrapper(config_manager, chromosome)
    (cyb_width, cyb_height) = cytoband_annotation_obj.get_dimensions()
    cyb_svg = svg_handler(cytoband_annotation_obj.return_svg(), cyb_width, cyb_height)

    chromosomeSvgObject.group(translate=(cyb_width, 0))
    chromosomeSvgObject.mergeSvg(cyb_svg)

    if args.geneFile:
        # Get centromere position:
        centromerePos = get_centromere_position(
            config_manager.get_cytoband_file(), chromosome
        )

        # Get plot dimension:
        plot_height = chromosomeSvgObject.getHeight()

        # Create gene annotator object:
        gene_annot = genes_annotation_wrapper(
            config_manager, chromosome, plot_height, args.geneFile
        )

        # Generate annotation:
        gene_annot.generate_gene_annotation()

        dimensions = gene_annot.get_dimensions()

        gene_svg = svg_handler(
            width=800,
            height=abs(dimensions[0]) + dimensions[1],
            svg_string=gene_annot.get_annotation(),
        )

        # Move svg object if there are negative values:
        if abs(dimensions[0]) > 0:
            gene_svg.group(translate=(0, abs(dimensions[0])))
            chromosomeSvgObject.group(translate=(0, abs(dimensions[0])))

        # Merging together with the chromosome:
        gene_svg.group(translate=(chromosomeSvgObject.getWidth(), 0))
        chromosomeSvgObject.mergeSvg(gene_svg)

    # Save file:
    logger.info(f"Saving image: {output_filename}")
    chromosomeSvgObject.savePng(output_filename)

    if args.textFile:
        logger.info(f'Saving svg file: {output_filename.replace("png","svg")}')
        chromosomeSvgObject.saveSvg(output_filename.replace("png", "svg"))

    logger.info("All done.")
