"""A set of functions to integrate parsed input data."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import TYPE_CHECKING, Any

import pandas as pd
import pybedtools
from loguru import logger

if TYPE_CHECKING:
    from genome_plotter.functions.ColorFunctions import ColorPicker


def read_data(file: str, types: dict[str, Any]) -> pd.DataFrame:
    """Read the input data.

    Args:
        file (str): The file to read.
        types (typedict): The types of the columns.

    Returns:
        pd.DataFrame: The data.
    """
    return pd.read_csv(file, compression="gzip", sep="\t", header=0, dtype=types)


def _integrate_one_chromosome(
    chromosome: str,
    output_dir: str,
    gencode_file: str,
    cytoband_file: str,
    dummy: bool,
) -> str:
    """Integrate data for a single chromosome.

    Reads all input files independently so this function is safe to run in a
    worker process without sharing any state with the parent or other workers.

    Args:
        chromosome (str): Chromosome identifier (e.g. '1', 'X', 'MT').
        output_dir (str): Directory containing per-chromosome input files and
            where the integrated output file will be written.
        gencode_file (str): Path to the processed GENCODE annotation file.
        cytoband_file (str): Path to the processed cytoband file.
        dummy (bool): When True, skip gene annotation and use a dummy label.

    Returns:
        str: The chromosome identifier, for use in completion logging.
    """
    gencode_df = read_data(
        gencode_file, {"chr": str, "start": int, "end": int, "type": str}
    )
    cytoband_df = read_data(
        cytoband_file,
        {"chr": str, "start": int, "end": int, "name": str, "type": str},
    )
    chr_df = read_data(
        f"{output_dir}/processed_chr{chromosome}.bed.gz",
        {"chr": str, "start": int, "end": int, "GC_ratio": float},
    )

    integrator = DataIntegrator(chr_df)

    if dummy:
        integrator.add_centromere(cytoband_df)
        integrator.add_dummy()
    else:
        integrator.add_genes(gencode_df)
        integrator.add_centromere(cytoband_df)
        integrator.assign_hetero()

    integrator.save_table(f"{output_dir}/integrated_chr{chromosome}.bed.gz")
    return chromosome


def integrate_data(
    output_dir: str,
    chromosomes: list[str],
    cytoband_file: str,
    gencode_file: str,
    dummy: bool = False,
    max_workers: int | None = None,
) -> None:
    """Integrate the parsed input data for all chromosomes in parallel.

    Args:
        output_dir (str): The directory containing per-chromosome files and
            where integrated outputs will be saved.
        chromosomes (list[str]): Chromosome identifiers to process.
        cytoband_file (str): Path to the processed cytoband file.
        gencode_file (str): Path to the processed GENCODE file.
        dummy (bool): When True, skip gene annotation and fill with a dummy
            label. Defaults to False.
        max_workers (int | None): Number of parallel worker processes. Defaults
            to the number of available CPU cores, capped at len(chromosomes).
    """
    n_workers = min(len(chromosomes), max_workers or os.cpu_count() or 1)
    logger.info(
        f"Integrating {len(chromosomes)} chromosomes using {n_workers} parallel workers."
    )

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_chrom = {
            executor.submit(
                _integrate_one_chromosome,
                chromosome,
                output_dir,
                gencode_file,
                cytoband_file,
                dummy,
            ): chromosome
            for chromosome in chromosomes
        }
        for future in as_completed(future_to_chrom):
            chromosome = future_to_chrom[future]
            future.result()  # re-raises any exception from the worker
            logger.info(f"Integration for chromosome {chromosome} complete.")

    logger.info("Integration complete.")


class DataIntegrator:
    """Assign color for each chunk based on genomic features.

    This class assigns color for each chunk in the chromosome
    based on GC content + GENCODE annotation + darkened fraction.
    """

    __required_columns = ["chr", "start", "end"]

    def __init__(self: DataIntegrator, genome_df: pd.DataFrame) -> None:
        """Initialize the data integrator."""
        self.__genome__ = genome_df.copy()
        self.chromosome_name = genome_df.iloc[0]["chr"]

        logger.info(f"Integrating data on chromosome: {self.chromosome_name}")
        logger.info(
            f"Number of chunks on chromosome {self.chromosome_name}: {len(self.__genome__):,}"
        )

        # Testing columns:
        for col in self.__required_columns:
            if col not in genome_df.columns:
                raise ValueError(
                    f"Manadatory colum: {col} is not found in the provided dataframe."
                )

    def add_xy_coordinates(self: DataIntegrator, width: int) -> None:
        """Calculate chunk coordinates.

        Args:
            width (int): Number of chunks in one row.
        """
        # If the data is already there, we might re-compute:
        if "x" in self.__genome__.columns:
            self.__genome__.drop(columns=["x", "y"], inplace=True)

        # By default, all chunks are written into the same row:
        if not width:
            width = len(self.__genome__)
            self.__genome__.reset_index(drop=True, inplace=True)

        self.__width__ = width

        self.__genome__ = (
            self.__genome__.assign(
                # Get x position of the chunk:
                x=self.__genome__.index.astype(int) % width,
                # Set y position of the chunk:
                y=self.__genome__.index.astype(int) / width,
            )
            # Set proper type:
            .astype({"y": "int32"})
        )
        logger.info(f"Number of chunks in one row: {width:,}")
        logger.info(f"Number of rows: {self.__genome__.y.max():,}")

    def add_genes(self: DataIntegrator, gencode_df: pd.DataFrame) -> None:
        """Add GENCODE annotation for each chunk in the genome.

        Args:
            gencode_df (pd.DataFrame): GENCODE annotation data.
        """
        # If the dataset has already been annotated, we return:
        if "GENCODE" in self.__genome__.columns:
            return None

        # Creating bedtools objects:
        gencode_bed = pybedtools.bedtool.BedTool.from_dataframe(
            (
                gencode_df.loc[gencode_df.chr == self.chromosome_name].rename(
                    columns={"chr": "chrom"}
                )
            )
        )
        chrom_bed = pybedtools.bedtool.BedTool.from_dataframe(
            self.__genome__.rename(columns={"chr": "chrom"})
        )

        # Build column names dynamically so the list matches the actual number of
        # columns produced by the intersection regardless of whether x/y coordinates
        # have been added to the genome DataFrame yet.
        a_names = list(self.__genome__.columns)
        b_names = ["chr_2", "start_2", "end_2", "gene_id", "gene_name", "transcript_id", "type"]
        intersect_names = a_names + b_names

        # Run intersectbed and extract result as dataframe:
        try:
            gencode_intersect = chrom_bed.intersect(gencode_bed, wa=True, wb=True)
            intersect_df = gencode_intersect.to_dataframe(
                header=None,
                names=intersect_names,
            )
        # At this point we should have a dataframe with all the data:
        except Exception as e:
            logger.error("Gencode data:")
            logger.error(gencode_bed.head())
            logger.error("Chromosome data:")
            logger.error(chrom_bed.head())

            raise e

        # pybedtools returns a bare DataFrame() with no columns when the
        # intersection output file is empty (0 bytes), ignoring the names
        # argument.  This happens for chromosomes with no protein-coding
        # GENCODE entries (e.g. MT).  Treat every chunk as intergenic.
        if "start" not in intersect_df.columns:
            logger.info(
                f"No GENCODE features found for chromosome {self.chromosome_name}; "
                "marking all chunks as intergenic."
            )
            self.__genome__["GENCODE"] = "intergenic"
            return

        # Parse out results:
        gencode_chunks = intersect_df.groupby("start").apply(
            lambda x: "exon" if "exon" in x.type.unique() else "gene"
        )
        gencode_chunks.name = "GENCODE"

        # Updating index:
        genome_df = self.__genome__.merge(
            gencode_chunks, left_on="start", right_index=True, how="left"
        ).fillna({"GENCODE": "intergenic"})

        # Adding annotation to df:
        self.__genome__ = genome_df

    def get_data(self: DataIntegrator) -> pd.DataFrame:
        """Get a copy of the integrated data.

        Returns:
            pd.DataFrame: The integrated data.
        """
        return self.__genome__.copy()

    def add_centromere(self: DataIntegrator, cytoband_df: pd.DataFrame) -> None:
        """Add centromere annotation to the genome.

        Args:
            cytoband_df (pd.DataFrame): Cytoband data.
        """
        chromosome = self.__genome__["chr"].iloc[0]
        centromer_loc = cytoband_df.loc[
            (cytoband_df.chr == str(chromosome)) & (cytoband_df.type == "acen"),
            ["start", "end"],
        ]
        centromer_loc = (int(centromer_loc.start.min()), int(centromer_loc.end.max()))

        # If GENCODE column is missing, let's initialize:
        if "GENCODE" not in self.__genome__.columns:
            self.__genome__["GENCODE"] = None

        # Assigning centromere:
        self.__genome__.loc[
            (self.__genome__.end > centromer_loc[0])
            & (self.__genome__.start < centromer_loc[1]),
            "GENCODE",
        ] = "centromere"

    def assign_hetero(self: DataIntegrator) -> None:
        """Assign heterochromatic regions."""
        self.__genome__.loc[self.__genome__.GC_ratio.isnull(), "GENCODE"] = (
            "heterochromatin"
        )

    def add_colors(self: DataIntegrator, color_picker: ColorPicker) -> None:
        """Assigned color to each chunk in the genome.

        Args:
            color_picker (ColorPicker): The color picker.
        """
        self.__genome__["color"] = color_picker.pick_colors_vectorized(self.__genome__)

    def save_table(self: DataIntegrator, file_name: str) -> None:
        """Save the integrated data as a table.

        Args:
            file_name (str): The file name.
        """
        self.__genome__.to_csv(file_name, sep="\t", index=False, compression="infer")

    def add_dummy(self: DataIntegrator) -> None:
        """Assume GENCODE annotation is dummy for all chunks."""
        if "GENCODE" not in self.__genome__.columns:
            self.__genome__ = self.__genome__.assign(GENCODE="dummy")
        else:
            self.__genome__["GENCODE"] = self.__genome__.GENCODE.fillna("dummy")
