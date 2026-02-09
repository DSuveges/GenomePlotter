"""A set of functions to integrate parsed input data."""

from __future__ import annotations

import logging
import pickle
from typing import TYPE_CHECKING, Any

import pandas as pd
import pybedtools

if TYPE_CHECKING:
    from genome_plotter.functions.ColorFunctions import ColorPicker


logger = logging.getLogger(__name__)


def read_data(file: str, types: dict[str, Any]) -> pd.DataFrame:
    """Read the input data.

    Args:
        file (str): The file to read.
        types (typedict): The types of the columns.

    Returns:
        pd.DataFrame: The data.
    """
    return pd.read_csv(file, compression="gzip", sep="\t", header=0, dtype=types)


def integrate_data(
    output_dir: str,
    chromosomes: list[str],
    cytoband_file: str,
    gencode_file: str,
    dummy: bool = False,
) -> None:
    """Integrate the parsed input data.

    Args:
        output_dir (str): The directory to save the data to.
        chromosomes (list[str]): The chromosomes.
        cytoband_file (str): The cytobands file.
        gencode_file (str): The gencode file.
        dummy (bool, optional): Whether to use dummy data. Defaults to False.
    """
    logger.info("Integrating parsed data.")
    # Read cytobands:
    gencode_df = read_data(
        gencode_file, {"chr": str, "start": int, "end": int, "type": str}
    )
    cyb_df = read_data(
        cytoband_file,
        {"chr": str, "start": int, "end": int, "name": str, "type": str},
    )
    logger.info(f"Number of GENCODE annotations in the genome: {len(gencode_df):,}")
    logger.info(f"Number of cytological bands in the genome: {len(cyb_df):,}")

    # Iterate over chromosomes:
    for chromosome in chromosomes:
        logger.info(f"Integrating data for chromosome: {chromosome}")
        # Reading chromosome data:
        chr_df = read_data(
            f"{output_dir}/processed_chr{chromosome}.bed.gz",
            {"chr": str, "start": int, "end": int, "GC_ratio": float},
        )

        # Initialize data integrator:
        integrator = DataIntegrator(chr_df)

        # Downstream processing depends on dummy status:
        if dummy:
            # Adding cytological band information to the data:
            integrator.add_centromere(cyb_df)

            # Adding dummy GENCODE annotation to genomic data:
            integrator.add_dummy()

        else:
            # Adding GENCODE annotation to genomic data:
            integrator.add_genes(gencode_df)

            # Adding cytological band information to the data:
            integrator.add_centromere(cyb_df)

            # Assigning heterocromatic regions:
            integrator.assign_hetero()

        # Integration is complete:
        integrator.save_table(f"{output_dir}/integrated_chr{chromosome}.bed.gz")
        logger.info(f"Integration for chromosome {chromosome} is complete.")

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

        # Run intersectbed and extract result as dataframe:
        try:
            gencode_intersect = chrom_bed.intersect(gencode_bed, wa=True, wb=True)
            intersect_df = gencode_intersect.to_dataframe(
                header=None,
                names=[
                    "chr",
                    "start",
                    "end",
                    "GC_ratio",
                    "x",
                    "y",
                    "chr_2",
                    "start_2",
                    "end_2",
                    "gene_id",
                    "gene_name",
                    "transcript_id",
                    "type",
                ],
            )
        # At this point we should have a dataframe with all the data:
        except Exception as e:
            logger.error("Gencode data:")
            logger.error(gencode_bed.head())
            logger.error("Chromosome data:")
            logger.error(chrom_bed.head())

            raise e

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
        chromosome = self.__genome__.chr[1]
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
        self.__genome__["color"] = self.__genome__.apply(
            color_picker.pick_color, axis=1
        )

    def save_pkl(self: DataIntegrator, file_name: str) -> None:
        """Save the integrated data as a pickle file.

        Args:
            file_name (str): The file name.
        """
        pickle.dump(self.__genome__, open(file_name, "wb"))

    def save_table(self: DataIntegrator, file_name: str) -> None:
        """Save the integrated data as a table.

        Args:
            file_name (str): The file name.
        """
        self.__genome__.to_csv(file_name, sep="\t", index=False, compression="infer")

    def add_dummy(self: DataIntegrator) -> None:
        """Assume GENCODE annotation is dummy for all chunks."""
        if "GENCODE" not in self.__genome__.columns:
            self.__genome__ = self.__genome__.assign("GENCODE", "dummy")
        else:
            self.__genome__["GENCODE"] = self.__genome__.GENCODE.fillna("dummy")
