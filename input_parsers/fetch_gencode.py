"""Module to fetch and parse GENCODE data."""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING

import pandas as pd

from functions.FetchFromFtp import FetchFromFtp

if TYPE_CHECKING:
    from functions.ConfigManager import SourcePrototype

logger = logging.getLogger(__name__)


# Fetch and process Gencode data:
class FetchGencode(FetchFromFtp):
    """Function to fetch gene data from GENCODE ftp"""

    def __init__(self, gencode_parameters: SourcePrototype) -> None:
        """Based on the source related parameters, fetches, parses and saves GENCODE gene info"""

        # Extract parameters:
        self.host = gencode_parameters.host
        self.path = gencode_parameters.path
        self.source_file = gencode_parameters.source_file
        self.processed_file = gencode_parameters.processed_file
        self.arrow_file = gencode_parameters.arrow_file
        self.release_date = None
        self.release = None

        # Initialize host:
        super().__init__(self.host)

    def retrieve_data(self):
        """
        Retrieving release data and association table.
        Then the connection is closed.
        """

        # Get last release:
        file_list = self.fetch_file_list(self.path)
        releases = [
            int(x.replace("release_", ""))
            for x in file_list
            if re.match(r"release_\d+$", x)
        ]
        self.release = max(releases)

        # Get date of last release:
        path_to_last_release = "{}/release_{}/".format(self.path, self.release)
        self.release_date = self.fetch_last_update_date(path_to_last_release)

        # Fetch data:
        self.fetch_tsv(
            path_to_last_release,
            self.source_file.format(self.release),
            skiprows=5,
            header=None,
        )

        # Parse data:
        columns = {
            0: "chr",
            2: "type",
            3: "start",
            4: "end",
            6: "strand",
            8: "annotation",
        }
        self.gencode_raw = self.tsv_data.rename(columns=columns)[columns.values()]

        # Close connection:
        self.close_connection()

        logger.info("GENCODE data retrieved.")
        logger.info(f"Number of rows in the table: {len(self.gencode_raw):,}")

        gene_count = len(self.gencode_raw.loc[self.gencode_raw.type == "gene"])
        logger.info(f"Number of genes in {self.release} release: {gene_count:,}")

    # Processing gwas data:
    def process_gencode_data(self):
        # Parsing gtf annotation:
        logger.info("Parsing GTF annotation.")
        parsed_annotation = self.gencode_raw.annotation.apply(
            lambda annotation: {
                x.strip().split(" ", 1)[0]: x.strip().split(" ", 1)[1].replace('"', "")
                for x in annotation.split(";")
                if x != ""
            }
        )

        # Merging annotation with coordinates:
        gencode_df_updated = self.gencode_raw.merge(
            pd.DataFrame(parsed_annotation.tolist()), left_index=True, right_index=True
        )

        # Drop unparsed annotation column:
        gencode_df_updated.drop(["annotation"], axis=1, inplace=True)

        # Removing gene identifier version:
        gencode_df_updated = gencode_df_updated.assign(
            gene_id=gencode_df_updated.gene_id.str.split(".").str[0]
        )

        # Filtering for protein coding genes:
        gencode_df_updated = gencode_df_updated.loc[
            gencode_df_updated.gene_type == "protein_coding"
        ]
        protein_coding_gene_count = len(
            gencode_df_updated.loc[gencode_df_updated.type == "gene"]
        )
        logger.info(f"Number of protein coding genes: {protein_coding_gene_count:,}")

        # Updating types:
        gencode_df_updated["start"] = gencode_df_updated.start.astype(int)
        gencode_df_updated["end"] = gencode_df_updated.end.astype(int)

        # Adding length to all features:
        gencode_df_updated = gencode_df_updated.assign(
            length=lambda row: row["end"] - row["start"]
        )

        # Initialize empty dataframes:
        processed = pd.DataFrame(
            columns=[
                "chr",
                "start",
                "end",
                "type",
                "gene_id",
                "gene_name",
                "transcript_id",
            ]
        )
        arrow_data = pd.DataFrame(
            columns=["chr", "start", "end", "strand", "type", "gene_id", "gene_name"]
        )

        logger.info(
            "Generate exon/intron annotations for the canonical transcripts for each gene... (it will take a while.)"
        )

        for (gene_id,), features in gencode_df_updated.groupby(["gene_id"]):
            # Selecting protein coding transcript identifiers:
            transcripts = features.loc[
                (features.type == "transcript")
                & (features.transcript_type == "protein_coding")
            ]

            # If no protein coding transcript is found, we skip gene:
            if len(transcripts) == 0:
                continue

            # Adding the length of the CDS to all transcripts:
            cds_length = transcripts.transcript_id.apply(
                lambda t_id: features.loc[
                    (features.type == "CDS") & (features.transcript_id == t_id)
                ].length.sum()
            )
            transcripts.insert(2, "cds_length", cds_length)

            # Get canonical transcript and properties:
            canonical_transcript_id = self.get_canonical_transcript(transcripts)
            [start, end] = (
                transcripts.loc[
                    transcripts.transcript_id == canonical_transcript_id,
                    ["start", "end"],
                ]
                .iloc[0]
                .tolist()
            )

            # Get data for the arrow plot:
            arrow_part = features.loc[
                (features.transcript_id == canonical_transcript_id)
                & (features.type.isin(["CDS", "UTR"])),
                ["chr", "start", "end", "strand", "type", "gene_id", "gene_name"],
            ]
            arrow_data = pd.concat([arrow_data, arrow_part])

            # Generate exon-intron splice:
            gene_df = self.generate_exon_intron_structure(
                gene_id,
                canonical_transcript_id,
                start,
                end,
                features.loc[
                    (features.transcript_id == canonical_transcript_id)
                    & (features.type == "exon")
                ],
            )
            # appending to existing data:
            processed = pd.concat([processed, gene_df])

        # Saving data:
        self.processed = processed
        self.arrow_data = arrow_data

    # Save data
    def save_gencode_data(self, data_dir):
        self.processed = self.processed[
            ["chr", "start", "end", "gene_id", "gene_name", "transcript_id", "type"]
        ]

        gencode_output_filename = f"{data_dir}/{self.processed_file}"
        self.processed.to_csv(
            gencode_output_filename, sep="\t", compression="infer", index=False
        )

        gencode_arrow_filename = f"{data_dir}/{self.arrow_file}"
        self.arrow_data.to_csv(
            gencode_arrow_filename, sep="\t", compression="infer", index=False
        )
        logger.info(f"GENCODE data saved into {gencode_arrow_filename}.")

    # Extract release date:
    def get_release_date(self):
        return self.release_date

    # Extract release date:
    def get_release(self):
        return self.release

    @staticmethod
    def get_canonical_transcript(transcripts):
        # Selecting canonical transcript following Ensembl guidelines:
        # http://www.ensembl.org/Help/Glossary

        # First choice: selecting the longest ccds transcript:
        if transcripts.ccdsid.notna().any():
            longest_CDS = transcripts.loc[~transcripts.ccdsid.isna()].cds_length.max()
            canonical_transcript_id = transcripts.loc[
                longest_CDS == transcripts.cds_length
            ].transcript_id.tolist()[0]

        # Secong choice: selecting the longest havana transcript:
        elif transcripts.havana_transcript.notna().any():
            longest_CDS = transcripts.loc[
                ~transcripts.havana_transcript.isna()
            ].cds_length.max()
            canonical_transcript_id = transcripts.loc[
                longest_CDS == transcripts.cds_length
            ].transcript_id.tolist()[0]

        # Third choice: selecting the longest protein coding transcript:
        else:
            longest_CDS = transcripts.cds_length.max()
            canonical_transcript_id = transcripts.loc[
                longest_CDS == transcripts.cds_length
            ].transcript_id.tolist()[0]

        return canonical_transcript_id

    @staticmethod
    def generate_exon_intron_structure(gene_id, transcript_id, start, end, exons):
        # Looping through all exons generate coordinates of the corresponding intron(s)
        new_features = []
        for _, values in (
            exons[["start", "end"]]
            .sort_values(by="start", axis=0, ascending=True)
            .iterrows()
        ):
            es = values["start"]
            ee = values["end"]

            # Adding intron:
            if es > start:
                new_features.append({"start": start, "end": es, "type": "intron"})

            # Adding exon to feature:
            new_features.append({"start": es, "end": ee, "type": "exon"})

            # updating start position:
            start = ee

        # Adding final exon if exists:
        if start < end:
            new_features.append({"start": start, "end": end, "type": "intron"})

        # Generate dataframe + add extra features:
        df = pd.DataFrame(new_features)
        df["gene_id"] = gene_id
        df["transcript_id"] = transcript_id
        df["chr"] = exons.chr.tolist()[0].replace("chr", "")
        df["gene_name"] = exons.gene_name.tolist()[0]

        return df
