from __future__ import annotations

import logging

import pandas as pd

from .DataIntegrator import DataIntegrator


class CustomGeneIntegrator(object):
    def __init__(self, query, config_manager) -> None:
        self.width = config_manager.get_width()
        self.gene_window = config_manager.get_custom_gene_window()

        logging.info(f"Generating integrated dataset for gene: {query}")

        # Reading gencode data:
        self.gencode_df = pd.read_csv(
            config_manager.get_gencode_file(),
            compression="infer",
            sep="\t",
            header=0,
            dtype={"chr": str, "start": int, "end": int, "type": str},
        )

        # Filter gencode dataset for a given gene:
        filtered_gencode = self.__filter_gencode_data(query, self.gencode_df)

        if len(filtered_gencode) == 0:
            raise ValueError(
                f"The gene {query} cound not be found in the GENCODE database."
            )

        # Report what we have:
        self.gene_name = filtered_gencode.iloc[0]["gene_name"]
        self.gene_id = filtered_gencode.iloc[0]["gene_id"]
        logging.info(
            f"Gene name: {self.gene_name }, Ensembl gene identifier: {self.gene_id}"
        )
        logging.info(
            f"Number of gencode feature for this gene: {len(filtered_gencode):,}"
        )

        # Extract gene coordinates:
        self.chromosome = filtered_gencode.iloc[0]["chr"]
        self.start = filtered_gencode.start.min()
        self.end = filtered_gencode.end.max()
        self.filtered_gencode = filtered_gencode

        logging.info(f"Genomic coordinates: {self.chromosome}:{self.start}-{self.end}")

        # Get the relevant genome file:
        genome_file = config_manager.get_chromosome_file(self.chromosome)
        genome_df = self.__load_genome(genome_file)
        logging.info(f"Number of genomic chunks for this gene: {len(genome_df)}")

        self.genome_df = genome_df

    @staticmethod
    def __filter_gencode_data(gene_name, gencode_df):
        if gene_name.startswith("ENSG"):
            gencode_filtered = gencode_df.loc[gencode_df.gene_id.str.match(gene_name)]
        else:
            gencode_filtered = gencode_df.loc[gencode_df.gene_name == gene_name]

        return gencode_filtered

    def get_gencode_data(self):
        return self.filtered_gencode

    def __load_genome(self, genome_file):
        """This function loads genome for the gene"""

        genome_df = pd.read_csv(
            genome_file,
            sep="\t",
            compression="infer",
            quotechar='"',
            header=0,
            dtype={"chr": str, "start": int, "end": int, "GC_ratio": float},
        )

        # Filter genome for the given gene:
        genome_filtered = genome_df.loc[
            (genome_df.start >= self.start - self.gene_window)
            & (genome_df.end <= self.end + self.gene_window)
        ]

        return genome_filtered

    def integrate(self, color_picker) -> None:
        # Initialize data integrator:
        integrator = DataIntegrator(self.genome_df)

        # Convert genomic coordinates to plot coordinates:
        integrator.add_xy_coordinates()

        # Adding GENCODE annotation to genomic data:
        integrator.add_genes(self.filtered_gencode)

        # Assigning heterocromatic regions:
        integrator.assign_hetero()

        # Assigning colors to individual regions:
        integrator.add_colors(color_picker)

        # Extract integrated data:
        integrated_data = integrator.get_data()
        self.integrated_data = integrated_data.assign(
            x=lambda df: df.x - integrated_data.x.min(),
            y=lambda df: df.y - integrated_data.y.min(),
        ).reset_index(drop=True)

    def get_integrated_data(self) -> pd.DataFrame:
        return self.integrated_data

    def get_data_width(self) -> int:
        return self.integrated_data.x.max()


class GenerateArrowPlot(object):
    CHUNK_SVG = '<rect x="{}" y="0" width="{}" height="{}" style="stroke-width:1;stroke:{};fill:{}" />\n'
    LINE_SVG = '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="1" />\n'
    ARROW_SVG = '<polygon points="{}" style="fill:{};stroke:{};stroke-width:1" />\n'

    def __init__(self, config_manager, arrow_width):
        self.arrow_data = pd.read_csv(
            config_manager.get_gencode_arrow_file(), sep="\t", compression="infer"
        )

        # Extract colors:
        arrow_colors = config_manager.get_arrow_colors()
        self.line_color = arrow_colors["line_color"]
        self.utr_color = arrow_colors["utr_color"]
        self.cds_color = arrow_colors["cds_color"]

        # Extract other values
        self.chunk_size = config_manager.basic_parameters.chunk_size
        self.pixel_size = config_manager.plot_parameters.pixel_size

        self.arrow_width = arrow_width

    def generate_arrow_polt(self, gene_name):
        # Get CDS and UTR for a given gene:
        gene_df = (
            self.arrow_data.loc[lambda df: df.gene_name == gene_name]
            .sort_values("start")
            .reset_index(drop=True)
            .drop(["chr", "gene_name", "gene_id"], axis=1)
        )

        # If nothing found raise:
        assert len(gene_df) > 1, f"Gene name ({gene_name}) was not found in the data."

        # Scaling genomic coordinates to screen coordinates:
        gene_df = gene_df.assign(
            relative_start=lambda row: (row["start"] - gene_df.start.min())
            / self.chunk_size
            * self.pixel_size,
            relative_end=lambda row: (row["end"] - gene_df.start.min())
            / self.chunk_size
            * self.pixel_size,
        ).assign(length=lambda row: row["relative_end"] - row["relative_start"])

        # Calculate intron midpoints:
        gene_df["next_start"] = gene_df.relative_start.drop(0).reset_index(drop=True)

        gene_df["midpoint"] = (
            gene_df.next_start - gene_df.relative_end
        ) / 2 + gene_df.relative_end
        self.df = gene_df.head()

        # Generate boxes:
        boxes = gene_df.apply(self.draw_box, axis=1).dropna()

        # Generate lines connecting boxes:
        lines = (
            gene_df.loc[gene_df.midpoint.notna()]
            .apply(self.draw_lines, axis=1)
            .dropna()
        )

        # Generate closing arrow:
        arrow = self.draw_arrow(gene_df)

        # Concatenate into single string:
        svg_string = "".join(lines.to_list() + boxes.to_list()) + arrow
        # svg_object = svg_handler(svg_string, gene_df.relative_end.max(), self.arrow_width)

        return svg_string, gene_df.relative_end.max()

    def draw_arrow(self, df):
        arrow_width = self.arrow_width

        # Get orientation:
        row = df.iloc[0] if df.iloc[0]["strand"] == "-" else df.tail(1).iloc[0]

        if row["strand"] == "+":
            coordinates = [
                (row["relative_start"], 0),
                (row["relative_end"], 0),
                (row["relative_end"], -arrow_width / 2),
                (row["relative_end"] + arrow_width, arrow_width / 2),
                (row["relative_end"], arrow_width * 1.5),
                (row["relative_end"], arrow_width),
                (row["relative_start"], arrow_width),
            ]
        else:
            coordinates = [
                (row["relative_start"], 0),
                (row["relative_end"], 0),
                (row["relative_end"], arrow_width),
                (row["relative_start"], arrow_width),
                (row["relative_start"], 1.5 * arrow_width),
                (row["relative_start"] - arrow_width, arrow_width / 2),
                (row["relative_start"], -arrow_width / 2),
            ]

        coordinate_str = " ".join([f"{x[0]},{x[1]}" for x in coordinates])

        return self.ARROW_SVG.format(coordinate_str, self.utr_color, self.line_color)

    def draw_box(self, row):
        if row["type"] == "UTR":
            fill_color = self.utr_color
        elif row["type"] == "CDS":
            fill_color = self.cds_color
        else:
            fill_color = "black"

        return self.CHUNK_SVG.format(
            row["relative_start"],
            row["length"],
            self.arrow_width,
            self.line_color,
            fill_color,
        )

    def draw_lines(self, row):
        if row["next_start"] - row["relative_end"] < 2:
            return None

        # Two bits are required:
        connector = self.LINE_SVG.format(
            row["relative_end"],
            0,
            row["midpoint"],
            -self.arrow_width / 2,
            self.line_color,
        )
        connector += self.LINE_SVG.format(
            row["midpoint"],
            -self.arrow_width / 2,
            row["next_start"],
            0,
            self.line_color,
        )

        return connector
