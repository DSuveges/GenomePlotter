import logging
import pandas as pd

from .DataIntegrator import DataIntegrator

class CustomGenePlotter(object):

    GENE_WINDOW = 10_000  # Langth of the flanking on each end of the gene

    def __init__(self, gene_name, config_manager) -> None:
        self.width = config_manager.get_width()
        # Reading gencode data:
        self.gencode_df = pd.read_csv(
            config_manager.get_gencode_file(), compression='infer', sep='\t',
            header=0, dtype={'chr': str, 'start': int, 'end': int, 'type': str}
        )

        # Gene name:
        self.gene_name = gene_name

        # Filter gencode dataset for a given gene:
        self.filtered_gencode = self.filter_gencode_data(gene_name, self.gencode_df)

        # Extract gene coordinates:
        self.chromosome = self.filtered_gencode.iloc[0]['chr']
        self.start = self.filtered_gencode.start.min()
        self.end = self.filtered_gencode.end.max()

        # Get the relevant genome file:
        genome_file = config_manager.get_chromosome_file(self.chromosome)
        self.genome_df = self.load_genome(genome_file)

    @staticmethod
    def filter_gencode_data(gene_name, gencode_df):

        if gene_name.startswith('ENSG'):
            gencode_filtered = gencode_df.loc[gencode_df.gene_id.str.match(gene_name)]
        else:
            gencode_filtered = gencode_df.loc[gencode_df.gene_name == gene_name]

        return gencode_filtered

    def load_genome(self, genome_file):
        """This function loads genome for the gene"""

        genome_df = pd.read_csv(
            genome_file, sep='\t', compression='infer', quotechar='"',
            header=0, dtype={'chr': str, 'start': int, 'end': int, 'GC_ratio': float}
        )

        # Filter genome for the given gene:
        genome_filtered = (
            genome_df
            .loc[
                (genome_df.start >= self.start - self.GENE_WINDOW)
                & (genome_df.end <= self.end + self.GENE_WINDOW)
            ]
        )

        return genome_filtered

    def integrate(self, color_picker) -> None:
        # Initialize data integrator:
        integrator = DataIntegrator(self.genome_df)

        # Convert genomic coordinates to plot coordinates:
        integrator.add_xy_coordinates(self.width)

        # Adding GENCODE annotation to genomic data:
        integrator.add_genes(self.filtered_gencode)

        # Assigning heterocromatic regions:
        integrator.assign_hetero()

        # Assigning colors to individual regions:
        integrator.add_colors(color_picker)

        # Extract integrated data:
        self.integrated_data = integrator.get_data()

    def get_integrated_data(self) -> pd.DataFrame:
        return self.integrated_data
