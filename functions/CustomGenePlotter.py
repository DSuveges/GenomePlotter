import logging
import pandas as pd

from .DataIntegrator import DataIntegrator

class CustomGeneIntegrator(object):

    GENE_WINDOW = 10_000  # Langth of the flanking on each end of the gene

    def __init__(self, query, config_manager) -> None:
        self.width = config_manager.get_width()

        logging.info(f'Generating integrated dataset for gene: {query}')

        # Reading gencode data:
        self.gencode_df = pd.read_csv(
            config_manager.get_gencode_file(), compression='infer', sep='\t',
            header=0, dtype={'chr': str, 'start': int, 'end': int, 'type': str}
        )

        # Filter gencode dataset for a given gene:
        filtered_gencode = self.filter_gencode_data(query, self.gencode_df)

        if len(filtered_gencode) == 0:
            raise ValueError(f'The gene {query} cound not be found in the GENCODE database.')

        # Report what we have:
        self.gene_name = filtered_gencode.iloc[0]['gene_name']
        self.gene_id = filtered_gencode.iloc[0]['gene_id']
        logging.info(f'Gene name: {self.gene_name }, Ensembl gene identifier: {self.gene_id}')
        logging.info(f'Number of gencode feature for this gene: {len(filtered_gencode)}')

        # Extract gene coordinates:
        self.chromosome = filtered_gencode.iloc[0]['chr']
        self.start = filtered_gencode.start.min()
        self.end = filtered_gencode.end.max()
        self.filtered_gencode = filtered_gencode

        logging.info(f'Genomic coordinates: {self.chromosome}:{self.start}-{self.end}')

        # Get the relevant genome file:
        genome_file = config_manager.get_chromosome_file(self.chromosome)
        genome_df = self.load_genome(genome_file)
        logging.info(f'Number of genomic chunks for this gene: {len(genome_df)}')

        self.genome_df = genome_df

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
