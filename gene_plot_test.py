import sys
import logging

from functions.CustomGenePlotter import CustomGeneIntegrator
from functions.ConfigManager import ConfigManager
from functions.ColorFunctions import ColorPicker


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout)]
)

cm = ConfigManager('config.updated.json')

# Get color schema:
colorScheme = {k: v for k, v in cm.get_color_scheme().items() if isinstance(v, str) and v.startswith('#')}

cp = ColorPicker(
    colorScheme,
    width=cm.get_width(),
    dark_threshold=cm.get_dark_start(),
    dark_max=cm.get_dark_max(),
    count=30
)

gene_name = sys.argv[1]
cgp = CustomGeneIntegrator(gene_name, cm)
cgp.integrate(cp)

df = cgp.get_integrated_data()
df.to_csv('custom_gene.tsv', sep='\t', index=False)

gencode_df = cgp.get_gencode_data()
gencode_df.to_csv('custom_gene_gencode.tsv', sep='\t', index=False)
