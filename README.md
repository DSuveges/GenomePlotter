## GenomePlotter

[![codecov](https://codecov.io/gh/DSuveges/GenomePlotter/branch/master/graph/badge.svg)](https://codecov.io/gh/DSuveges/GenomePlotter)

The motivation behind this project was to create a scientifically correct visualization of the human genome with showing the location of genes, exons, already published genome-wide associations and many more. I chose a heat-map kind of visualization that allows the representation of the underlying chemical properties of the DNA ([GC content](https://en.wikipedia.org/wiki/GC-content)) as well.

### Approach

The project is divided into several parts:

1. A script (`Prepare_data.py`) that downloads and pre-processes the source data files with genome sequence, GWAS signals, gene annotation and the [cytological bands](https://en.wikipedia.org/wiki/G_banding).
2. A script (`plot_chromosome.py`), that integrates sequence data and gene annotation into a single svg (and png).

### Requirements

**Required command line tools:**

- [cairo graphics library](https://www.cairographics.org/download/)
- [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) v2.27 or above
- [uv](https://docs.astral.sh/uv/)

### Installing package:

```bash
# Cloning repo:
git clone https://github.com/DSuveges/GenomePlotter

# Creating environment:
cd GenomePlotter
uv sync --all-extras
```

### Source data:

As input, the most recent datasets are pulled from the respected sources as part of the data preparation process. All source data are mapped to the GRCh38 build of the human genome.

* **The sequence of the human genome** is dowloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html) (checking for the most recent version).
* **Genome wide association signals** most recent version of the NHGRI-EBI [GWAS catalog](https://www.ebi.ac.uk/gwas/) (checking the most recent version).
* **Gene annotation** the most recent gene coordinates are downloaded from [GENCODE](http://www.gencodegenes.org/releases/current.html).
* **Canonical transcripts** of protein coding genes are defined according to [Ensembl](https://www.ensembl.org/Help/Glossary).
* **Cytological bands** coordinates fetched from the Ensembl [REST API](http://rest.ensembl.org/).

More information on the sources can be found in the `config.json` configuration file.

### Step 1 - Pre-processing

```bash
python Prepare_data.py -d data_folder/ -c config.json -s 450 -t 0.5
```

help output:

```
usage: prepare_data.py [-h] -d DATADIR -c CONFIG -s CHUNKSIZE -t TOLERANCE

This script fetches and parses input data for the genome plotter project

optional arguments:
  -h, --help            show this help message and exit
  -d DATADIR, --dataDir DATADIR
                        Folder into which the input data and the temporary files will be saved
  -c CONFIG, --config CONFIG
                        JSON file with configuration data
  -s CHUNKSIZE, --chunkSize CHUNKSIZE
                        Chunk size to pool genomic sequence in basepairs.
  -t TOLERANCE, --tolerance TOLERANCE
                        Fraction of a chunk that cannot be N.
```

* *<DATADIR>* folder into which the files are going to be saved.
* *<CONFIG>* JSON file containing the project level configuration. Will be used for multiple scripts
* *<LOGFILE>* information on the run is saved here.
* *<CHUNKSIZE>* the length of non-overlapping window used to pool together to calculate [GC content](https://en.wikipedia.org/wiki/GC-content). In basepairs.
* *<TOLERANCE>* Ns are discarded from the GC content calculation. This float (ranging from 0-1) shows the maximum of Ns in a chunk tolerated. Chunks with too high N ratio is considered as heterochromatic region on the plot.


### Step 2 - Generate chromosome plot

```bash
./plot_chromosome.py --help
```

```
usage: plot_chromosome.py [-h] -c CHROMOSOME [-w WIDTH] [-p PIXEL] [-s DARKSTART] [-m DARKMAX] -f FOLDER [--textFile]
                          [-g GENEFILE] [-t TEST] [--dummy] --config CONFIG [-l LOGFILE]

Script to plot genome chunks colored based on GC content and gene annotation.
See github: https://github.com/DSuveges/GenomePlotter

optional arguments:
  -h, --help            show this help message and exit
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Selected chromosome to process
  -w WIDTH, --width WIDTH
                        Number of chunks in one row.
  -p PIXEL, --pixel PIXEL
                        The size of a plotted chunk in pixels (default: 3).
  -s DARKSTART, --darkStart DARKSTART
                        Fraction of the width from where the colors start
                        getting darker (default: 0.75).
  -m DARKMAX, --darkMax DARKMAX
                        How dark a pixel can get at the right end of the plot
                        (default: 0.15).
  -f FOLDER, --folder FOLDER
                        Folder into which the plots are saved.
  --textFile            Flag to indicate if svg file should also be saved.
  -g GENEFILE, --geneFile GENEFILE
                        A .bed file with genes to add to the chromosome.
  -t TEST, --test TEST  The number of chunks to be read (by default the whole
                        chromosome is processed.)
  --dummy               If instead of the chunks, a dummy is drawn with
                        identical dimensions
  --config CONFIG       Specifying json file containing custom configuration
  -l LOGFILE, --logFile LOGFILE
                        File into which the logs are generated.
```

The script at first assigns GENCODE feature to each chunk as follows: the default value is intergenic, if a chunk has at least one base overlap with a gene then the chunk is considered to be gene, unless the chunk has at least one basepair overlap with an exon in which case the cunk is consideret to be exon, unless the GC content is NA, in which calse the chunk is considered to be heterochromatin. Based on cyto-band annotation, chunks overlapping with centromeres will be colored accordingly. The default color is adjusted based on the GC content.

Then genome-wide association signals are added as black dots. The size of the dots depends on the number of independent associations on a given chunk. Then cytological bands are added on the left side of the chromosome. Finally, if a gene set is given, the genes on the given chromosome are marked on the right side of the chromosome.

Finally the `.png` file is saved (and `.svg` file if required).

### Gene sets

The `gene_sets/` folder contains a set of files that can be used as gene annotation:

* *kinases_Hs.bed.gz*: list of kinase genes in the human genome, generated by `get_kinases.sh` script.
* *gene_w_drugs.tsv.gz*: list of genes for which approved drugs exists (where the mechanism of action is known) based on [OpenTargets tractability](https://docs.targetvalidation.org/getting-started/target-tractability) data.

### Result

The following image was created based on the data of chromosome 20, where 450 bp-s were averaged to get GC content, and 200 of these chunks were plotted in each row. The list of kinase genes are shown on the right side of the chromosome.

<img src="plots/chr20.png" alt="Chromosome 20" height="1500"/>


The combined plot showing the entire genome with all the protein kinases highlighted:

<img src="plots/kinases.png" alt="Protein kinases in the genome" width="1000">

### TO-DOs

1. Extra script to generate figure legend also as a csv and png.
2. Script to create the whole composition with all the chromosomes. Resolving the alignment is not trivial.
