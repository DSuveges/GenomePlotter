## GenomePlotter

The motivation behind this project was to visualize how the genes and already known genome wide association signals are 
distributed across the chromosomes. I chose a heat-map kind of visualization that allows the representation of the underlying chemical
properties of the DNA (GC content). 

### Approach

The project is divided into two parts: 

1. A shell script that downloads and pre-processes the source data files with genome sequence, GWAS signals, gene annotation. 
2. A python workbook, that combines all source data into a single data-frame and creates the final plots.

### Requirements

**Required bash tools:**
* [cairo graphics library](https://www.cairographics.org/download/)
* [tabix](http://www.htslib.org/download/)
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)

**Required Python packages:**

* [pandas](https://pandas.pydata.org/)
* [numpy](http://www.numpy.org/)
* [cairosvg](http://cairosvg.org/)
* [pybedtools](https://pypi.python.org/pypi/pybedtools)

### Source data:

All applies source data is mapped to the GRCh38 build of the human genome.

* **The sequnce of the human genome** was dowloaded from [Ensembl](http://www.ensembl.org/info/data/ftp/index.html) checked for the most recent version. 
* **Genome wide association signals** most recent version of the NHGRI-EBI [GWAS catalog](https://www.ebi.ac.uk/gwas/) downloading the most recent version.
* **Gene annotation** gene coordinates were downloaded from [GENCODE](http://www.gencodegenes.org/releases/current.html) checking for the most recent version. 

### Step 1 - Pre-processing.

#### GWAS signals

Data is download from the GWAS catalog following the fixed link that is being updated regularly. Only those signals are kept that properly mapped with genomic coordinates.An indexed bedfile is created for further processing.

#### GENCODE annotation

Gencode website is checked for the most recent version. Then the compressed gtf file is downloaded. The coordinates of protein coding genes and the corresponding exons are extracted then merged using [mergeBed](http://bedtools.readthedocs.io/en/latest/content/tools/merge.html). Exons and genes are pooled together and a sorted, indexed bedfile is generated.

#### Sequences


### Step 2 - Plotting

#### Saving processed data

#### Results:

The plot of the 11th chromosomes of the human genome:
![chr11 - genes, GWAS signals and CG content](https://github.com/DSuveges/GenomePlotter/blob/master/plot_chr11.png)

