## GenomePlotter

The motivation behind this project was to visualize how the genes and already known genome wide association signals are 
distributed across the chromosomes. I chose a heat-map kind of visualization that allows the representation of the underlying chemical
properties of the DNA (GC content). 

### Approach

The project is divided into two parts: 

1. A shell script that downloads and pre-processes the source data files with genome sequence, GWAS signals, gene annotation. 
2. A python workbook, that combines all source data into a single data-frame and creates the final plots.

### Source data:

All applies source data is mapped to the GRCh38 build of the human genome.

* **Genome sequence** dowloaded from [Ensembl](ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/) version 84. 
* **Genome wide association signals** most recent version of the NHGRI-EBI [GWAS catalog](https://www.ebi.ac.uk/gwas/) was downloaded at 2016.04.10.
* **Gene annotation** gene coordinates were downloaded from [GENCODE](ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/) version 24. 

### Step 1 - Pre-processing.

#### GWAS signals

A bedfile is created with the 
