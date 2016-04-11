#!/usr/bin/bash

# v.1.0 Last modified: 2016.04.09

# Description:
### This script was written to download and pre-process all necessary data
### for the genome plotter script.

# Method:
### Downloads human genome
### Downloads the GENCODE gene database
### Downloads the NHGRI-EBI GWAS catalog

# Setting working dir:
workingDir=$(pwd)

# Setting up source directories:
GWAS_file="www.ebi.ac.uk/gwas/api/search/downloads/full"
GENCODE_file="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
Ensembl_folder="ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/"

# Creating data folder:
mkdir -p ${workingDir}/source_data # Downloaded sources prior to pre-processing.
mkdir -p ${workingDir}/data # Processed source data ready to use by the plotting script
mkdir -p ${workingDir}/processed_data # Saved dataframes ready to print, used when only minor changes are about to made

# Status update:
echo "[INFO] Source files are generated for genome plotter. Date: " $(date "+%Y-%m-%d")
echo "[INFO] All genomic coordinates are mapped to GRCh38 build."
echo "[INFO] All files will be saved in a sorted bed format. More info: https://genome.ucsc.edu/FAQ/FAQformat.html#format1"

# Download GWAS file:
echo "[Info] Downloading the most recent NHGRI-EBI GWAS catalog."
wget -q ${GWAS_file} -O ${workingDir}/source_data/GWAS_catalog.tsv

# Checking download:
if [[ $? -ne 0 ]]; then
    echo "[Error] Failed to download the GWAS catalog. Exiting"
    exit
fi

# Creating proper sorted bedfile.
# rsID = 30
# chr = 12
# pos = 13
# trait = 8
echo "[Info] Processing GWAS file: creating bed file and keeping only rsID and trait."
cat ${workingDir}/source_data/GWAS_catalog.tsv | awk 'BEGIN{FS=OFS="\t"}{if (NR != 1 && $13 ~ /^[0-9]+$/ ){print $12, $13 - 1, $13, $30, $8}}' \
    | sort -k1,1 -k2,2n | gzip > ${workingDir}/data/processed_GWAS.bed.gz


# Downloading GENCODE file and process
echo "[Info] Downloading version 24. GENCODE file."
wget -q ${GENCODE_file} -O ${workingDir}/source_data/gencode.v24.annotation.gtf.gz

# Checking download:
if [[ $? -ne 0 ]]; then
    echo "[Error] Failed to download the GENCODE file. Exiting"
    exit
fi

# Processing gencode data:
echo "[Info] Processing gencode data: creating bedfile."
zcat source_data/gencode.v24.annotation.gtf.gz | perl -F"\t" -lane '
    next unless $F[2] eq "gene" and $F[8] =~ /gene_type "protein_coding"/;
    $F[0] =~ s/chr//i;
    $F[8] =~ /gene_id "(ENSG.+?)";.+gene_name "(.+?)";/;
    print join("\t", ($F[0], $F[3], $F[4], $1, $2))' | sort -k1,1 -k2,2n | gzip > ${workingDir}/data/processed_GENCODE.bed.gz

# Download and process Ensembl genomic data:
echo "[Info] Downloading the non-masked sequence of human genome. It may take a while..."
for chr in {1..22}; do
    echo -e "\tDownloading chromosome ${chr}"
    wget -q ${Ensembl_folder}/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz -O ${workingDir}/source_data/human_genome_GRCh38_chr${chr}.fa.gz

    # Checking download:
    if [[ $? -ne 0 ]]; then
        echo "[Error] Failed to download the chromosome file. Exiting"
        exit
    fi
done

# The threshold for chunking the genome:
export ChunkSize=450

# Chunking the human genome, calculating GC content.
echo "[Info] Processing chrmosomes. Chunks size: ${ChunkSize}bp. Sorted bed file will be saved."
# Downloading all chromosomes:
for chr in {1..22}; do
    export chr
    # Processing downloaded data:
    gzcat ${workingDir}/source_data/human_genome_GRCh38_chr${chr}.fa.gz \
            | perl -lane '
                $threshold = $ENV{ChunkSize}; # Importing chunk size
                $chr = $ENV{chr}; # Importing chromosome name

                die "[Error] No chromosome number exported! Exiting." unless $chr;
                die "[Error] No threshold has been set and exported! Exiting." unless $threshold;

                next if $_ =~ />/; # Excluding header line

                chomp $_; # Remove whitespace
                $chunkEnd += length($_); # Counting position throughout the entire datafile.

                $chunkStart = $chunkEnd - length( $_ ) unless $chunk; # Updated when starting new chunk.

                my $N_count = () = $_ =~ /N/gi;
                next if $N_count > length($_)/ 2; # Excluding lines with too many N-s

                $chunk .= $_; # Adding current line to chunk.

                if ( length($chunk) >= $threshold ){
                    $chunkCount += 1;

                    my $CG_count = () = $chunk =~ /[CG]/gi;
                    my $AT_count = () = $chunk =~ /[AT]/gi;

                    my $CG_content = $CG_count / ($CG_count + $AT_count);

                    printf "%s\t%s\t%s\t%s\t%s\n", $chr, $chunkStart, $chunkEnd, $CG_content, $chunkCount;

                    $chunk = ""; # Emptying chunk.
                }' \
            | gzip > ${workingDir}/data/Processed_chr${chr}.bed.gz

    chunk_no=$( gzcat ${workingDir}/data/Processed_chr${chr}.bed.gz | wc -l )
    echo -e "\tChromosome ${chr} is split into ${chunk_no} chunks."
done

echo "[Info] All source files have been downloaded and pre-processed for the plots."



