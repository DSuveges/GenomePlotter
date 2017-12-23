#!/usr/bin/env bash

# Description:
### This script was written to download and pre-process all necessary data
### for the genome plotter script.

# Method:
# 1. Downloading source data: human genome from Ensembl, gene definitions from GENCODE, signals from the GWAS catalog
# 2. Process downloaded data: format and calculate GC content for the whole genome.
# 3. Saving all outputs required by the plotter.

# v.1.2 Last modified: 2017.12.22
    # Before downloading, check the most recent GENCODE and Ensembl version.
    # Also output reports about the available verions.
    # Heuristic approach to find relevant columns in the GWAS file.
    # Both exons and genes are handled. There will be three colors on the plot.
    # Bedtools are required and tabix.... I'm not sure if we need tabix for this.
    # All the resulting files are indexed by tabix.

# Checking dependencies:
if [[ -z $(which tabix ) ]]; then "[Error] Tabix is required. Exiting"; exit; fi
if [[ -z $(which mergeBed ) ]]; then "[Error] Bedtools is required. Exiting"; exit; fi

# Setting working dir:
workingDir=$(pwd)

# Setting up source directories:
GWAS_file="www.ebi.ac.uk/gwas/api/search/downloads/full"
Ensembl_base_URL="ftp://ftp.ensembl.org/pub"
GENCODE_base_URL="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human"

# Determine the most recent GENCODE release:
GENCODE_release=$(curl -s --list-only ${GENCODE_base_URL}/ | grep -i release | sed -e 's/release_//' | sort -n | tail -n1)
if [[ -z ${GENCODE_release} ]]; then 
    echo "[Warning] Failed to determine the most recent GENCODE version. Using version 27 (released 22/08/2017)."
    GENCODE_release=24
fi

# Determine the most recent Ensembl release:
ENSEMBL_release=$(curl -s --list-only ${Ensembl_base_URL}/ | grep -i release | sed -e 's/release-//' | sort -n | tail -n1)
if [[ -z ${ENSEMBL_release} ]]; then 
    echo "[Warning] Failed to determine the most recent Ensembl release. Using version 91 (released 11/12/2017)."
    ENSEMBL_release=91
fi

# Creating data folder:
mkdir -p ${workingDir}/source_data # Downloaded sources prior to pre-processing.
mkdir -p ${workingDir}/data # Processed source data ready to use by the plotting script
mkdir -p ${workingDir}/processed_data # Saved dataframes ready to print, used when only minor changes are about to made

# Status update:
echo "[Info] Source files are generated for genome plotter. Date: " $(date "+%Y-%m-%d")
echo "[Info] All genomic coordinates are mapped to GRCh38 build."
echo "[Info] All files will be saved in a sorted bed format. More info: https://genome.ucsc.edu/FAQ/FAQformat.html#format1"

# Download GWAS file:
echo "[Info] Downloading the most recent NHGRI-EBI GWAS catalog."
wget -q ${GWAS_file} -O ${workingDir}/source_data/GWAS_catalog.tsv

# Checking download:
if [[ $? -ne 0 ]]; then
    echo "[Error] Failed to download the GWAS catalog. Exiting."
    exit
fi

# Extracting columns of interest:
columns=$(head -n1 ${workingDir}/source_data/GWAS_catalog.tsv | tr " " "_" | awk 'BEGIN{FS="\t"}{ for(i = 1; i <= NF; i++) { print i, $i; } }' | grep -w -e CHR_ID -e CHR_POS -e DISEASE/TRAIT -e SNPS  | cut -d" " -f1 )

# If the extraction has been failed, we exit from the script.
if [[ -z ${columns} ]]; then
    echo "[Error] The relevant fields from the GWAS file could not be extracted. Exiting."
    exit;
fi

# Extracting columns and properly format:
echo "[Info] Processing GWAS file: creating bed file and keeping only rsID and trait."
cat <(echo -e "#chr\tstart\tend\trsID\ttrait") <(cat ${workingDir}/source_data/GWAS_catalog.tsv | cut -f$(echo ${columns} | sed -e 's/ /,/g') | awk 'BEGIN{FS=OFS="\t"} NR != 1 && $3 ~ /^[0-9]+$/ { printf "%s\t%s\t%s\t%s\t%s\n", $2, $3 - 1, $3, $4, $1 }' | sort -k1,1 -k2,2n) | bgzip > ${workingDir}/data/processed_GWAS.bed.gz
tabix -p bed ${workingDir}/data/processed_GWAS.bed.gz

# Checking if the process was successful:
if [[ ! -e ${workingDir}/data/processed_GWAS.bed.gz ]]; then
    echo "[Error] Creation of the processed gwas file was failed. Exiting."
    exit;
else
    hitcount=$(zcat ${workingDir}/data/processed_GWAS.bed.gz | grep -v "#" | wc -l )
    snpcount=$(zcat ${workingDir}/data/processed_GWAS.bed.gz | grep -v "#" | cut -f4 | sort -u | wc -l )
    echo "[Info] Number of processed associations in the GWAS file: ${hitcount}, number of snps: ${snpcount}"
fi;

# Downloading GENCODE file and process
echo "[Info] Downloading GENCODE file (v. ${GENCODE_release})."
wget -q ${GENCODE_base_URL}/release_${GENCODE_release}/gencode.v${GENCODE_release}.annotation.gtf.gz -O ${workingDir}/source_data/gencode.v${GENCODE_release}.annotation.gtf.gz

# Checking download:
if [[ $? -ne 0 ]]; then
    echo "[Error] Failed to download the GENCODE file. Exiting"
    exit
fi

echo "[Info] Processing gencode data: extracting gene and exon positions. Then creating bedfile with the merged coordinates."
# GENCODE data is split into two parts: exons and genes.
zcat ${workingDir}/source_data/gencode.v${GENCODE_release}.annotation.gtf.gz | perl -F"\t" -lane '
    next unless $F[2] eq "gene" and $F[8] =~ /gene_type "protein_coding"/;
    $F[0] =~ s/chr//i;
    $F[8] =~ /gene_id "(ENSG.+?)";.+gene_name "(.+?)";/;
    print join("\t", ($F[0], $F[3], $F[4], $1, $2))' | sort -k1,1 -k2,2n | gzip > ${workingDir}/data/genes_GENCODE.bed.gz

mergeBed -i ${workingDir}/data/genes_GENCODE.bed.gz | awk 'BEGIN{FS=OFS="\t"}{print $0, "gene"}' | gzip > ${workingDir}/data/genes_GENCODE.merged.bed.gz

zcat ${workingDir}/source_data/gencode.v${GENCODE_release}.annotation.gtf.gz | perl -F"\t" -lane '
    next unless $F[2] eq "exon" and $F[8] =~ /gene_type "protein_coding"/;
    $F[0] =~ s/chr//i;
    $F[8] =~ /gene_id "(ENSG.+?)";.+gene_name "(.+?)";/;
    print join("\t", ($F[0], $F[3], $F[4], $1, $2))' | sort -k1,1 -k2,2n | gzip > ${workingDir}/data/exons_GENCODE.bed.gz

mergeBed -i  ${workingDir}/data/exons_GENCODE.bed.gz | awk 'BEGIN{FS=OFS="\t"}{print $0, "exon"}' | gzip > ${workingDir}/data/exons_GENCODE.merged.bed.gz

# Combining the exon and the gene datasets together:
zcat ${workingDir}/data/exons_GENCODE.merged.bed.gz ${workingDir}/data/genes_GENCODE.merged.bed.gz | sort -k1,1 -k2,2n | bgzip > ${workingDir}/data/GENCODE.merged.bed.gz
tabix -p bed ${workingDir}/data/GENCODE.merged.bed.gz

# Download and process Ensembl genomic data:
echo "[Info] Downloading the non-masked sequence of human genome (Ensembl release: ${ENSEMBL_release}). It may take a while..."
for chr in {1..22} X Y ; do
    echo -e "\tDownloading chromosome ${chr}"
    wget -q ${Ensembl_base_URL}/release-${ENSEMBL_release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz -O ${workingDir}/source_data/human_genome_GRCh38_chr${chr}.fa.gz

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
for chr in {1..22} X Y ; do
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
            | sort -k1,1 -k2,2n | bgzip > ${workingDir}/data/Processed_chr${chr}.bed.gz
    tabix -p bed ${workingDir}/data/Processed_chr${chr}.bed.gz # Indexing.
    
    chunk_no=$( gzcat ${workingDir}/data/Processed_chr${chr}.bed.gz | wc -l )
    echo -e "\tChromosome ${chr} is split into ${chunk_no} chunks."
done

echo "[Info] All source files have been downloaded and pre-processed for the plots."



