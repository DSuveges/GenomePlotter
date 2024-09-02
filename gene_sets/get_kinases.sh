#!/usr/bin/env bash

# This script prepares input file for gene annotation
# This time all the kinases of the human genome will be added to the annotation.

# Steps:
# 1. Downloading list of kinases
# 2. Parse kinase name and class
# 3. Based on the name the coordinates are fetched.
# 4. Prepare sorted bed file.

# Prepare environment:
scriptDir=$(pwd)
mkdir -p "${scriptDir}/gene_sets"

##
## Get all protein kinases from the kinase.com website:
##

# Download kinases:
kinasesURL="http://kinase.com/web/current/kinbase/genes/SpeciesID/9606/fasta/protein/"

# Save kinase names and classes:
curl -s $kinasesURL | grep ">" | perl -lane '($name, $class) = $_ =~ /Hs_(.+?).AA class=(.+?):/; print "$name $class"' > "${scriptDir}/gene_sets/prot_kinase_names.lst"

cat <(echo -e "chr\tstart\tend\tensemblID\tname\tclass") \
<(cat "${scriptDir}/gene_sets/prot_kinase_names.lst"  | while read kinaseName class; do

	# Fetch Ensembl ID based on kinase name:
	ensemblID=$(curl -s "http://rest.ensembl.org/xrefs/symbol/homo_sapiens/${kinaseName}?content-type=application/json" --max-time 3 --retry-delay 2 --retry 5 --user-agent "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.79 Safari/537.36" | jq -r '.[0].id')
	if [[ -z ${ensemblID} ]]; then
		echo "[Warning] Ensembl ID for ${kinaseName} was not found. Skipping." >&2
		continue
	fi

	# Fetch genomic location based on Ensembl ID:
	read chr start end <<<$(curl -s "http://rest.ensembl.org/lookup/id/${ensemblID}?content-type=application/json;expand=0" --max-time 3 --retry-delay 2 --retry 5 --user-agent "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.79 Safari/537.36" | jq -r '"\(.seq_region_name) \(.start) \(.end)"')
	if [[ -z ${start} ]]; then
		echo "[Warning] Genomic location could not be retrieved for ${ensemblID} ($kinaseName). Skipping" >&2
		continue
	fi

	# Save bedfile:
	echo -e "${chr:-NA}\t${start:-NA}\t${end:-NA}\t${ensemblID:-}\t${kinaseName}\t${class}"
	echo -ne "." >&2

	# Wait a bit:
	sleep 1
done | sort -k1,1 -k2,2n ) | bgzip > "${scriptDir}/gene_sets/prot_kinase.bed.gz"


##
## Get all kinases (not only protein kinases) from the Uniprot:
##

# Retrieve all kinases:
curl -s \
    --data-urlencode 'query=reviewed:yes organism:"Homo sapiens (Human) [9606]" goa:("kinase activity [16301]")' \
    --data-urlencode 'columns=id,entry name,protein names,genes,go,go(molecular function)' \
    --data-urlencode 'format=tab' \
    'https://www.uniprot.org/uniprot/' > "${scriptDir}/gene_sets/all_kinases.tsv"

# Extract gene names, look up on Ensembl, fetch coordinates for the genes:
cat <(echo -e "chr\tstart\tend\tensemblID\tname\tuniprot_entry") \
<(tail -n+2 "${scriptDir}/gene_sets/all_kinases.tsv" | cut -f1,4 | while read uniprotID names; do
	read kinaseName <<<$(echo $names | cut -f1 -d " " )

	# Fetch Ensembl ID based on kinase name:
	read IDs <<<$( curl -s "http://rest.ensembl.org/xrefs/symbol/homo_sapiens/${kinaseName}?content-type=application/json" --max-time 3 --retry-delay 2 --retry 5 --user-agent "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.79 Safari/537.36" | jq -r '.[].id'| tr '\n' ' ' )
	if [[ -z ${IDs} ]]; then
		echo "[Warning] Ensembl ID for ${kinaseName} was not found. Skipping." >&2
		continue
	fi

	# Empty
	ensemblID='';

	# Loop through all the fetched IDs and check if the assigned HGNC name is indeed the same as the gene name:
	for id in $IDs; do

		# Fetch genomic location based on Ensembl ID:
		read chr start end display_name <<<$(curl -s "http://rest.ensembl.org/lookup/id/${id}?content-type=application/json;expand=0" --max-time 3 --retry-delay 2 --retry 5 --user-agent "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.79 Safari/537.36" | jq -r '"\(.seq_region_name) \(.start) \(.end) \(.display_name)"' )
		if [[ -z ${start} ]]; then
			echo "[Warning] Genomic location could not be retrieved for ${id} ($kinaseName). Skipping" >&2
			continue
		fi

		# Test if the corresponding name is the same:
		if [[ $display_name == $kinaseName ]]; then
			ensemblID=$id;
			# Save bedfile:
			echo -e "${chr:-NA}\t${start:-NA}\t${end:-NA}\t${ensemblID:-}\t${kinaseName}\t${uniprotID}"
			echo -ne "." >&2

			# Wait a bit (help with the rest queries):
			sleep 1
		else
			echo -ne 'x' >&2
		fi
	done
done | sort -k1,1 -k2,2n ) | bgzip > "${scriptDir}/gene_sets/all_kinases.bed.gz"


##
## Run annotator script to add GWAS signals + cytobands + all kinases:
##

# Re-run all annotation:
for chr in {1..22} X Y; do
	if [[ -e "${scriptDir}/plots/chr${chr}_annotated.png" ]]; then
		echo "[Info] Skipping $chr";
		continue;
	fi
	python "${scriptDir}/chromosome_annotator.py" --chromosome ${chr} --geneFile "${scriptDir}/gene_sets/all_kinases.bed.gz"
	sleep 120
done
