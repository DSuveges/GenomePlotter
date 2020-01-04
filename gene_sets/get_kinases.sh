#!/usr/bin/env bash

# This script prepares input file for gene annotation
# This time all the kinases of the human genome will be added to the annotation.

# Steps:
# 1. Downloading list of kinases 
# 2. Parse kinase name and class
# 3. Based on the name the coordinates are fetched.
# 4. Prepare sorted bed file.


# Download kinases:
kinasesURL="http://kinase.com/web/current/kinbase/genes/SpeciesID/9606/fasta/protein/"

# Save kinase names and classes:
curl -s $kinasesURL | grep ">" | perl -lane '($name, $class) = $_ =~ /Hs_(.+?).AA class=(.+?):/; print "$name $class"' > kinase_names.lst

cat <(echo -e "chr\tstart\tend\tensemblID\tname\tclass") \
<(cat kinase_names.lst  | while read kinaseName class; do

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
done | sort -k1,1 -k2,2n ) | bgzip > kinases_Hs.bed.gz


# 