#!/usr/bin/env bash

echo "This is a temporary solution to fetch input data for the drugs module."
echo "If no access to GCP storeage, files can be fetched from EBI FTP: ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/"

TARGET_DIR='source_data'
OT_RELEASE='22.11'

# Fetch evidence:
gsutil cp -r gs://open-targets-data-releases/${OT_RELEASE}/input/evidence-files/genetics-portal-evidences.json.gz ${TARGET_DIR}/

# Fetch targets:
gsutil cp -r gs://open-targets-data-releases/${OT_RELEASE}/output/etl/parquet/targets ${TARGET_DIR}/

# Fetch molecules:
gsutil cp -r gs://open-targets-data-releases/${OT_RELEASE}/output/etl/parquet/molecule ${TARGET_DIR}/

# Fetch disease:
gsutil cp -r gs://open-targets-data-releases/${OT_RELEASE}/output/etl/parquet/diseases ${TARGET_DIR}/

