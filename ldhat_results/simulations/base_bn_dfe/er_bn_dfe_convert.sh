#!/bin/bash

# Execution script for running LDhat convert on numerous simulation models.
# Usage:
#	./er_bn_dfe_convert.sh [modified_fasta_file]

# Load stashcache
module load stashcache

# Extract and build LDHat software locally
tar xf LDhat-master.tar.gz
cd LDhat-master
make
cd ../
# Copy archive of fastas from the public directory to the host directory
stashcp /osgconnect/public/sjsabin/er_bn_dfe_sub_slim.tar.gz ./
# Extract the archive of fastas
tar xvf er_bn_dfe_sub_slim.tar.gz
# Chop up the input fasta file name into different variables
REP=$(echo "$1" | cut -f1-2 -d".")
FILE=$(echo "$1" | cut -f1-5 -d".")
# Extract the fasta file
gunzip er_bn_dfe_sub_slim/$1
# Run LDhat convert on the modified fasta file
LDhat-master/convert -seq er_bn_dfe_sub_slim/${FILE} -prefix ${REP}. > ${REP}.convert.log
# Remove the large input file from the host directory
rm -r er_bn_dfe_sub_slim/
