#!/bin/bash

# Execution script for running LDhat pairwise on numerous simulation models.
# Usage:
#	./er_psi_bn_dfe_pw.sh [0-99]

# Load stashcache
#module load stashcache

# Extract and build LDHat software locally
tar xf LDhat-master.tar.gz
cd LDhat-master
make
cd ../
# Copy archive of fastas from the public directory to the host directory
#stashcp /osgconnect/public/sjsabin/er_psi_bn_dfe_convert_out.tar.gz ./
# Extract the fasta files
tar xf er_psi_bn_dfe_convert_out.tar.gz

REP=er_psi_bn_dfe.$1

echo "0.000539 100 201 0 0 1 1 1" | LDhat-master/pairwise -seq er_psi_bn_dfe_convert_out/${REP}.sites.txt -loc er_psi_bn_dfe_convert_out/${REP}.locs.txt -prefix ${REP}. >> ${REP}.pw.log

# Remove the large input file from the host directory
rm -r er_psi_bn_dfe_convert_out/
rm -r LDhat-master/