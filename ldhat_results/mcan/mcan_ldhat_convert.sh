#!/bin/bash

# LDHat "convert" run on empirical data
# Usage
	# script.sh [split] [element]

# Extract and build LDHat software locally
tar xf LDhat-master.tar.gz
cd LDhat-master
make
cd ../
LDhat-master/convert -seq full.multi.fasta
