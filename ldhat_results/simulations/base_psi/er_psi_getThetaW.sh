#!/bin/bash

BASES=413587
THETA=$(grep "Watterson" er_psi.${1}.cons.convert.log | awk '{split($0,a," ") ; print a[4]}')
THETA_pb=$(awk -v NUM=${THETA} -v BASES=${BASES} 'BEGIN { print (NUM/BASES)}')
echo -e $1"\t"$THETA_pb >> er_psi.consensus_pairwise_args.txt
