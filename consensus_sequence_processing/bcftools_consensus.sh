#!/bin/bash
#SBATCH -n 1
#SBATCH -J bcftools
#SBATCH -o bcftools.${1}.out
#SBATCH -e bcftools.${1}.err

module purge
module load bcftools/1.9
bcftools consensus -H 1 -f /home/sjsabin/mcan/ref/GCA_000253375.1_ASM25337v1_genomic.fasta -s $1 -o ${1}_consensus.fa /scratch/sjsabin/freebayes/${1}.subsamp25.sitesExcl.snps.recode.filt.name.vcf.gz
