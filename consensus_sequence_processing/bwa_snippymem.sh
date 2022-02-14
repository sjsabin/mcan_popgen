#!/bin/bash

LIST=$(cat /home/sjsabin/mcan/data/trim_fq/list.txt)

refdir=/home/sjsabin/mcan/ref

for FILE in $LIST; do
	# Load bwa and samtools modules for agave.
	module load bwa/0.7.17
	module load samtools/1.9	
	# Make output directory for each sample.
	mkdir -p /scratch/sjsabin/bwa_mem/${FILE}
	dir=/scratch/sjsabin/bwa_mem/${FILE}
	# Extract zipped fastq.gz files for alignment.
	gunzip /home/sjsabin/mcan/data/trim_fq/$FILE/$FILE_*.trim.fastq.gz
	# Align samples to the MCAN reference using Snippy settings. Pipe alignment to be sorted, then go through mate fixing, then sorted again. Save this file prior to removing duplicates.
	bwa mem -Y -M -R @RG\tID:ILLUMINA-${FILE}\tSM:${FILE}\tPL:illumina ${refdir}/GCA_000253375.1_ASM25337v1_genomic.fasta /home/sjsabin/mcan/data/trim_fq/$FILE/$FILE_1.trim.fastq /home/sjsabin/mcan/data/trim_fq/$FILE/$FILE_2.trim.fastq | samtools sort -n -l 0 | samtools fixmate -m | samtools sort -O bam -o ${dir}/${FILE}-pe.sorted.bam
	echo "Finished bwa mem pipe."
	# Remove duplicates (-r option)!
	samtools markdup -r -s ${dir}/${FILE}-pe.sorted.bam ${dir}/${FILE}-pe.sorted.dedup.bam
	echo "Finished samtools markdup."
	# Index the bam file that will be used for variant calling.
	samtools index ${dir}/${FILE}-pe.sorted.dedup.bam
	echo "Finished samtools index."
	# Rezip the fastq files for keeping everything tidy...
	gzip /home/sjsabin/mcan/data/trim_fq/$FILE/$FILE_*.trim.fastq

done
