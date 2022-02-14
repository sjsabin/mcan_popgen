#!/bin/bash
#!/bin/perl

# by Susanna Sabin
# Jan 2021
# Based on scripts provided by Trauner et al. on Github (https://github.com/SwissTPH/TBRU_serialTB/blob/master/scripts/variant_extraction/support_scripts/0_sanger_triming_paired_ended.sh)
# This script will process raw PE fastq files with the programs scythe and sickle, which will perform adapter trimming and quality filtering.
# Paths and module loading are specified for use on ASU Agave HPC

# Usage:
#	./traunerPipe_pe_scythe_sickle.sh [INDIR] [FILE] [OUTDIR]
# Usage for slurm submission, where pe_list.txt contains list of names of paired-end samples:
#	LIST=$(cat pe_list.txt)
#	for i in $LIST; do
#		sbatch -n 4 -t 4-00:00 -o ${i}.out -e ${i}.err -J scy_sick_${i} --wrap="./traunerPipe_pe_scythe_sickle.sh [INDIR] ${i} [OUTDIR]"
#	done

# Provide the date and time to a variable
ltime=`date +%Y-%m-%d-%H%M%S`

# Assign the arguments to variables
INDIR=$1
FILE=$2
OUTDIR=$3

# Load scythe and sickle modules
module load scythe/0.994
module load sickle/1.33

#
echo "Trauner SNP Calling Pipeline: Scythe and Sickle fastq processing"
echo "Log "$ltime 
echo ""
echo "Sample: "$FILE
echo "Input directory: "$INDIR
echo "Output directory: "$OUTDIR

# Based on scripts provided by Trauner et al. on Github (https://github.com/SwissTPH/TBRU_serialTB/blob/master/scripts/variant_extraction/support_scripts/0_sanger_triming_paired_ended.sh)

# Define the path to a TruSeq adapter file on Agave.
adapters=/packages/7x/scythe/0.994/illumina_adapters.fa
# Unzip the gzipped fastq file(s)
gunzip $INDIR/${FILE}*.gz


# Run scythe on paired-end, unzipped fastq files
#	-a : path to adapter fasta file
#	--quiet : 
#	-q : type of quality encoding, either "sanger" or "solexa"
#	-o : path to output file
echo "Beginning scythe..."
echo ""
scythe -a $adapters -q sanger -o $OUTDIR/${FILE}_1.scy.fastq $INDIR/${FILE}_1.fastq
scythe -a $adapters -q sanger -o $OUTDIR/${FILE}_2.scy.fastq $INDIR/${FILE}_2.fastq
echo ""
echo "Scythe complete."

# Run sickle on paired end fastq files (scythe output)
#	pe	: either "pe" or "se" to denote paired-end or single-end read files respectively
#	-f	: forward input fastq file (*_1)
#	-r	: reverse input fastq file (*_2)
#	-t	: quality type
#	-o	: forward output fastq file (*_2)
#	-p	: reverse output fastq file (*_1)
#	-s	: "singles" output file (reads that passed quality filtering in one direction, but not the other)
echo "Beginning sickle..."
sickle pe -f $OUTDIR/${FILE}_1.scy.fastq -r $OUTDIR/${FILE}_2.scy.fastq -t sanger -o $OUTDIR/${FILE}_1.scy.sick.fastq -p $OUTDIR/${FILE}_2.scy.sick.fastq -s $OUTDIR/${FILE}_singles.scy.sick.fastq
echo "Sickle complete."

echo "Beginning cleanup..."
# Clean up intermediate files
rm $OUTDIR/${FILE}*scy.fastq

# Zip terminal files
gzip $OUTDIR/${FILE}*.scy.sick.fastq
echo "Cleanup complete."
