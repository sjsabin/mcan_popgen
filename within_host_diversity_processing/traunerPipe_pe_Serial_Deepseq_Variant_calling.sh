#!/bin/bash

# by Susanna Sabin
# Jan 2021
# Based on scripts provided by Trauner et al. on Github (https://github.com/SwissTPH/TBRU_serialTB/blob/master/scripts/variant_extraction/0_Serial_DeepSeq_Variant_calling.sh)
# This script will process BAM files previously generated and processed by traunerPipe_pe_bwa.sh.
#
# Paths and module loading are specified for use on ASU Agave HPC
#
# Usage:
#	./traunerPipe_pe_bwa.sh [INDIR] [FILE] [OUTDIR]
# Usage for slurm submission, where list.txt contains list of names of paired-end samples:
#	LIST=$(cat list.txt)
#	for i in $LIST; do
#		sbatch -n 4 -t 4-00:00 -o ${i}.out -e ${i}.err -J scy_sick_${i} --wrap="./traunerPipe_pe_Serial_Deepseq_Variant_calling.sh [INDIR] ${i} [OUTDIR]"
#	done

# Provide the date and time to a variable
ltime=`date +%Y-%m-%d-%H%M%S`

# Load relevant modules
module load samtools/1.9
module load perl/5.26.0
module load python/2.7.11
module load varscan/2.3.9
module load lofreq_star/2.1.3.1
module load r/3.2.3

INDIR=$1
FILE=$2
OUTDIR=$3

REF=/scratch/sjsabin/MTB_emp_dataProc/ref/NC_000962.3.fasta 
SCRIPTS=/scratch/sjsabin/traunerPipe_scripts

VARSCAN=/packages/7x/varscan/2.3.9/VarScan.v2.3.9.jar

# Start output report
echo "Trauner SNP Calling Pipeline: LoFreq and Varscan SNP calling and filtering"
echo "Log "$ltime 
echo ""
echo "Sample: "$FILE
echo "Input directory: "$INDIR
echo "Output directory: "$OUTDIR



echo -e "Processing $FILE\n"

samtools sort $INDIR/${FILE}.sort.rmdup.realn.bam -o $OUTDIR/${FILE}.sort.bam

echo -e "Generating mpileup for $FILE\n"


#This step generates a mpileup and filters out base_q<30, mapping_q<20, and outputs position data
#the pileup file is stored in temporary storage to be deleted at the end of the run.

samtools mpileup -q 30 -Q 20 -BOf $REF $OUTDIR/${FILE}.sort.bam > $OUTDIR/${FILE}.pileup
echo -e "mpileup for $FILE made\n"

#Use LoFreq to filter out SNPs that:
# don't comply with the default quality parameters of LoFreq.
# have a coverage of greater than 50, or less than 3000

echo -e "Calling variants in $FILE with LoFreq\n"
#call low frequency SNPs using lofreq
lofreq call -f $REF -o $OUTDIR/${FILE}.lofreq.snp $OUTDIR/${FILE}.sort.bam
#filter false positive SNPs with lofreq
lofreq filter --no-defaults --sb-mtc holm --cov-min 50 --cov-max 3000 -i $OUTDIR/${FILE}.lofreq.snp -o $OUTDIR/${FILE}.lofreq.filt.snp
echo -e "LoFreq finished\n\n"

#Use VarScan to filter out SNPs that have:
# fewer than 4 reads supporting them,
# a coverage of less than 50,
# min average quality <20,
# min frequency <0.5%
# max frequency <90%

echo -e "Calling variants in $FILE with VARSCAN\n"
java -jar $VARSCAN mpileup2snp $OUTDIR/${FILE}.pileup --min-coverage 50 --min-reads2 4 --min-avg-qual 20 --min-var-freq 0.005 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > $OUTDIR/${FILE}.varscan.var
echo -e "VARSCAN finished\n\n"

#Merge the LoFreq and Varscan output

echo -e "Merge congruent variants in $FILE\n"
perl $SCRIPTS/varscan_lofreq_compare.pl $OUTDIR/${FILE}.lofreq.filt.snp $OUTDIR/${FILE}.varscan.var > $OUTDIR/${FILE}.var
echo -e "Merge done\n\n"

#Re-format the file for downstream processing

echo -e "Re-formatting the file\n"
perl  $SCRIPTS/1_format_trans.pl $OUTDIR/${FILE}.var > $OUTDIR/${FILE}.for
echo -e "Re-formatting done\n\n"

#Get heterozygous SNP candidates:

echo -e "Extracting low frequency SNP candidates\n"
perl $SCRIPTS/3_mix_extract.pl $OUTDIR/${FILE}.for > $OUTDIR/${FILE}.lofreq
echo -e "Extraction done\n\n"

#Merge info to produce the data to use for a Kolmogorov-Smirnov test of distributions:

echo -e "Merging the pileups and low freq annotations\n"
perl $SCRIPTS/3.1_mix_pileup_merge.pl $OUTDIR/${FILE}.pileup  $OUTDIR/${FILE}.lofreq >  $OUTDIR/${FILE}.merge
echo -e "Merging done\n\n"

#R script to remove SNPs whose support comes from the tail-end of the distribution.

echo -e "Running R: Kolmogorov–Smirnov test of SNP support\n"
Rscript $SCRIPTS/tp_fp_test_new_151_201503.R $OUTDIR/${FILE}.merge
echo -e "KS test done\n\n"

#Python script to remove SNPs at specific sites, if frequency is less than 0.5% and the KS
#test is positive (p<0.05).

echo -e "Running Python: filter out SNPs from problematic regions, KS+ SNPs\n"
python $SCRIPTS/5_lable_filter.py $OUTDIR/${FILE}.merge.lable $SCRIPTS/4_H37Rv_PPE_INSERT_loci_list
echo -e "$FILE sample fully processed\n\n"

echo -e "Annotating filtered SNP calls\n"
#annotate SNPs using H37Rv template(NC_000962.2)
perl $SCRIPTS/6_H37Rv_annotate.pl $SCRIPTS/1_Tuberculist.info $SCRIPTS/2_genetic_codes $REF $OUTDIR/${FILE}.filter.snp > $OUTDIR/${FILE}.ano.snp
echo -e "Annotation done\n"
