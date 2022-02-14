#!/bin/bash

# Based on scripts provided by Trauner et al. on Github (https://github.com/SwissTPH/TBRU_serialTB/blob/master/scripts/variant_extraction/support_scripts/1_Mapping_bwa_paired_ended.sh)

# Usage:
#	./traunerPipe_pw_bwa.sh [INDIR] [FILE] [OUTDIR]
# Usage for slurm submission, where se_list.txt contains list of names of paired-end samples:
#	LIST=$(cat se_list.txt)
#	for i in $LIST; do
#		sbatch -n 4 -t 4-00:00 -o ${i}.out -e ${i}.err -J ${i}_bwa --wrap="./traunerPipe_pe_bwa.sh [INDIR] ${i} [OUTDIR]"
#	done

# Provide the date and time to a variable
ltime=`date +%Y-%m-%d-%H%M%S`

# Load relevant modules
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.9.2
module load gatk/3.7.0
module load bcftools/1.10.2
module load perl/5.26.0

# Assign the arguments to variables
INDIR=$1
FILE=$2
OUTDIR=$3

# Write the reference genome directory to a variable
REFDIR=/scratch/sjsabin/MTB_emp_dataProc/ref

# Define variables to represent intermediate and output files


echo "Trauner SNP Calling Pipeline: BWA alignment"
echo "Log "$ltime 
echo ""
echo "Sample: "$FILE
echo "Input directory: "$INDIR
echo "Output directory: "$OUTDIR

# Unzip the input files
gunzip $INDIR/${FILE}_1.scy.sick.fastq.gz
gunzip $INDIR/${FILE}_2.scy.sick.fastq.gz
gunzip $INDIR/${FILE}_singles.scy.sick.fastq.gz

#perform bwa paired end mapping
echo -e "\nBwa mapping started ...\n"
echo -e "1. Initial mapping of ${FILE} by bwa ... \n"

bwa aln -R 1 $REFDIR/NC_000962.3.fasta $INDIR/${FILE}_1.scy.sick.fastq > $OUTDIR/${FILE}_1.sai
echo "bwa aln #1 done"
bwa aln -R 1 $REFDIR/NC_000962.3.fasta $INDIR/${FILE}_2.scy.sick.fastq > $OUTDIR/${FILE}_2.sai
echo "bwa aln #2 done"

bwa aln -R 1 $REFDIR/NC_000962.3.fasta $INDIR/${FILE}_singles.scy.sick.fastq > $OUTDIR/${FILE}_singles.sai
echo "bwa aln singles done"

bwa sampe -a 1000 -n 1 -N 1 $REFDIR/NC_000962.3.fasta $OUTDIR/${FILE}_1.sai $OUTDIR/${FILE}_2.sai $INDIR/${FILE}_1.scy.sick.fastq $INDIR/${FILE}_1.scy.sick.fastq > $OUTDIR/${FILE}.paired.sam
echo "sampe done"
bwa samse -n 1 $REFDIR/NC_000962.3.fasta $OUTDIR/${FILE}_singles.sai $INDIR/${FILE}_singles.scy.sick.fastq > $OUTDIR/${FILE}.single.sam
echo "samse done"
echo -e "\n2. SAM processing of strain"

echo -e "a) SAM -> BAM\n"
#SAM -> BAM
samtools view -bhSt $REFDIR/NC_000962.3.fasta $OUTDIR/${FILE}.paired.sam -o $OUTDIR/${FILE}.paired.bam
samtools view -bhSt $REFDIR/NC_000962.3.fasta $OUTDIR/${FILE}.single.sam -o $OUTDIR/${FILE}.single.bam
samtools merge $OUTDIR/${FILE}.merged.bam $OUTDIR/${FILE}.paired.bam $OUTDIR/${FILE}.single.bam
samtools sort $OUTDIR/${FILE}.merged.bam -o $OUTDIR/${FILE}.merged.sort.bam
samtools index $OUTDIR/${FILE}.merged.sort.bam

echo -e "b) Mark duplicates\n"
#Mark duplicates
java -Xmx2g -jar /packages/7x/picard/2.9.2/picard.jar MarkDuplicates INPUT=$OUTDIR/${FILE}.merged.sort.bam OUTPUT=$OUTDIR/${FILE}.merged.sort.rmdup.temp.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=10 METRICS_FILE=$OUTDIR/${FILE}.merged.sort.rmdup.metrics.txt VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
java -Xmx2g -jar  /packages/7x/picard/2.9.2/picard.jar AddOrReplaceReadGroups INPUT=$OUTDIR/${FILE}.merged.sort.rmdup.temp.bam OUTPUT=$OUTDIR/${FILE}.merged.sort.rmdup.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT RGLB=8 RGPL=illumina RGPU=1 RGSM=${FILE}
samtools index $OUTDIR/${FILE}.merged.sort.rmdup.bam

echo -e "c) Indel realignment\n"
#Indel realignment
java -Xmx2g -jar /packages/7x/gatk/3.7.0/GenomeAnalysisTK.jar -I $OUTDIR/${FILE}.merged.sort.rmdup.bam -R $REFDIR/NC_000962.3.fasta -T RealignerTargetCreator -maxInterval 20000 -o $OUTDIR/${FILE}.sort.rmdup.indel.intervals
java -Xmx4g -jar /packages/7x/gatk/3.7.0/GenomeAnalysisTK.jar -I $OUTDIR/${FILE}.merged.sort.rmdup.bam -R $REFDIR/NC_000962.3.fasta -T IndelRealigner -targetIntervals $OUTDIR/${FILE}.sort.rmdup.indel.intervals -o $OUTDIR/${FILE}.sort.rmdup.realn.bam

echo -e "d) Base quality recalibration\n"
#Base quality recalibration
bcftools mpileup -Ov -f $REFDIR/NC_000962.3.fasta $OUTDIR/${FILE}.sort.rmdup.realn.bam -o $OUTDIR/${FILE}.sort.rmdup.realn.vcf
#samtools mpileup -uf ./template/fai/${faifile%.fai} samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam > samtools/bcf/${samfile[$c]%.sam}.sort.rmdup.realn.bcf
#bcftools view -vI samtools/bcf/${samfile[$c]%.sam}.sort.rmdup.realn.bcf > samtools/vcf/${samfile[$c]%.sam}.sort.rmdup.realn.vcf
perl /packages/7x/bcftools/1.10.2/bin/vcfutils.pl varFilter -D1000 $OUTDIR/${FILE}.sort.rmdup.realn.vcf > $OUTDIR/${FILE}.sort.rmdup.realn.filt.vcf
java -Xmx4g -jar /packages/7x/gatk/3.7.0/GenomeAnalysisTK.jar -T BaseRecalibrator -R MTB_emp_dataProc/ref/NC_000962.3.fasta --knownSites ERR027294.sort.rmdup.realn.filt.vcf -I ERR027294.sort.rmdup.realn.bam -fP illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ERR027294.sort.rmdup.realn.csv
java -Xmx4g -jar /packages/7x/gatk/3.7.0/GenomeAnalysisTK.jar -R $REFDIR/NC_000962.3.fasta -I $OUTDIR/${FILE}.sort.rmdup.realn.bam -T PrintReads -BQSR $OUTDIR/${FILE}.sort.rmdup.realn.csv -o $OUTDIR/${FILE}.sort.rmdup.recal.bam
samtools index $OUTDIR/${FILE}.sort.rmdup.recal.bam

echo -e "\n3. Statistics for strain ${FILE}\n"
reads_sum=(`wc -l ${OUTDIR}/${FILE}.paired.sam ${OUTDIR}/${FILE}.paired.sam`)
reads_sum[0]=`expr ${reads_sum[4]} - 4`

reads_sum_mapped_paired=(`cut -f5 ${OUTDIR}/${FILE}.paired.sam | sed -e '1,2d;/^0/d' | wc -l`)
reads_sum_mapped_single=(`cut -f5 ${OUTDIR}/${FILE}.single.sam | sed -e '1,2d;/^0/d' | wc -l`)
reads_sum_mapped=`expr $reads_sum_mapped_paired + $reads_sum_mapped_single`

percent_mapping=`echo "scale=2;$reads_sum_mapped / $reads_sum * 100" | bc`%


#export summary files
sumfile[$FILE]=${FILE}.summary
touch $OUTDIR/${sumfile[$FILE]}

echo "RESULT SUMMARY FOR STRAIN ${FILE}: "
echo "RESULT SUMMARY FOR STRAIN ${FILE}:" > $OUTDIR/${sumfile[$FILE]}

echo "Total Reads Number= ${reads_sum[0]}"
echo "Total Reads Number= ${reads_sum[0]}" >>  $OUTDIR/${sumfile[$FILE]}

echo "Initially Mapped Reads Number= $reads_sum_mapped"
echo "Initially Mapped Reads Number= $reads_sum_mapped" >> $OUTDIR/${sumfile[$FILE]}

echo "Initially Mapped Reads Percentage= $percent_mapping"
echo "Initially Mapped Reads Percentage= $percent_mapping" >> $OUTDIR/${sumfile[$FILE]}


rm $OUTDIR/${FILE}_1.sai
rm $OUTDIR/${FILE}_2.sai
rm $OUTDIR/${FILE}_singles.sai
rm $OUTDIR/${FILE}.paired.sam
rm $OUTDIR/${FILE}.single.sam

rm $OUTDIR/${FILE}.paired.bam
rm $OUTDIR/${FILE}.single.bam
rm $OUTDIR/${FILE}.merged.bam

rm $OUTDIR/${FILE}.merged.sort.bam
rm $OUTDIR/${FILE}.merged.sort.bam.bai
rm $OUTDIR/${FILE}.merged.sort.rmdup.temp.bam
rm $OUTDIR/${FILE}.merged.sort.rmdup.bam
rm $OUTDIR/${FILE}.merged.sort.rmdup.bam.bai
#rm $OUTDIR/${FILE}.merged.sort.rmdup.metrics
rm $OUTDIR/${FILE}.sort.rmdup.indel.intervals

echo -e "\nAlignment and SNP calling completed for sample ${FILE}"