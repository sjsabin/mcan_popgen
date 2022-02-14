#!/bin/bash
#SBATCH -n 4
#SBATCH -p debug
#SBATCH -q wildfire
#SBATCH -t 1
#SBATCH -J freebayes
#SBATCH -o freebayes_01.out
#SBATCH -e freebayes_01.err


LIST=$(cat /home/sjsabin/mcan/data/trim_fq/list.txt)
refdir=/home/sjsabin/mcan/ref 
#Define output path for vcf files.
fbout=/scratch/sjsabin/freebayes
#Set up your modules.
module purge
module load freebayes/1.1.0
module load vcftools/0.1.12b

#freebayes options from Snippy:
	# -p [INTEGER] sets ploidy
		# Here: -p 1 - We are calling SNPs from a haploid organism (but pooled! see below).
	# -P "Report sites if the probability that there is a polymorphism" [?] Default 0.0
		# Here: using default
	# -C [INTEGER] number of supporting observations required to consider a variant
		# Here: we are not using this filter, but sticking with the frequency to do this work for us.
	# -F [FREQUENCY] frequency required to call a variant
		# Here: -F 0.01
	# --min-coverage [INT] " Require at least this coverage to process a site. default: 0"
	# --min-repeat-entropy [0<<1] "haplotype which is called has Shannon entropy less than --min-repeat-entropy which is off by default but can be set to ~1 for optimal genotyping of indels in lower-complexity sequence."
	# --min-mapping-quality [INT] default in Snippy is 60
	# --min-base-quality [INT] default in Snippy is 13
	# --strict-vcf
	
#freebayes options to add:
	# --pooled-continuous (for when you don't know the number of samples in a pool)
	# --pooled-discrete (for when you know the # of samples in a pool)
	

mkdir -p ${fbout}

for FILE in $LIST; do
	echo "###########"
	echo "Starting ${FILE} variant calling..."
	#Define path to bam directory.
	dir=/scratch/sjsabin/bwa_mem/${FILE}
	#Run freebayes on the full bam, but subsample coverage within freebayes.
	echo "### Full coverage variant calling with freebayes downsampling of coverage ###"
	freebayes -f ${refdir}/GCA_000253375.1_ASM25337v1_genomic.fasta -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 13 -m 60 --strict-vcf --pooled-continuous --limit-coverage 25 ${dir}/${FILE}-pe.sorted.dedup.bam > ${fbout}/${FILE}.vcf
	echo "Full coverage variant calling with freebayes downsampling of coverage completed."
	#Run freebayes on the 25x cov subsampled bam
	echo "### 25x coverage subsampled variant calling ###"
	freebayes -f ${refdir}/GCA_000253375.1_ASM25337v1_genomic.fasta -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 13 -m 60 --strict-vcf --pooled-continuous ${dir}/${FILE}-pe.sorted.dedup.subsamp25.bam > ${fbout}/${FILE}.subsamp25.vcf
	echo "25x coverage subsampled variant calling completed."
	#Run vcftools
	module load vcftools/0.1.12b
	# Create VCF file without the sites to exclude listed in the MCAN_regionsToExclude.bed file.
	echo "### Problematic site removal: full (25x)  ###"
	vcftools --vcf ${fbout}/${FILE}.vcf --exclude-bed ${refdir}/MCAN_regionsToExclude_new.bed --out ${fbout}/${FILE}.sitesExcl.vcf --recode
	echo "Problematic site removal: full (25x) completed."
	echo "### Problematic site removal: 25x ###"
	vcftools --vcf ${fbout}/${FILE}.subsamp25.vcf --exclude-bed ${refdir}/MCAN_regionsToExclude_new.bed --out ${fbout}/${FILE}.subsamp25.sitesExcl.vcf --recode
	echo "Problematic site removal: 25x completed."

done

