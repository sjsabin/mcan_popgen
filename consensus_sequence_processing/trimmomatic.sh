#!/bin/bash

# Trimmomatic to trim adaptors and low quality bases from NGS read files

# loops over all fastq files sored in variable $read_directory
read_directory=/home/sjsabin/mcan/data/raw_fq/

for forward in $(find "${read_directory}" -name "*_1.fastq")
do
# extract read ID and store in variable $stem
stem=$(basename "${forward}" _1.fastq)
# extract directory of forward reads and store in variable $directory
directory=$(dirname "${forward}")/
# build path to reverse reads and store in variable $reverse
reverse="${directory}${stem}_2.fastq"
mkdir -p ~/mcan/trimmomatic_out_TruSeq3-PE-2-min75/"${stem}"
# output according to path set in variable $out_path
out_path=/home/sjsabin/mcan/data/trim_fq/"${stem}"/

# trimmomatic settings
# adding the appropriate suffix to the basename var for each of 4 output files.
# 'p' for paired, 'u' for unpaird, 'f' for foward, 'r' for reverse
  p_output_f="${out_path}${stem}_1.trim.fastq.gz"
  u_output_f="${out_path}${stem}_1_unpaired.trim.fastq.gz"
  p_output_r="${out_path}${stem}_2.trim.fastq.gz"
  u_output_r="${out_path}${stem}_2_unpaired.trim.fastq.gz"
# trimmomatic arguments
# number of bases to be analyzed at a time (the window)
  window=4
# the average quality required across window to be retained
  min_quality=20
# the minimum length of a read to be kept
  min_length=25
# adapters used to create the PE sequences (file should be included with trimmomatic install)
  adapters=/packages/7x/trimmomatic/0.33/adapters/
# how many mismatches are tolerated to allow a seed
  mismatch=2
# minimum match for palindrome alignment
  palindrome=40
# minimum match to remove adapter sequence
  adapter_trim=15

# running Trimmomatic based on the variables created above
 trimmomatic PE $forward $reverse $p_output_f $u_output_f $p_output_r $u_output_r\
 ILLUMINACLIP:${adapters}:${mismatch}:${palindrome}:${adapter_trim}\
 SLIDINGWINDOW:${window}:${min_quality}\
 MINLEN:${min_length}
done
