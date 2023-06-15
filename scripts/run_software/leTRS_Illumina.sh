#!/bin/bash

# import path to the different software
source percorsi.sh
#conda activate LeTRS

sample=$1

cd $sample

fastq1=$(realpath trimmed*R1*.fastq.gz);
fastq2=$(realpath trimmed*R2*.fastq.gz);

echo "$(pwd)"
echo "fastq1: $fastq1 ; fastq2: $fastq2"

# run LeTRS using illumina mode 
# CHANGE TO PRIMER PATH
perl $LeTRS -t 4 -extractfasta -pool 0 -mode illumina -fq $fastq1:$fastq2 -o LeTRS_output_cov0 -primer_bed /opt/references/nCoV-2019/V3/nCoV-2019.primer.bed -covcutoff 0 

cd LeTRS_output/results

perl $LeTRS_plot -count 1 -i known_junction.tab
mv leader-TRS.pdf leader-TRS_known.pdf

perl $LeTRS_plot -count 1 -i novel_junction.tab
mv leader-TRS.pdf leader-TRS_novel.pdf
