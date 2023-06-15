#!/bin/bash

# import path to the different software
source percorsi.sh

folder=$1

fastq1=$folder/trimmed.R1.fastq.gz
fastq2=$folder/trimmed.R2.fastq.gz

mkdir -p $folder/sgDI_results/R1
mkdir -p $folder/sgDI_results/R2

reference="/opt/references/SARS-CoV2/SARS-CoV-2.fasta"
annotation_file="/opt/references/SARS-CoV2/NC_045512.2_sgRNA.fasta"

# run sgDI on the two fastq separatly
python $sgDI -r $annotation_file -o $folder/sgDI_results/R1 $reference $fastq1
python $sgDI -r $annotation_file -o $folder/sgDI_results/R2 $reference $fastq2 
