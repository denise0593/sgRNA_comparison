#!/bin/bash

# author: Denise Lavezzari
# 07/06/2022

#conda activate periscope
# import path to the different software
source percorsi.sh

cd $1

fastq1=$(realpath *1.fastq.gz);
fastq2=$(realpath *2.fastq.gz);

name=$(basename $(pwd))

# run periscope
periscope --fastq $fastq1 $fastq2 --output-prefix $name --artic-primers V3 --resources $periscope_resource --technology illumina --threads 3

