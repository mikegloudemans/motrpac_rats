#!/bin/bash

module load STAR/2.5.3a

# First pass alignment
STAR --runThreadN 20 \
     --genomeDir $1 \
     --readFilesIn $2 $3 \
     --outFileNamePrefix $4 \
     --outSAMtype BAM Unsorted


