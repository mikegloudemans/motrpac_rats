#!/bin/bash

module load STAR/2.5.3a

junctions=`echo $1 | sed s/_space_/ /g`

STAR --runThreadN 20 \
     --genomeDir $2 \
     --readFilesIn $3 $4 \
     --outFileNamePrefix $5 \
     --outSAMtype BAM Unsorted \
     --sjdbFileChrStartEnd $junctions \
     --quantMode GeneCounts

