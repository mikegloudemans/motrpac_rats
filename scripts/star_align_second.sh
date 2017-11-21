#!/bin/bash

module load STAR/2.5.3a

junctions=`echo $1 | sed s/_space_/\ /g`

echo $1
echo $2
echo $3
echo $4
echo $5

STAR --runThreadN 20 \
     --genomeDir $2 \
     --readFilesIn $3 $4 \
     --outFileNamePrefix $5 \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbFileChrStartEnd $junctions \
     --quantMode GeneCounts

