#!/bin/bash

module load STAR/2.5.3a

STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir $1 \
     --genomeFastaFiles $2 \
     --sjdbGTFfile $3 \
     --sjdbOverhang 100  # readLength - 1
