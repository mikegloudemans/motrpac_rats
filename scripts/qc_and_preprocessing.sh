#!/bin/bash
# Author: Mike Gloudemans
# Date created: 11/9/2017

module load STAR/2.5.3a

# Generate genome index
genomeDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/STAR_index
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir $genomeDir \
     --genomeFastaFiles  \
     --sjdbGTFfile annotation/gencode.v19.annotation_chr10.gtf \
     --sjdbOverhang 100  # readLength - 1



