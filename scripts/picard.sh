#!/bin/bash

module load picard-tools/2.13.2


java -jar /srv/gsfs0/software/picard-tools/2.13.2/picard.jar CollectRnaSeqMetrics \
     REF_FLAT=$1 \
     INPUT=$2 \
     OUTPUT=$3 \
     STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND


