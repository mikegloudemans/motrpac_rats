#!/bin/bash

java -jar /srv/gsfs0/software/picard-tools/2.8.0/picard.jar CollectRnaSeqMetrics \
     REF_FLAT=$1 \
     INPUT=$2 \
     OUTPUT=$3 \
     STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND


