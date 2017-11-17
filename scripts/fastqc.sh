#!/bin/bash
module load fastqc/0.11.5
fastqc --outdir $1 --format fastq $2 $3

