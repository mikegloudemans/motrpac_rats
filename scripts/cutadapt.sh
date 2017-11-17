#!/bin/bash

module load python/2.7
module load cutadapt/1.8.1

cutadapt -q 10 --minimum-length 30 \
	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o $1 \
	-p $2 \
	$3 $4

