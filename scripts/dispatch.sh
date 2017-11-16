#!/bin/bash

# Fastq conversion
qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=24:00:00 -m ea -M mgloud@stanford.edu fastq_conversion.sh
