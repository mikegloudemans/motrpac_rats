#!/bin/bash
# Author: Mike Gloudemans
# Date created: 11/9/2017

# Pipeline adapted from https://github.com/zaczap/bios201/tree/master/Workshop2
# Look more into RNA-seq best practices. Find pipeline for quantification

# NOTE: I've currently outlined a sketch of how to make this work in sequence,
# having qsub jobs diverge and then waiting for all of them to complete after
# each step. However, if scaling, this would be faster if we just went straight
# from fastq creation to trimming to alignment, just launching all the QC 
# and other stuff independently instead of waiting for every job to complete at
# the end of each step

dataDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data
logDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/log

##########################
# Part -1: Indexing genome
##########################

gunzip $dataDir/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna.gz
gunzip $dataDir/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz

# Generate genome index
qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_index.sh $dataDir/genome/STAR_index $dataDir/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna $dataDir/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf

# (index_job_id=get job id from qsub)

##########################
# Part 0: Fastq conversion
##########################

# De-multiplexing
runs="171103_D00125R_0249_AHYFKVBCXY 171106_D00125R_0250_BHYFLGBCXY"
for runName in $runs; do
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=24:00:00 -m ea -M mgloud@stanford.edu -o $logDir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastq_conversion.sh $runName
	# (track list of jobs submitted here for multiplexing)
done;

# qstat repeatedly until all multiplexing jobs are done

# Merging same samples across runs/lanes
#for f in `ls -d $dataDir/fastq/171103*/*`; do
for f in `ls -d $dataDir/fastq/17*/* | grep -v trimmed`; do
	sample=`basename $f | sed s/.gz//g | sed s/L001_//g | sed s/L002_//g`
	zcat $f >> $dataDir/fastq/combined/pre_qc/$sample
done

# Count reads per sample
bash /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/read_count.sh > /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/read_counts.txt
# (need to wait until previous step is done first)


##########################
# Part 1: QC
##########################

qcDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/fastqc

# For each sample, for each lane:
for f1 in `ls $dataDir/fastq/1*/* | grep -v trimmed | grep R1_001`; do

	f2=`echo $f1 | sed s/_R1_001/_R2_001/g`

	# Run fastQC on every lane/sample combination
	mkdir -p $qcDir/pre_cutadapt
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 -o $logDir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastqc.sh $qcDir/pre_cutadapt $f1 $f2

	# Run cutadapt to trim adapters and filter short/low-quality reads
	# TODO: Verify that these adapters worked correctly
	# Currently chosen according to https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

	trimmedF1=`echo $f1 | sed s/\\\.fastq/\\.trimmed.fastq/g`
	trimmedF2=`echo $f2 | sed s/\\\.fastq/\\.trimmed.fastq/g`
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 -o $logDir cutadapt.sh $trimmedF1 $trimmedF2 $f1 $f2

	# (get job ids for cutadapt and store in list)

	# See what fastQC did for us
	mkdir -p $qcDir/post_cutadapt
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 -o $logDir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastqc.sh $qcDir/post_cutadapt $trimmedF1 $trimmedF2

	# We're choosing not to remove PCR duplicates, for reasons described in
	# https://www.nature.com/articles/srep25533

	# Unzip FastQC output for processing with MultiQC
	for z in `ls $qcDir/*.zip`; do unzip $z; done

done

# Now, after QCing and running cutadapt, recombine all samples

for f in `ls -d $dataDir/fastq/17*/* | grep trimmed`; do
	sample=`basename $f | sed s/\\\.gz//g | sed s/L001_//g | sed s/L002_//g`
	zcat $f >> $dataDir/fastq/combined/$sample
done



##########################
# Part 2: Alignment
##########################

# (qstat and loop until all cutadapt jobs have finished)

# First pass alignment, to infer unknown splice junctions
# For each sample...
for f1 in `ls $dataDir/fastq/combined/*trimmed* | grep R1_001`; do

	f2=`echo $f1 | sed s/R1_001/R2_001/g`
	
        trimmedF1=`basename $f1 | sed s/R1_001.trimmed.fastq//g`
	qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu -o $logDir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_align_first.sh $dataDir/genome/STAR_index $f1 $f2 $dataDir/bam/bam_pass1/$trimmedF1
done

# Second pass alignment, using additional inferred splice junctions
# from all samples
junctions=`ls $dataDir/bam/bam_pass1/*_SJ.out.tab | xargs | sed s/\ /_space_/g`
# For each sample...
for f1 in `ls $dataDir/fastq/combined/*trimmed* | grep R1_001`; do

        f2=`echo $f1 | sed s/R1_001/R2_001/g`

        trimmedF1=`basename $f1 | sed s/R1_001.trimmed.fastq//g`
	qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu -o $logDir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_align_second.sh $junctions $dataDir/genome/STAR_index $f1 $f2 $dataDir/bam/bam_pass2/$trimmedF1

	# TODO:
	# Run Picard to generate mapping metrics
	# java -jar $PICARD CollectRnaSeqMetrics \
	#     REF_FLAT=annotation/refFlat.chr10.txt \
	#     INPUT=bam_pass2/Norm1_Aligned.out.sorted.bam \
	#     OUTPUT=bam_pass2/Norm1_Aligned.out.sorted.metrics.txt \
	#     STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
done

# Then we're ready for exploration with R/DEseq2.
