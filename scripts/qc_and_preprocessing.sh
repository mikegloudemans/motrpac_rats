#!/bin/bash
# Author: Mike Gloudemans
# Date created: 11/9/2017

# Pipeline adapted from https://github.com/zaczap/bios201/tree/master/Workshop2
# Look more into RNA-seq best practices. Find pipeline for quantification

dataDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data

##########################
# Part -1: Indexing genome
##########################

gunzip $dataDir/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna.gz
gunzip $dataDir/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.chr.gtf.gz

# Generate genome index
qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_index.sh $dataDir/genome/STAR_index $dataDir/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna $dataDir/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.chr.gtf.gz


##########################
# Part 0: Fastq conversion
##########################

# De-multiplexing
runs="171103_D00125R_0249_AHYFKVBCXY 171106_D00125R_0250_BHYFLGBCXY"
for runName in $runs; do
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=24:00:00 -m ea -M mgloud@stanford.edu /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastq_conversion.sh $runName
done;

# TODO: Add waiting step here, for qsub to finish. Not
# sure what's the clean way to do this

# Merging like samples
#for f in `ls -d $dataDir/fastq/171103*/*`; do
for f in `ls -d $dataDir/fastq/17*/*`; do
	sample=`basename $f | sed s/.gz//g | sed s/L001_//g | sed s/L002_//g`
	zcat $f >> $dataDir/fastq/combined/$sample
done


##########################
# Part 1: QC
##########################

qcDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/fastqc

# For each sample
for f1 in `ls $dataDir/fastq/combined | grep R1_001`; do

	f2=`echo $f1 | sed s/_R1_001/_R2_001/g`

	# Run fastQC
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastqc.sh $qcDir/pre_cutadapt $dataDir/fastq/combined/$f1 $dataDir/fastq/combined/$f2

	# Run cutadapt to trim adapters and filter short/low-quality reads
	# TODO: Verify that these adapters worked correctly
	# Currently chosen according to https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
	trimmedF1=`echo $f1 | sed s/fastq/trimmed.fastq/g`
	trimmedF2=`echo $f2 | sed s/fastq/trimmed.fastq/g`
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 cutadapt.sh $dataDir/fastq/combined/$trimmedF1 $dataDir/fastq/combined/$trimmedF2 $dataDir/fastq/combined/$f1 $dataDir/fastq/combined/$f2

	# We're choosing not to remove PCR duplicates, for reasons described in
	# https://www.nature.com/articles/srep25533

	# Run fastQC again and see if anything changed this time
	qsub -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastqc.sh $dataDir/fastq/combined/$f1 $dataDir/fastq/combined/$f2

done


##########################
# Part 2: Alignment
##########################

# First pass alignment, to infer unknown splice junctions
# For each sample...
for f1 in `ls $dataDir/fastq/combined/*trimmed* | grep R1_001`

	f2=`echo $f1 | sed s/R1_001/R2_001/g`
	
        trimmedF1=`echo $f1 | sed s/R1_001.trimmed.fastq//g`
	# First pass alignment
	STAR --runThreadN 20 \
	     --genomeDir $dataDir/genome/STAR_index \
	     --readFilesIn $f1 $f2 \
	     --outFileNamePrefix $dataDir/bam/bam_pass1/$trimmedF1 \
	     --outSAMtype BAM Unsorted

	qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu star_align_first.sh $dataDir/genome/STAR_index $f1 $f2 $dataDir/bam/bam_pass1/$trimmedF1

# Second pass alignment, using additional inferred splice junctions
# from all samples
junctions=`ls $dataDir/bam/bam_pass1/*_SJ.out.tab | xargs | sed s/\ /_space_/g`
# For each sample...
for f1 in `ls $dataDir/fastq/combined/*trimmed* | grep R1_001`

        f2=`echo $f1 | sed s/R1_001/R2_001/g`

        trimmedF1=`echo $f1 | sed s/R1_001.trimmed.fastq//g`
	STAR --runThreadN 20 \
	     --genomeDir $dataDir/genome/STAR_index \
	     --readFilesIn $f1 $f2 \
	     --outFileNamePrefix $dataDir/bam/bam_pass2/$trimmedF1 \
	     --outSAMtype BAM Unsorted \
	     --sjdbFileChrStartEnd $junctions \
	     --quantMode GeneCounts

	qsub -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu star_align_second.sh $junctions $dataDir/genome/STAR_index $f1 $f2 $dataDir/bam/bam_pass2/$trimmedF1

	# Run Picard to generate mapping metrics
	java -jar $PICARD CollectRnaSeqMetrics \
	     REF_FLAT=annotation/refFlat.chr10.txt \
	     INPUT=bam_pass2/Norm1_Aligned.out.sorted.bam \
	     OUTPUT=bam_pass2/Norm1_Aligned.out.sorted.metrics.txt \
	     STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND

# At this point, we're ready for exploration with R/DEseq2.
