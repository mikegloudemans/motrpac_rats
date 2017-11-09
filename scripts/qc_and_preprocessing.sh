#!/bin/bash
# Author: Mike Gloudemans
# Date created: 11/9/2017

# Pipeline adapted from https://github.com/zaczap/bios201/tree/master/Workshop2
# Look more into RNA-seq best practices. Find pipeline for quantification

module load STAR/2.5.3a

##########################
# Part 0: Fastq conversion
##########################

# TODO: If not done already

##########################
# Part 1: QC
##########################

# For each sample

	# Run fastQC
	fastqDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/fastq
	fastqc --outdir /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/fastqc --format fastq $fastqDir/R1 $fastqDir/R2

	# Run cutadapt to trim adapters and filter short/low-quality reads
	# TODO

	# We're choosing not to remove PCR duplicates, for reasons described in
	# https://www.nature.com/articles/srep25533


##########################
# Part 2: Alignment
##########################

gunzip /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna.gz
gunzip /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.chr.gtf.gz

# Generate genome index
homeDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats
genomeDir=$homeDir/data/genome/STAR_index
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir $genomeDir \
     --genomeFastaFiles /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr12.fna \
     --sjdbGTFfile /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.chr.gtf.gz \
     --sjdbOverhang 100  # readLength - 1 TODO TODO TODO fix based on the actual read length we find

#     --genomeFastaFiles /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna \

# For each sample...

	# Map reads with STAR

	# First pass alignment
	STAR --runThreadN 20 \
	     --genomeDir $genomeDir \
	     --readFilesIn $homeDir/data/fastq/path/to/read1 $homeDir/data/fastq/path/to/read2 \
	     --outFileNamePrefix $homeDir/data/bam/bam_pass1/samplename_ \
	     --outSAMtype BAM Unsorted \

# For each sample...

	# Second pass alignment, using additional inferred splice junctions
	# from all samples
	junctions=`ls $homeDir/data/bam/bam_pass1/*_SJ.out.tab`
	STAR --runThreadN 20 \
	     --genomeDir $genomeDir \
	     --readFilesIn $homeDir/data/fastq/path/to/read1 $homeDir/data/fastq/path/to/read2 \
	     --outFileNamePrefix bam_pass2/Norm1_ \
	     --outSAMtype BAM Unsorted \
	     --sjdbFileChrStartEnd $junctions \
	     --quantMode GeneCounts

	# Run Picard to generate mapping metrics
	java -jar $PICARD CollectRnaSeqMetrics \
	     REF_FLAT=annotation/refFlat.chr10.txt \
	     INPUT=bam_pass2/Norm1_Aligned.out.sorted.bam \
	     OUTPUT=bam_pass2/Norm1_Aligned.out.sorted.metrics.txt \
	     STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND

# At this point, we're ready for exploration with R/DEseq2.
