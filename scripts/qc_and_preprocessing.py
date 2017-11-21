#!/usr/bin/python
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

##########
# Settings
##########

create_index = True	# Create fresh index of the genome?

##############################

import subprocess
import time
import glob

dataDir = "/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data"
logDir = "/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/log"

def wait_for_jobs(job_list):
# qstat repeatedly until all de-multiplexing jobs are done
	jobs_left = [j.strip() for j in job_list[:]]
	while len(jobs_left) > 0:
		status = subprocess.check_output("qstat | awk '{{print $1}}' | tail -n +3".format(**locals()), shell=True).split()
		remaining = []
		for job in jobs_left:
			if job in status:
				remaining.append(job)

		jobs_left=remaining
		time.sleep(5)

##########################
# Part -1: Indexing genome
##########################

# Unzip these files if necessary
subprocess.call("gunzip {dataDir}/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna.gz > /dev/null".format(**locals()), shell=True)
subprocess.call("gunzip {dataDir}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf.gz > /dev/null".format(**locals()), shell=True)

# Generate genome index
index_job_id = []

if create_index:
	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu -o {logDir} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_index.sh {dataDir}/genome/STAR_index {dataDir}/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna {dataDir}/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.90.gtf".format(**locals()), shell=True)
	index_job_id.append(id)

##########################
# Part 0: Fastq conversion
##########################

# De-multiplexing
runs=["171103_D00125R_0249_AHYFKVBCXY", "171106_D00125R_0250_BHYFLGBCXY"]
jobs_left = []
for runName in runs:
	id=subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=32G -l h_rt=24:00:00 -m ea -M mgloud@stanford.edu -o {logDir} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastq_conversion.sh {runName}".format(**locals()), shell=True)
	jobs_left.append(id)

wait_for_jobs(jobs_left)

##########################
# Part 1: QC
##########################

qcDir = "/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/preprocessing/fastqc"
subprocess.call("mkdir -p {qcDir}".format(**locals()), shell=True)

samples = glob.glob("{dataDir}/fastq/1*/*".format(**locals()))
samples = [s for s in samples if "trimmed" not in s and "R1_001" in s]

# Run fastQC on every lane/sample combination
fastqc_jobs = []
for f1 in samples:
	f2 = f1.replace("_R1_001", "_R2_001")

	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 -o {logDir} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/fastqc.sh {qcDir} {f1} {f2}".format(**locals()), shell=True)
	fastqc_jobs.append(id)

# Run cutadapt on individual samples to trim adapters and filter short/low-quality reads
# Adapters chosen according to https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html
cutadapt_jobs = []
for f1 in samples:
        f2 = f1.replace("_R1_001", "_R2_001")
	trimmedF1 = f1.replace(".fastq", ".trimmed.fastq")	
	trimmedF2 = f2.replace(".fastq", ".trimmed.fastq")

	base=f1.strip().split("/")[-1]

	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=32G -l h_rt=1:00:00 -o /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/preprocessing/cutadapt/{base}.cutadapt.log cutadapt.sh {trimmedF1} {trimmedF2} {f1} {f2}".format(**locals()), shell=True)
	cutadapt_jobs.append(id)

# Combine samples across runs/lanes
wait_for_jobs(cutadapt_jobs)
samples = glob.glob("{dataDir}/fastq/1*/*".format(**locals()))
samples = [s for s in samples if "trimmed" in s and "R1_001" in s]
subprocess.call("mkdir -p {dataDir}/fastq/combined_qced".format(**locals()), shell=True)
subprocess.call("rm -f {dataDir}/fastq/combined_qced/*fastq*", shell=True)
for sample in samples:
	sample_mod = sample.strip().split("/")[1].replace("L001_", "").replace("L002_", "")
	subprocess.check_call("zcat {sample} >> {dataDir}/fastq/combined_qced/{sample_mod}".format(**locals()), shell=True)

# We're choosing not to remove PCR duplicates, for reasons described in
# https://www.nature.com/articles/srep25533

# Count reads per sample
subprocess.check_call("bash /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/read_count.sh > /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/preprocessing/read_counts.txt", shell=True)

##########################
# Part 2: Alignment
##########################

subprocess.check_call("mkdir -p {dataDir}/bam/bam_pass1/".format(**locals()), shell=True)
subprocess.check_call("mkdir -p {dataDir}/bam/bam_pass2/".format(**locals()), shell=True)

wait_for_jobs(index_job_id)

# First pass alignment, to infer unknown splice junctions
# For each sample...
align_first_jobs = []
samples = glob.glob("{dataDir}/fastq/combined/*trimmed*".format(**locals()))

for f1 in samples:
        f2 = f1.replace("_R1_001", "_R2_001")
	abbrev_f1 = f1.replace("R1_001.trimmed.fastq", "")

	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu -o /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/preprocessing/star/pass1 /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_align_first.sh {dataDir}/genome/STAR_index {f1} {f2} {dataDir}/bam/bam_pass1/{abbrev_f1}")
	align_first_jobs.append(id)

# Second pass alignment, using additional inferred splice junctions
# from all samples
wait_for_jobs(align_first_jobs)
junctions = glob.glob("{dataDir}/bam/bam_pass1/*_SJ.out.tab").format(**locals()).strip().split()
junctions = "_space_".join(junctions)
# For each sample...
align_second_jobs = []
for f1 in samples:
        f2 = f1.replace("_R1_001", "_R2_001")
        abbrev_f1 = f1.replace("R1_001.trimmed.fastq", "")

	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=80G -l h_rt=12:00:00 -m ea -M mgloud@stanford.edu -o {logDir} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/star_align_second.sh {junctions} {dataDir}/genome/STAR_index {f1} {f2} {dataDir}/bam/bam_pass2/{abbrev_f1}".format(**locals()), shell=True)

# Sort BAMs

# Run Picard to generate mapping metrics
wait_for_jobs(align_second_jobs)
picard_jobs = []
#module load picard-tools/2.8.0

samples = glob.glob("{dataDir}/bam/bam_pass2/*".format(**locals()))
for f in samples:
	base = f.strip().split("/")[-1]
	id = subprocess.check_output("qsub -j y -terse -cwd -S /bin/bash -l h_vmem=20G -l h_rt=4:00:00 -m ea -M mgloud@stanford.edu -o {logDir} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/scripts/picard.sh {dataDir}/refFlat.txt {f} /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/output/preprocessing/picard/{base}.picard.log".format(**locals()), shell=True)
	picard_jobs.append(id)

# Run MultiQC
wait_for_jobs(picard_jobs)
wait_for_jobs(fastqc_jobs)

# Unzip FastQC output for processing with MultiQC
#for z in `ls $1/*.zip`; do unzip -d $1 $z; done

# Then we're ready for exploration with R/DEseq2.
