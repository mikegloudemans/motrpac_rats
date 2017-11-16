# De-multiplexing https://www.biostars.org/p/205712/

module load bcl2fastq2/2.17.1.14

# Example adapted from Laure Fresard:
set -o nounset -o pipefail
# SCG3

runs="171103_D00125R_0249_AHYFKVBCXY 171106_D00125R_0250_BHYFLGBCXY" 

for runName in $runs; do
	runDirectory=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/bcl/rsync

	user=mgloud@stanford.edu
	seqDir=$runDirectory/$runName
	outDir=/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/fastq/$runName
	JOB=${runName}_bcl2fastq
	logOUT=$JOB.out
	logError=$JOB.error
	rm -f $logOUT $logError

	# 32G to 12G
	bcl2fastq \
	    --runfolder-dir $seqDir \
	    --output-dir $outDir \
	    --minimum-trimmed-read-length 50 \
	    --stats-dir ./stats \
	    --reports-dir ./reports
done
