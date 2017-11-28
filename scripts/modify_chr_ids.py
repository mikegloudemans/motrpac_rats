#!/usr/bin/python

# Author: Mike Gloudemans

# Modify the IDs of the rat chromosomes to make them numbers, so that they are
# more readable and match with gene annotations

import subprocess

subprocess.check_call('grep chrom /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna | grep CM | cut -d " " -f1,7 | sed s/,/\ /g > /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/chr_ids.txt', shell=True)

id_map = {}
with open("/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/chr_ids.txt") as f:
	for line in f:
		data = line.strip()[1:].split()
		id_map[data[0]] = data[1]

print id_map

with open("/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.numeric_chr_ids.fna", "w") as w:
	with open("/srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats/data/genome/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_genomic.fna") as f:
		for line in f:
			if not line.startswith(">"):
				w.write(line)
				continue
			for im in id_map:
				if im in line:
					line2 = line.replace(im, id_map[im])
					w.write(line2)
					break
