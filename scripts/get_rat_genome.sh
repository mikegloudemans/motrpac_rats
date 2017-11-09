#!/bin/bash

# Download rat genome
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Rattus_norvegicus/latest_assembly_versions/GCA_000001895.4_Rnor_6.0 .

# Check MD5 sums. The only difference should be that our directory has
# a few files that weren't part of the original checksummed directory.
find -type f -exec md5sum '{}' \; > md5sum.txt
diff <(sort md5sum.txt) <(sort md5checksums.txt)
