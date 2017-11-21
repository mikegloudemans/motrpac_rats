#!/usr/bash

for f in `ls -d ../data/fastq/combined/*`; do 
	s=`wc -l $f | awk '{print $1}'`
	c=`expr $s / 4`
	echo $c $f
done
