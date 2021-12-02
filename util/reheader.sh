#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 2 ]
then
	echo "Usage: $0 <in.bam> <out.bam> [header=auto]" > /dev/stderr
	echo "This program is designed to reheader the BAM file (for genomes where chromosome ids do not start with 'chr')."
	echo "You can provide the header file in SAM format, or I will automatically build it using the raw BAM file."
	exit 2
fi

PRG=`dirname $0`

if [ $# -gt 2 ]	## the user provide the sam header
then
	samtools reheader $3 $1 > $2
else	## generate header by removing chr prefix
	samtools view -H $1 > raw.header.sam
	cat raw.header.sam | perl -ne 's/SN:chr/SN:/; print' > alt.header.sam
	rawSize=`stat -c %s raw.header.sam`
	altSize=`stat -c %s alt.header.sam`
	if [ $rawSize -ne $altSize ]
	then
		samtools reheader modified.header.sam $1 > $2
	fi
fi

