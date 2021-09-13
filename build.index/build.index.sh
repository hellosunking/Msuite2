#!/bin/bash
#
# Author: Kun Sun @ SZBL (hellosunking@foxmail.com)
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 3 ]
then
	echo -e "\nUsage: $0 <genome folder|fa> <refSeq.txt|gene_anno.gff> <index.identifier> [thread=auto]\n" >/dev/stderr
	echo -e "This is a utility program of Msuite, designed to build genome references." >/dev/stderr
	echo -e "Please refer to README file for more information.\n" >/dev/stderr
	exit 2
fi

PRG=`dirname $0`
indexDIR=$PRG/../index
#indexDIR=$PRG/index
lambdaGenome=$PRG/lambda.genome.fa

id=$3
if [ $# -gt 3 ]
then
	thread=$4
else
	thread=`grep processor /proc/cpuinfo | wc -l`
fi
echo "=> INFO: Threads used: $thread"

bb=`which bowtie2-build 2>/dev/null`
if [ -z "$bb" ]; then
	echo -e "\e[31mFatal error: Could not find 'bowtie2-build' in your path!\e[39m" >/dev/stderr
	exit 1
fi

ver=`$bb --version | perl -ne 'print $1 if /bowtie2-build\S* version ([\d\.]+)/'`
if [ -z "$ver" ]; then
	ver=unknown
fi

echo -e "\n\e[32mINFO: bowtie2-build (version $ver) found at '$bb'.\e[39m\n"
echo -e "\e[34mInput genome:      $1"
echo -e "RefSeq annotation: $2"
echo -e "Lambda genome:     $lambdaGenome"
echo -e "Index identifier:  $3"
echo -e "Index location:    $indexDIR/$3\e[39m\n"

echo "Preprocessing genome ..."
## prepare directories
mkdir -p $indexDIR/$id
mkdir -p $indexDIR/$id/indices
mkdir -p $indexDIR/$id/fasta

perl $PRG/process.genome.pl $1 $lambdaGenome $indexDIR/$id/
##ln -s watson.fa $indexDIR/$id/genome.fa

echo "Building 4-letter indices ..."
$bb --threads $thread --large-index \
	$indexDIR/$id/CG2TG.fa $indexDIR/$id/indices/m4 >/dev/null

echo "Building 3-letter indices ..."
$bb --threads $thread --large-index \
	$indexDIR/$id/C2T.fa   $indexDIR/$id/indices/m3 >/dev/null

echo "Processing gene annotation ..."
if [[ $2 =~ gff ]]
then
	echo "Info: processing file in GFF format ..."
	perl $PRG/process.gff.pl $2 >$indexDIR/$id/tss.ext.bed
else
	echo "Info: processing file in UCSC's refGene format ..."
	perl $PRG/process.refGene.pl $2 >$indexDIR/$id/tss.ext.bed
fi

rm -f $indexDIR/$id/CG2TG.fa $indexDIR/$id/C2T.fa
echo
echo -e "\e[35mDone. Now you can use \"-x $3\" to use this genome in Msuite2.\e[39m"
echo

