#!/bin/bash
#
# Author: Kun Sun @ SZBL (hellosunking@foxmail.com)
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 3 ]
then
	echo -e "\nUsage: $0 <genome folder|fa> <refSeq.txt|gene_anno.gff|null> <index.identifier> [thread=auto]\n" >/dev/stderr
	echo -e "This is a utility program of Msuite, designed to build genome references." >/dev/stderr
	echo -e "Please refer to README file for more information.\n" >/dev/stderr
	exit 2
fi

PRG=`dirname $0`
indexDIR=$PRG/../index
#indexDIR=$PRG/index
lambdaGenome=$PRG/lambda.genome.fa
pUC19Genome=$PRG/pUC19.fa

echo -e "\e[34mInput genome:      $1"
echo -e "RefSeq annotation: $2"
echo -e "Lambda genome:     $lambdaGenome"
echo -e "pUC19 genome:      $pUC19Genome"
echo -e "Index identifier:  $3"
echo -e "Index location:    $indexDIR/$3\e[39m\n"

id=$3
if [ $# -gt 3 ]
then
	thread=$4
else
	thread=`grep processor /proc/cpuinfo | wc -l`
fi
echo "=> INFO: Threads used: $thread"

## check bowtie2
bb=`which bowtie2-build 2>/dev/null`
if [ $? != 0 ] || [ -z "$bb" ]; then
	echo -e "\e[31mFatal error: Could not find 'bowtie2-build' in your path!\e[39m" >/dev/stderr
	exit 1
fi
bbver=`$bb --version | perl -ne 'print $1 if /bowtie2-build\S* version ([\d\.]+)/'`
if [ -z "$bbver" ]; then
	ver=unknown
fi
echo -e "\n\e[32mINFO: bowtie2-build (version $bbver) found at '$bb'.\e[39m\n"

## check hisat2
hb=`which hisat2-build 2>/dev/null`
if [ $? != 0 ] || [ -z "$hb" ]; then
	echo -e "\e[31mFatal error: Could not find 'hisat2-build' in your path!\e[39m" >/dev/stderr
	exit 1
fi
hbver=`$hb --version | perl -ne 'print $1 if /hisat2-build\S* version ([\d\.]+)/'`
if [ -z "$hbver" ]; then
	ver=unknown
fi
echo -e "\n\e[32mINFO: hisat2-build (version $hbver) found at '$hb'.\e[39m\n"

## prepare directories
mkdir -p $indexDIR/$id
mkdir -p $indexDIR/$id/indices/bowtie2
mkdir -p $indexDIR/$id/indices/hisat2
mkdir -p $indexDIR/$id/fasta

echo "Processing gene annotation ..."
if [ $2 == "null" ]
then
	echo "WARNING: null is specified as gene annotation file, I will skip this step!"
else
	if [[ $2 =~ g[tf]f ]]
	then
		echo "Info: processing file in GFF format ..."
		perl $PRG/process.gff.pl $2 >$indexDIR/$id/tss.ext.bed &
	else
		echo "Info: processing file in UCSC's refGene format ..."
		perl $PRG/process.refGene.pl $2 >$indexDIR/$id/tss.ext.bed &
	fi
fi

echo "Preprocessing genome ..."
perl $PRG/process.genome.pl $1 $lambdaGenome,$pUC19Genome $indexDIR/$id/
##ln -s watson.fa $indexDIR/$id/genome.fa

echo "Building 4-letter indices for bowtie2 ..."
#$bb --threads $thread $indexDIR/$id/CG2TG.fa $indexDIR/$id/indices/bowtie2/m4 >/dev/null
echo "Building 3-letter indices for bowtie2 ..."
#$bb --threads $thread $indexDIR/$id/C2T.fa   $indexDIR/$id/indices/bowtie2/m3 >/dev/null
touch $indexDIR/$id/indices/bowtie2.ready

echo "Building 4-letter indices for hisat2 ..."
#$hb -p $thread $indexDIR/$id/CG2TG.fa $indexDIR/$id/indices/hisat2/m4 >/dev/null
echo "Building 3-letter indices for hisat2 ..."
$hb -p $thread $indexDIR/$id/C2T.fa   $indexDIR/$id/indices/hisat2/m3 >/dev/null
touch $indexDIR/$id/indices/hisat2.ready

wait
rm -f $indexDIR/$id/CG2TG.fa $indexDIR/$id/C2T.fa
echo
echo -e "\e[35mDone. Now you can use \"-x $3\" to use this genome in Msuite2.\e[39m"
echo

