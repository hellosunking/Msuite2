#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : May 6, 2020

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Msuite.final.bam> [output.bed=stdout]\n";
	print STDERR "\nThis program is designed to translate the bam file (paired-end mode) to bed file.\n\n";
	exit 1;
}

my $o = $ARGV[1] || '/dev/stdout';
open OUT, ">$o" or die( "$!" );

open IN, "samtools view $ARGV[0] |" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/;	##SRR8885067.sra.142787   163     chr17   35511165        3       50M     =       35512033        918	seqseq
	if( $l[6] eq '=' || $l[8]>=0 ) {
		-- $l[3];	## change to 0-base
		print join("\t", $l[2], $l[3], $l[3]+$l[8], $l[0], $l[4]), "\n";
	}
}
close IN;

close OUT;

