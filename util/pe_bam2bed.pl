#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : May 6, 2020

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Msuite.final.bam> [output.bed=stdout] [min.mapQ=0] [r1.head.cut=0] [r2.head.cut=0]\n";
	print STDERR "\nThis program is designed to translate the bam file (paired-end mode) to bed file.\n\n";
	exit 1;
}

my $mapQ = $ARGV[2] || 0;
my $r1   = $ARGV[3] || 0;
my $r2   = $ARGV[4] || 0;

my $o = $ARGV[1] || '/dev/stdout';
if( $o =~ /\.gz$/ ) {
	open OUT, "| gzip >$o" or die( "$!" );
} else {
	open OUT, ">$o" or die( "$!" );
}

open IN, "samtools view -@ 2 $ARGV[0] |" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/;	##SRR8885067.sra.142787   163     chr17   35511165        3       50M     =       35512033        918	seqseq
	if( $l[4]>=$mapQ && $l[6] eq '=' && $l[8]>=0 ) {
		-- $l[3];	## change the left-most end to 0-base
		my $strand = '.';
		my ( $s, $e );
		if( $l[-1] =~ /XG:Z:CT/ ) {	## watson strand
			$strand = '+';
			$s = $l[3] - $r1;
			$e = $l[3] + $l[8] + $r2;
		} elsif( $l[-1] =~ /XG:Z:GA/ ) {	## crick strand
			$strand = '-';
			$s = $l[3] - $r2;
			$e = $l[3] + $l[8] + $r1;
		}
		print OUT join("\t", $l[2], $s, $e, '.', $l[4], $strand), "\n";
	}
}
close IN;
close OUT;

