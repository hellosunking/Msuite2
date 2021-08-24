#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : May 6, 2020

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Msuite.final.bam> [frag.size=200] [output.bed=stdout]\n";
	print STDERR "\nThis program is designed to translate the bam file (Single-End mode) to bed file.\n\n";
	exit 1;
}

my $fs = $ARGV[1] || 200;

my $o = $ARGV[2] || '/dev/stdout';
open OUT, ">$o" or die( "$!" );

open IN, "samtools view $ARGV[0] |" or die( "$!" );
while( <IN> ) {
	my @l = split /\t/;	##SRR1045842.157 0 chr3 164603272 255 36M * 0 0 CCTCATTTGT
	my $pos = $l[3];
	my $chr = $l[2];
	if( $l[1] & 0x10 ) {	## reverse strand
		## SAM format records the leftmost position so I need to calculate the rightmost
		## position of this read then deduce the left position of the fragment
		my $trueReadLen = extract_size_from_CIGAR( $l[5] );
		print OUT "$chr\t", $pos+$trueReadLen-1-$fs, "\t", $pos+$trueReadLen-1, "\t$l[0]\t$l[4]\t-\n";
	} else {	## forward strand
		print OUT "$chr\t", $pos-1, "\t", $pos-1+$fs, "\t$l[0]\t$l[4]\t+\n";
	}
}
close IN;

close OUT;

sub extract_size_from_CIGAR {
	my $cigar = shift;
	
	if( $cigar =~ /^(\d+)M$/ ) {	## there is NO indels
		return $1;
	} else { ## there are indels here
		$cigar =~ s/([MID])/$1:/g;
		my @info = split /:/, $cigar;
		my $size = 0;
		foreach my $m ( @info ) {
			if( $m =~ /^(\d+)I$/ ) {        ## insertion: CAUTION!!! I WILL DISCARD IT!!!
			} elsif ( $m =~ /^(\d+)D$/ )    { ## deletion: CAUTION!!! I WILL ADD 'D' to indicate deletion !!!
				$size += $1;
			} elsif( $m =~ /^(\d+)M$/ ) {
				$size += $1;
			}
		}
		return $size;
	}
}


