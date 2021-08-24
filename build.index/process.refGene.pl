#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (hellosunnking@foxmail.com)
#

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.refGene.txt[.gz|bz2]>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

## prepare the [-4000, +4000] window sround TSS with 40 bins
my $step = 200;
my $bin  = 20;

my %tss;

if( $ARGV[0] =~ /\.bz2$/ ) {
	open IN, "bzip2 -cd $ARGV[0] |" or die( "$!" );
} elsif( $ARGV[0] =~ /\.gz$/ ) {
	open IN, "gzip  -cd $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##585	NR_024540	chr1	-	14361	29370	29370	29370	11	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,	0	WASH7P	unk	unk -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,

#	next unless $l[1] =~ /^NM_/;		## keep mRNAs only
	next unless $l[2] =~ /^chr[\dXY]+$/;	## discard those on chr6_random thing
	next if $l[12]=~/^MIR/ || $l[12]=~/^SNO/;	## discard miRNA and snoRNA
	next unless $l[5]-$l[4] >= 1000;		## discard genes that are too short
	if( $l[3] eq '+' ) {
		$tss{"$l[2]:$l[4]:+"} = $l[12];
	} else {
		$tss{"$l[2]:$l[5]:-"} = $l[12];
	}
}
close IN;

foreach my $site ( keys %tss ) {
	my ($chr, $locus, $strand) = split /:/, $site;
	my $gene = $tss{$site};

	if( $strand eq '+' ) {
		for(my $s=-$bin; $s<=$bin; ++$s ) {
			print join("\t", $chr, $locus+$s*$step, $locus+$s*$step+$step, $gene, $s*$step, '+'), "\n";
		}
	} else {
		for(my $s=-$bin; $s<=$bin; ++$s ) {
			print join("\t", $chr, $locus-$s*$step-$step, $locus-$s*$step, $gene, $s*$step, '-'), "\n";
		}
	}
}

