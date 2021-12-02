#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (hellosunnking@foxmail.com)
#

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <gene.anno.gff[.gz|bz2]>\n\n";
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
	next if /^#/;
	chomp;
	next unless /\S+/;
	my @l = split /\t/;
	## chr1	RefSeq	exon	11874	12227	.	+	.	transcript_id "NR_046018.2_1"; gene_id "DDX11L1"; gene_name "DDX11L1"
#	next unless $l[0] =~ /^chr[\dXY]+$/;	## discard those chr6_random thing; does not work for plants
	$l[0] = "chr$l[0]" unless $l[0]=~/^chr/i;

	$l[-1] =~ s/"//g;
	$l[-1] =~ s/;/; /g;
	$l[-1] =~ s/=/ /g;
	my $id = '';
	if( $l[-1] =~ /transcript_id (\S+);?/i ) {
		$id = $1;
	} else {
#		print STDERR "ERROR: there is NO transcript_id information on line $.!\n";
		next;
	}
	$id =~ s/;$//;

	if( $l[-1] =~ /gene_name (\S+);?/ ) {
		next if $1=~/^MIR/ || $1=~/^SNO/;	## discard miRNA and snoRNA
	}

	my $tid = "$id:$l[0]:$l[6]";
	if( $l[6] eq '+' ) {
		if( exists $tss{$tid} ) {
			$tss{$tid} = $l[3] if $l[3] < $tss{$tid};
		} else {
			$tss{$tid} = $l[3];
		}
	} else {
		if( exists $tss{$tid} ) {
			$tss{$tid} = $l[4] if $l[4] > $tss{$tid};
		} else {
			$tss{$tid} = $l[4];
		}
	}
}
close IN;

my %meet;
foreach my $site ( keys %tss ) {
	my ($id, $chr, $strand) = split /:/, $site;
	my $locus = $tss{$site};

	## for multiple trancripts that sharing 1 promoter, only output once
	next if exists $meet{"$chr:$locus:$strand"};
	$meet{"$chr:$locus:$strand"} = 1;

	if( $strand eq '+' ) {
		for(my $s=-$bin; $s<=$bin; ++$s ) {
			next if $locus+$s*$step < 0;
			print join("\t", $chr, $locus+$s*$step, $locus+$s*$step+$step, $id, $s*$step, '+'), "\n";
		}
	} else {
		for(my $s=-$bin; $s<=$bin; ++$s ) {
			next if $locus-$s*$step-$step < 0;
			print join("\t", $chr, $locus-$s*$step-$step, $locus-$s*$step, $id, $s*$step, '-'), "\n";
		}
	}
}

