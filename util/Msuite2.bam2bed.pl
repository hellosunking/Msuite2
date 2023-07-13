#!/usr/bin/perl
#
# Author : Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Version: Oct 2022
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <Msuite2.final.bam> <output.bed[.gz]> <out.size> [R1.head.cut=0] [R2.head.cut=0] [min.qual=0] [autosomal.only=No|Yes]\n\n",
				 "This program is designed to translate the Msuite2's output BAM file to bed format.\n\n";
	exit 2;
}

my $R1headcut = $ARGV[3] || 0;
my $R2headcut = $ARGV[4] || 0;
my $minqual   = $ARGV[5] || 0;
my $autoONLY  = $ARGV[6] || 0;

my $fixHead = $R1headcut + $R2headcut;

if( $ARGV[0] =~ /bam$/ ) {
	open IN, "samtools view -@ 2 $ARGV[0] |" or die("$!");
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}

my $outbed = $ARGV[1];
if( $outbed =~ /\.gz$/ ) {
	open OUT, "| gzip >$outbed" or die( "$!" );
} else {
	open OUT, ">$outbed" or die( "$!" );
}

my %size;
my %chrcount;
my ($chr, $pos, $strand);
my @l;
my $all = 0;
while( <IN> ) {
	@l = split /\t/;
	my $mate = $l[8];	## QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL EXTRA-MARKS
	next if $mate<=0 || $l[4]<$minqual;
	next if $l[2]=~/^chrM/;	## always discard chrM; BUT keep HBV, Lambda, EBV/E.coli
	next if $autoONLY && $l[2]!~/^chr\d+$/;

	## For BWA result
#	next unless $l[1] & 0x02;       ## the fragment is propoerly mapped to the reference genome
#	next if     $l[1] & 0x100;      ## discard secondary alignment
#	next unless $l[4] >= $minqual;  ## high quality mapping
#	next unless $l[6] eq '=';		## properly mapped pair
#	my $mate = $l[8];
#	next if $mate <= 0;

	++ $chrcount{$l[2]};

	$chr = $l[2];
	$pos = $l[3]-1;
	$strand = ($l[1] & 0x40) ? '+' : '-';
	## for the left-most read:
	## if it is the first template, then it is read1 and the fragment should be on watson chain
	## otherwise it is read2 and the fragment should be on crick chain
	## for sorted BAM files, you cannot ensure the order of READ1 and READ2,
	## therefore need to check this flag while NOT 0x10 for strand
	if( $strand eq '+' ) {
		if( $l[-1] =~ /CT$/ ) {
			print OUT "$chr\t", $pos-$R1headcut, "\t", $pos+$mate+$R2headcut, "\t$l[4]\t+\n";
		} else {
			print STDERR "ERROR $l[0]!\n";
			next;
		}
	} else {
		if( $l[-1] =~ /GA$/ ) {
			print OUT "$chr\t", $pos-$R2headcut, "\t", $pos+$mate+$R1headcut, "\t$l[4]\t-\n";
		} else {
			print STDERR "ERROR $l[0]!\n";
			next;
		}
	}

	if( $chr =~ /\d$/ ) {	## size pattern: always autosome only
		$mate += $fixHead;
		$size{$mate} ++;
		++ $all;
	}
}
close IN;
close OUT;

if( $all == 0 ) {
	print STDERR "ERROR: No valid reads loaded!\n";
	exit 1;
}

## size distribution
my $outsize = $ARGV[2];
exit 0 if $outsize eq "/dev/null";
my @fraglen = sort {$a<=>$b} keys %size;
my $maxLen = $fraglen[-1];

open OUT, ">$outsize" or die( "$!" );
print OUT "#Size\tCount\tPercent%\tCumulative\n";
my $cumu = 0;
for( my $i=1; $i<=$maxLen; ++$i ) {
	my $here = $size{$i} || 0;
	$cumu += $here;
	print OUT join("\t", $i, $here, $here/$all*100, $cumu/$all), "\n";
}
close OUT;

