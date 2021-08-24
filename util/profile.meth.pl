#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : May 6, 2020

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <Msuite.meth.call> <bin.size> [chr.mask=chr] [protocol=BS|TAPS]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $bin = $ARGV[1];
my $msk = $ARGV[2] || 'chr';
my $pro = $ARGV[3] || 'BS';

my (%C, %T);
open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## header
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr    Locus   Total   wC      wT      wOther  Context cC      cT      cOther
	my $index = int($l[1]/$bin);
	$C{$l[0]}->{$index} += $l[3] + $l[7];
	$T{$l[0]}->{$index} += $l[4] + $l[8];
}
close IN;

foreach my $chr ( sort keys %C ) {
	next unless $chr=~/[\dXY]$/ && $chr!~/_/;
	my $masked_chr = $chr;
	$masked_chr =~ s/^chr//;
	$masked_chr = "$msk$masked_chr";

	my @loci = sort {$a<=>$b} keys %{$C{$chr}};

	foreach my $i ( @loci ) {
		my $sigC = $C{$chr}->{$i};
		my $sigT = $T{$chr}->{$i};
		next if $sigC+$sigT==0;
		my $m;
		if( $pro =~ /BS/i ) {
			$m = $sigC/($sigC+$sigT)*100;
		} else {
			$m = $sigT/($sigC+$sigT)*100;
		}
		print join("\t", $masked_chr, $i*$bin, $i*$bin+$bin, $m), "\n";
	}
}

