#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <output.prefix> <in.call> [in.call ...]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $output = shift;

my %call;
foreach my $file ( @ARGV ) {
	print "Loading $file\n";
	open IN, "$file" or die( "$!" );
	<IN>;	## skip header
	while( <IN> ) {
		my @l = split /\t/, $_, 3;	##chr	Locus	Total	wC	wT	wOther	Context	cC	cT	cOther
		$call{$l[0]}->{$l[1]} .= $l[2];
	}
	close IN;
}

open CALL, ">$output.meth.call"     or die( "$!" );
open BED,  ">$output.meth.bedgraph" or die( "$!" );

foreach my $chr ( sort keys %call ) {
	my $here = $call{$chr};
	foreach my $locus ( sort {$a<=>$b} keys %$here ) {
		my @info = split /\n/, $here->{$locus};
		my ($total, $wC, $wT, $wO, $cC, $cT, $cO) = (0,0,0,0,0,0,0);
		my $context = "NA";
		foreach my $record ( @info ) {
			my @l = split /\t/, $record;	##Total   wC  wT  wOther  Context cC  cT  cOther
			$total += $l[0];
			$wC    += $l[1];
			$wT    += $l[2];
			$wO    += $l[3];
			$cC    += $l[5];
			$cT    += $l[6];
			$cO    += $l[7];
			$context = $l[4];
		}
		print CALL join("\t", $chr, $locus, $total, $wC, $wT, $wO, $context, $cC, $cT, $cO), "\n";

		if( $wC+$wT+$cC+$cT == 0 ) {
#			print STDERR "WARNING: $chr $locus\n";
		} else {
			my $m = ($wC+$cC)/($wC+$wT+$cC+$cT)*100;
			print BED join("\t", $chr, $locus-1, $locus, $m), "\n";
		}
	}
}
close CALL;
close BED;

