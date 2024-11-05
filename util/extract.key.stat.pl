#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <Msuite2.dir> [Msuite2.dir ...]\n\n";
	exit 2;
}

print "#Sid\tAll\tReported\t%Usable\t%DNAm\t%ConversionRate\n";
foreach my $sid ( @ARGV ) {
	( -s "$sid/Msuite2.report/index.html" ) || next;
	open IN, "$sid/Msuite2.report/index.html" or die( "$!" );

	my ($raw, $clean, $meth, $conversion) = (0, 0, 0, 0);
	while( <IN> ) {
		chomp;
		if( /Total input reads/ ) {
			if( /<b>([,\d]+)<\/b>/) {
				$raw = $1;
			}
		} elsif( /Reported alignments/ ) {
			if( /<b>([,\d]+)/) {
				$clean = $1;
			}
		} elsif( /Overall CpG methylation density/ ) {
			if( /<b>([\.\d]+)/) {
				$meth = $1;
			}
		} elsif( /conversion rate/ && /lambda/ ) {
			if( /<b>([\.\d]+)/) {
				$conversion = $1;
			}
		}
	}
	close IN;

	my $r = $raw;
	my $c = $clean;
	$r =~ s/,//g;
	$c =~ s/,//g;
	my $ratio = sprintf( "%.2f", $c/$r*100 );

	print join("\t", $sid, $raw, $clean, $ratio, $meth, $conversion), "\n";
}

