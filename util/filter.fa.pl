#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.fa[.gz]>\n";
	print STDERR "\nThis program is designed to discard unwanted scaffolds in the genome.";
	print STDERR "\nOnly the chromosomes named using digits, X, Y, and M/MT will be kept.\n\n";
	exit 2;
}

my $pass = 1;
open IN, "less $ARGV[0] | " or die( "$!" );
while( <IN> ) {
	if( /^>(\S+)/ ) {
		my $chr = $1;
		$chr =~ s/^chr//i;
		$chr =~ s/MT/M/;

		if( $chr=~/^\d+$/ || $chr=~/^[XYM]$/ ) {
			print ">chr$chr\n";
			$pass = 1;
		} else {
			print STDERR "Discard $chr\n";
			$pass = 0;
		}
	} else {
		print if $pass;
	}
}
close IN;


