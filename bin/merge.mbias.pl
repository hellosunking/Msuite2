#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.bias> [in.bias ...]\n\n";
	exit 2;
}

my %bias;
foreach my $file ( @ARGV ) {
	open MB, "$file" or die( "$!" );
	while( <MB> ) {
		next if /^#/;	#Locus	C	T	Z
		chomp;
		my @l = split /\t/;
		$bias{$l[0]}->{C} += $l[1];
		$bias{$l[0]}->{T} += $l[2];
	}
	close MB;
}

print "Cycle\tC\tT\n";
my @cycle = sort {$a<=>$b} keys %bias;
foreach my $i ( $cycle[0] .. $cycle[-1] ) {
	print join("\t", $i, $bias{$i}->{C}||0, $bias{$i}->{T}||0), "\n";
}

