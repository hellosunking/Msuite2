#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.size> [in.size ...]\n\n";
	exit 2;
}

my %size;
foreach my $file ( @ARGV ) {
	open S, "$file" or die( "$!" );
	while( <S> ) {
		next if /^#/;	#Locus	C	T	Z
		chomp;
		my @l = split /\t/;
		$size{$l[0]} += $l[1];
	}
	close S;
}

print "Size\tCount\n";
my @cycle = sort {$a<=>$b} keys %size;
foreach my $i ( $cycle[0] .. $cycle[-1] ) {
	print join("\t", $i, $size{$i}||0), "\n";
}

