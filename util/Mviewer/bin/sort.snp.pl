#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;
use List::Util qw/sum/;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.vis>\n";
	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

open IN, "$ARGV[0]" or die( "$!" );
my $header = <IN>;
print $header;
$header = <IN>;
print $header;
$header = <IN>;
print $header;

my @seq;
my @geno;
while( <IN> ) {
	chomp;
	push @seq, $_;
	foreach my $i ( 0..length($_) )
	{
		my $here = substr($_, $i, 1);
		++ $geno[$i]->{$here} unless $here=~/[\. UM]/;
	}
}
close IN;

my $min_depth = 10;	## minimum depth to call SNP
my $min_snp   = 3;	## minimum coverage for the SNP allele

my $index;
for($index=0; $index<=$#geno; ++$index ) {
	next unless exists $geno[$index];
	my @here = sort {$b<=>$a} values %{$geno[$index]};
#	print STDERR "$index: $#here\n";
	next if $#here <= 0;	## only one allele
	my $total = sum @here;
	last if $total>=$min_depth && $here[1]>=$min_snp;
}

if( $index > $#geno ) {	## no SNP found
#	print STDERR "No SNP found!\n";
	print join("\n", @seq),"\n";
} else {
#	print STDERR "SNP index: $index!\n";
	my $uncovered = '';
	my %sorted;
	foreach my $s (@seq) {
		if( length($s)<=$index ) {
			$uncovered .= "$s\n";
		} else {
			$sorted{ substr($s,$index,1) } .= "$s\n";
		}
	}
	my @geno = sort keys %sorted;
	print $uncovered, join("", map { $sorted{$_} } @geno), "\n";
}

