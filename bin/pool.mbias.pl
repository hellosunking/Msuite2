#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;
#use KSLIB::cigarUtil qw/fix_seq_from_CIGAR extract_size_from_CIGAR/;
#use KSLIB::Digitalize qw/digitalize/;

if( $#ARGV < 3 ) {
	print STDERR "\nUsage: $0 <chr.info> <dir> <cycle> <mode=BS|TAPS>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $bs;
if( $ARGV[3]=~/^BS$/i ) {
	$bs = 1;
} elsif( $ARGV[3]=~/^TAPS$/i ) {
	$bs = 0;
} else {
	print STDERR "ERROR: unsupported mode!\n";
	exit 1;
}

my (%R1w, %R2w, %R1c, %R2c);

open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	my $chr = $l[0];

	load_mbias( "$ARGV[1]/chr$chr.R1.mbias", \%R1w );
	load_mbias( "$ARGV[1]/chr$chr.R2.mbias", \%R2w );
	load_mbias( "$ARGV[1]/rhr$chr.R1.mbias", \%R1c );
	load_mbias( "$ARGV[1]/rhr$chr.R2.mbias", \%R2c );
}
close IN;

print "#Cycle\tR1w\tR2w\tR1c\tR2c\n";
foreach my $i ( 1 .. $ARGV[2] ) {
	my ($m1w, $m2w, $m1c, $m2c);
	if( $bs ) {
		$m1w = $R1w{$i}->{C} / ( $R1w{$i}->{C} + $R1w{$i}->{T} ) * 100;
		$m2w = $R2w{$i}->{C} / ( $R2w{$i}->{C} + $R2w{$i}->{T} ) * 100;
		$m1c = $R1c{$i}->{C} / ( $R1c{$i}->{C} + $R1c{$i}->{T} ) * 100;
		$m2c = $R2c{$i}->{C} / ( $R2c{$i}->{C} + $R2c{$i}->{T} ) * 100;
	} else {
		$m1w = $R1w{$i}->{T} / ( $R1w{$i}->{C} + $R1w{$i}->{T} ) * 100;
		$m2w = $R2w{$i}->{T} / ( $R2w{$i}->{C} + $R2w{$i}->{T} ) * 100;
		$m1c = $R1c{$i}->{T} / ( $R1c{$i}->{C} + $R1c{$i}->{T} ) * 100;
		$m2c = $R2c{$i}->{T} / ( $R2c{$i}->{C} + $R2c{$i}->{T} ) * 100;
	}
	print join( "\t", $i, $m1w, $m2w, $m1c, $m2c ), "\n";
}

sub load_mbias {
	my $file = shift;
	my $bias = shift;

	open MB, "$file" or die( "$!" );
	while( <MB> ) {
		next if /^#/;	#Locus	C	T	Z
		chomp;
		my @l = split /\t/;
		$bias->{$l[0]}->{C} = $l[1];
		$bias->{$l[0]}->{T} = $l[2];
	}
	close MB;
}

