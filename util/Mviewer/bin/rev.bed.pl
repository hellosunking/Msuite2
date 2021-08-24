#!/usr/bin/perl
use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <genome.info> <query.bed>\n\n";
	exit 2;
}

my %chr2len;
open CHRLEN, "$ARGV[0]" or die( "$!" );
while( <CHRLEN> )
{
	chomp;
	my ($chr, $len) = split /\s+/;
	$chr2len{$chr} = $len + 1;	## for easier computing
}
close CHRLEN;

open BED, "$ARGV[1]" or die "$!";
while( <BED> )
{
	chomp;
	my ( $chr, $start, $end, $extra ) = split /\t/, $_, 4;
	next unless exists $chr2len{$chr};
	my $rs = $chr2len{$chr} - $end;
	my $re = $chr2len{$chr} - $start;
	print "$chr\t$rs\t$re\t$extra\n";
}
close BED;

