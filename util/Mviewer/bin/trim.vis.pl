#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <in.vis> <targetLoci>\n\n";
	exit 2;
}

my ($targetChr, $targetStart, $targetEnd) = split /[:-]/, $ARGV[1];

open IN, "$ARGV[0]" or die( "$!" );
my $id = <IN>;	## gene=XXX
my $pos = <IN>;	## position=chrXX:XXX-XXX
$pos =~ /(chr\S+):(\d+)-(\d+)/;
my ( $chr, $start, $end ) = ( $1, $2, $3 );
my $seq = <IN>;	## genomic sequence
chomp($seq);
$end = $start + length($seq)-1;	## fix the real end position

if( $targetChr ne $chr || $start>$targetEnd || $end<$targetStart )
{
	print STDERR "Error: target loci cannot be implemented!\n";
	exit 100;
}

my $realStart = ($targetStart>$start) ? $targetStart : $start;
my $realEnd   = ($targetEnd<$end)     ? $targetEnd   : $end;

my $trimPos = $realStart - $start;
my $trimLen = $realEnd - $realStart + 1;


print "$id#position=$targetChr:$realStart-$realEnd\n",
	  substr($seq, $trimPos, $trimLen), "\n";

while( <IN> )
{
	chomp;
	next if length($_) <= $trimPos;
	my $here = substr($_, $trimPos, $trimLen);
	print "$here\n" unless $here=~/^[\s=-]+$/;
}
close IN;

