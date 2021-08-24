#!/usr/bin/perl

#
# Author: Ahfyth
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <in.vis> <level=1|2>\n\n";
	exit 1;
}

open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## gene=xxx
my $info = <IN>;	## pos=chrX:xx-xx
print $info;
my @l = split /[=:-]/, $info;
#print STDERR "Chr=$l[1];Start=$l[2];End=$l[3]\n";
my $ref = <IN>;		## reference sequence

my $cpg = get_CpG_site( $ref );
mark( $l[2], $cpg );
my $deal_ref = deal( $ref, $cpg );	## ref cannot be void
print "$deal_ref\n";

## light dense: only keep CpG sites, one read per line
if( $ARGV[1] == 1 ) {
	while( <IN> ) {
		chomp;
		my $result = deal( $_, $cpg );
		print "$result\n" if $result;
	}
	close IN;

	exit 0;
}

## heavy dense: only keep CpG sites, multiple reads per line
my @seq;
while( <IN> ) {
	chomp;
	my $result = deal( $_, $cpg );
	push @seq, $result if $result;
}
close IN;

my $DIST = 5;

my %drawed;
for( my $i=0; $i<=$#seq; ++$i )
{
	next if exists $drawed{$i};
	my $here = $seq[$i];
	for( my $j=$i+1; $j<=$#seq; ++$j )
	{
		next if exists $drawed{$j};
		my $s;
		for( $s=0; $s<length($seq[$j]); ++$s )
		{
			last if substr( $seq[$j], $s, 1 ) ne ' ';
		}
		if( $s >= length($here) + $DIST )	## far enough from the stop position, add it
		{
			$here .= substr( $seq[$j], length($here) );
			$drawed{$j} = 1;
		}
	}
	chomp($here);
	print "$here\n" if $here=~/\S+/;
}


### sub routines ###
sub mark
{
	my $offset = shift;
	my $cpg = shift;

	my $i = 0;
	foreach my $locus ( @$cpg ) {
		if( $i % 10 == 0 ) {
			my $here = $offset + $locus;
			print $here, ' ' x (20-length($here));
		}
		++ $i;
	}
	print "\n";
}

sub deal {
	my $fa  = shift;
	my $cpg = shift;

	my $here = '';
	foreach my $i ( @$cpg ) {
		last if $i > length($fa)-2;
		$here .= substr($fa, $i, 2);
	}
	$here =~ s/-$//;
	if( $here =~ /\S+/ ) {
		return $here;
	} else {
		return '';
	}
}

sub get_CpG_site
{
	my $fa = uc shift;
	my $i=0;
	my @CpG;
	while( ($i=index($fa, 'CG', $i))!=-1 ) {
		push @CpG, $i;
		++ $i;
	}

	return \@CpG;
}

