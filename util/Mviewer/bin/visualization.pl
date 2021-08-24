#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;
use FindBin qw/$Bin/;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <genome> <in.region> <in.file> [mode=TAPS|BS]\n\n";
	exit 2;
}

## m/u signals in the data
my ($wmSig, $wuSig, $cmSig, $cuSig);
if( defined $ARGV[3] && $ARGV[3]=~/BS/i ) {	## configuration for BS-seq
	$wmSig = 'C';
	$wuSig = 'T';
	$cmSig = 'G';
	$cuSig = 'A';
} else {	## configuration for TAPS
	$wmSig = 'T';
	$wuSig = 'C';
	$cmSig = 'A';
	$cuSig = 'G';
}

## m/u signals to be shown
my ($mTag, $uTag) = ('M', 'U');

##################################
my $g_file = $ARGV[0];
my $region = $ARGV[1];
my $ifile  = $ARGV[2];

## process region
my ($chr, $spos, $epos) = split /[:-]/, $region;

## load genome
my $fa = 'X';
open IN, "$g_file" or die( "$!" );
while( <IN> ) {
	last if s/^>$chr\s//;	## locate the target chr
}
while( <IN> ) {
	chomp;
	last if s/^>//;	## another chr meet, stop
	$fa .= uc $_;
}
close IN;

if( $fa eq "" ) {
	print STDERR "Error: Could not find $chr in the genome!\n";
	exit 1;
}

## load input file
open IN, "$ifile" or die( "$!" );
my %dat;
my ($lmost, $rmost) = (3e9, 0);
while( <IN> ) {
	chomp;
	my @l = split /\t/;	## pos seq strand
	$lmost = $l[0] if $l[0] < $lmost;
	$rmost = $l[0] if $l[0] > $rmost;

	push @{$dat{$l[0]}}, "$l[1]\t$l[2]";
}
close IN;

## plot the data, watson and crick reads are separated
## plot watson first
#print substr( $g->{$chr}, $pos[0], $pos[-1]-$pos[0]+150 ), "\n";
my $max = 0;
my $vis = '';
my @pos = sort {$a<=>$b} keys %dat;
foreach my $p ( @pos ) {
	#print STDERR "Dealing $p\n";
	my $h = $dat{$p};
	foreach my $info ( @$h ) {
		my ($r, $strand) = split /\t/, $info;
		#print STDERR "Dealing $r (", length($r), ")\n";

		my $visual = ' ' x ($p-$lmost);
		for( my $i=0; $i<length($r); ++$i ) {
			my $ref = substr( $fa, $p+$i, 1 );
			my $seq = substr(  $r,    $i, 1 );
			#print STDERR "$i\t$ref\t$seq\n";
			if( $ref eq 'C' ) {
				if( substr( $fa, $p+$i+1, 1) eq 'G' ) {	## CpG site
					if( $strand eq 'W' ) {	## watson chain, straightforward
						if( $seq eq $wmSig )	{ $visual .= $mTag; }
						elsif( $seq eq $wuSig )	{ $visual .= $uTag; }
						else					{ $visual .= $seq;  }
					} else {	## crick chain, check the G site for Meth signal while C for mutation
						++ $i;
						if( $i < length($r) ) {	## in case the C is at the END of the sequence
							my $meth = substr( $r, $i, 1 );
							#$meth =~ tr/AG/TC/;
							if( $meth eq $cmSig )	{ $visual .= $mTag; }
							elsif( $meth eq $cuSig ){ $visual .= $uTag; }
							else					{ $visual .= $meth; }

							## check the mutation signal for the C site, BUT this signal will be shown on G site
							## to keep the M signal well synchronized in the visualization
							if( $seq eq 'C' )	{ $visual .=  '-'; }
							else				{ $visual .= $seq; }
						} else {
							## discard this site
						}
					}
				} else {
					if( $seq eq $wmSig || $seq eq $wuSig ) {
						$visual .= '-';
					} else {
						$visual .= $seq;
					}
				}
			} elsif( $ref eq 'G' && $strand eq 'C' ) {	## G->A issue in crick chain
				if( $seq eq $cmSig || $seq eq $cuSig ) {
					$visual .=  '-';
				} else {
					$visual .= $seq;
				}
			} else {
				if( $seq eq $ref ) {
					$visual .=  '-';
				} else {
					$visual .= $seq;
				}
			}
		}
		$vis .= "$visual\n";
		$max = length($visual) if length ($visual) > $max;
	}
}

print "#region=$region\n",
	  "#position=$chr:$lmost-", $lmost+$max, "\n",
	  substr( $fa, $lmost, $max+2 ), "\n", $vis;

