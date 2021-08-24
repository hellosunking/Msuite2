#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <in.sam> <use.watson> <use.crick>\n\n";
	exit 1;
}

my $watson = $ARGV[1];
my $crick  = $ARGV[2];

## PE/SE will be automatically detected
## load data
my %raw;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##555626_chr1:54339-54614_R1	99	chr1	54339	23	36M	=	54579	276	TCAATTAAGAGAAACCGTACCTATGCTATTTTGTCC	HHHHHHHHHHGGGGFFEEDDCBA@?>=<;986420.	XG:Z:GA
	$l[0] =~ /[#\s].*$/;

	my $strand;
	if( $l[-1] =~ /GA$/ ) {
		next unless $crick;
		$strand = 'C';
	} else {
		next unless $watson;
		$strand = 'W';
	}

	my $pos = $l[3];
	my ($mseq, $mqual);

	my $cigar = $l[5];
	if( $cigar !~ /^\d+M$/ ) {	## there are indels here
		$cigar =~ s/([MID])/$1:/g;
		my @info = split /:/, $cigar;
		($mseq, $mqual) = ('', '');	## modified seq and qual
		my $curr = 0;
		foreach my $m ( @info ) {
			if( $m =~ /^(\d+)I$/ ) {	## insertion: CAUTION!!! I WILL DISCARD IT!!!
				$curr += $1;
			} elsif ( $m =~ /^(\d+)D$/ )	{ ## deletion: CAUTION!!! I WILL ADD 'D' to indicate deletion !!!
				$mseq .= 'D' x $1;
				$mqual.= '!' x $1;
			} elsif( $m =~ /^(\d+)M$/ ) {
				$mseq .= substr($l[9],  $curr, $1);
				$mqual.= substr($l[10], $curr, $1);
				$curr += $1;
			}
		}
	} else {
		$mseq .=  $l[9];
		$mqual.= $l[10];
	}

	push @{$raw{$l[0]}}, "$pos\t$mseq\t$mqual\t$strand";
}
close IN;

## process data, deal with PE/SE
my ($p1, $s1, $q1, $strand);
my ($p2, $s2, $q2);
foreach my $k ( keys %raw ) {
	my $sam = $raw{$k};
	if( $#$sam == 0 ) {	## SE or only 1 read is in the given region
		my ($pos, $mseq, $mqual, $strand) = split /\t/, $sam->[0];
		print "$pos\t$mseq\t$strand\n";
	} elsif ( $#$sam == 1 ) {	## PE
		my @info1 = split /\t/, $sam->[0];
		my @info2 = split /\t/, $sam->[1];

		if( $info1[0] <= $info2[0] ) {
			($p1, $s1, $q1, $strand) = @info1;
			($p2, $s2, $q2, $strand) = @info2;
		} else {
			($p1, $s1, $q1, $strand) = @info2;
			($p2, $s2, $q2, $strand) = @info1;
		}

		if( $p1 + length($s1) <= $p2 ) {   # no overlap between read1 and read2
			my $merge = $s1;
			$merge .= '.' x ($p2-$p1-length($s1));
			$merge .= $s2;

			print "$p1\t$merge\t$strand\n";
		} else {    # merge read1 and read2
			my $merge = '';
			if( $p1 + length($s1) < $p2+length($s2) ) {
				#my $frag_len = $p2 - $p1 + length($s2);
				my $k;
				#for( $k=0; $k!=$p2-$p1; ++$k ) { # read1 only
				#	$merge .= substr($s1,$k,1);  # can also use substr
				#}
				$merge = substr($s1, 0, $p2-$p1);
				for( $k=$p2-$p1; $k!=length($s1); ++$k ) {   # overlapped region
					if( substr($q1,$k,1) ge substr($q2,$k+$p1-$p2,1) ) {
						$merge .= substr($s1, $k, 1);
					} else {
						$merge .= substr($s2, $k+$p1-$p2, 1);
					}
				}
				#for( ; $k!=$frag_len; ++$k ) {   # read2 only
				#	$merge .= substr($s2, $k+$p1-$p2, 1);
				#}
				$merge .= substr($s2, length($s1)+$p1-$p2);
			} else {	## R1 completely contains R2, very rare
				$merge = $s1;
			}
#			if( $strand eq 'C' ) {
#				$merge = reverse $merge;
#				$merge =~ tr/ACGT/TGCA/;
#			}
			print "$p1\t$merge\t$strand\n";
		}
	}
}

