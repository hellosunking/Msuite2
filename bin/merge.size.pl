#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <w|c> <in.size> [in.size ...]\n\n";
	exit 2;
}

my $strand = shift;

my (%auto_size, %lambda_size, %pUC19_size);
foreach my $file ( @ARGV ) {
	if( $file =~ /[cr]hr\d+/ ) {
		open S, "$file" or die( "$!" );
		while( <S> ) {
			next if /^#/;	#size count
			chomp;
			my @l = split /\t/;
			$auto_size{$l[0]} += $l[1];
		}
		close S;
	} elsif( $file =~ /[cr]hrL\./ ) {
		open S, "$file" or die( "$!" );
		while( <S> ) {
			next if /^#/;	#size count
			chomp;
			my @l = split /\t/;
			$lambda_size{$l[0]} += $l[1];
		}
		close S;
	} elsif( $file =~ /[cr]hrP\./ ) {
		open S, "$file" or die( "$!" );
		while( <S> ) {
			next if /^#/;	#size count
			chomp;
			my @l = split /\t/;
			$pUC19_size{$l[0]} += $l[1];
		}
		close S;
	}
}

open OUT, ">Msuite2.$strand.size" or die( "$!" );
print OUT "Size\tCount\n";
my @cycle = sort {$a<=>$b} keys %auto_size;
foreach my $i ( $cycle[0] .. $cycle[-1] ) {
	print OUT join("\t", $i, $auto_size{$i}||0), "\n";
}

open OUT, ">Msuite2.$strand.lambda.size" or die( "$!" );
print OUT "Size\tCount\n";
@cycle = sort {$a<=>$b} keys %lambda_size;
foreach my $i ( $cycle[0] .. $cycle[-1] ) {
	print OUT join("\t", $i, $lambda_size{$i}||0), "\n";
}

open OUT, ">Msuite2.$strand.pUC19.size" or die( "$!" );
print OUT "Size\tCount\n";
@cycle = sort {$a<=>$b} keys %pUC19_size;
foreach my $i ( $cycle[0] .. $cycle[-1] ) {
	print OUT join("\t", $i, $pUC19_size{$i}||0), "\n";
}

