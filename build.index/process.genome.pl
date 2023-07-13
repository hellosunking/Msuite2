#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (hellosunking@foxmail.com)
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <in.genome> <lambda.genome[,pUC.genome]> <out.dir>\n";
	print STDERR "\nThis program is designed to treat CpG sites.";
	print STDERR "\nPlease note that chrL and chrP are reserved, and they are NOT allowed to present in your genome.\n\n";
	exit 2;
}

my @falist;

if( -d "$ARGV[0]" ) {	## parameter 1 is a directory, then load all FASTA files under this directory
	print "INFO: The given parameter is a directory! I will load all the fasta files in it!\n";
	opendir DIR, "$ARGV[0]" or die( "$!" );
	while( my $fa = readdir(DIR) ) {
		push @falist, "$ARGV[0]/$fa" if $fa =~ /\.fa(sta)?$/;
	}
	closedir DIR;
} else {	## parameter 1 is a file, which I assume should contain all sequences
	print "INFO: The given parameter is a file and I will load it directly!\n";
	push @falist, $ARGV[0];
}

my %g;
my $chr;
print "Loading genome files ...\n";
foreach my $fa ( @falist ) {
	open IN, "less $fa |" or die( "$!" );

	while( <IN> ) {
		chomp;
		if( s/^>(\S+)// ) {
			$chr = $1;
			$chr =~ s/^chr//i;
			if( $chr eq 'L' || $chr eq 'P' ) {	## chrL is reserved for Lambda genome
				print STDERR "ERORR: Your genome contains chrL or chrP, which is NOT permitted! Please rename these chromosomes and try again.\n";
				close IN;
				exit 100;
			}
			$g{$chr} = '';
		} else {
			$g{$chr} .= uc $_;
		}
	}
	close IN;
}

## add lambda,pUC19 genome
my @spikeinfiles = split /,/, $ARGV[1];
my $lambda_fa = $spikeinfiles[0];
print "Loading lambda genome: $lambda_fa\n";
open IN, "less $lambda_fa |" or die( "$!" );
<IN>;	# chrL header
$g{'L'} = '';
while( <IN> ) {
	chomp;
	if( s/^>// ) {
		print STDERR "ERROR: Your lambda genome file is incorrect, please check it.\n";
		exit 10;
	} else {
		$g{'L'} .= uc $_;
	}
}
close IN;

if( $#spikeinfiles > 0 ) {	## pUC19 is provided
	my $pUC19_fa = $spikeinfiles[1];
	print "Loading pUC19 genome: $pUC19_fa\n";
	open IN, "less $pUC19_fa |" or die( "$!" );
	<IN>;	# header
	$g{'P'} = '';
	while( <IN> ) {
		chomp;
		if( s/^>// ) {
			print STDERR "ERROR: Your pUC19 genome file is incorrect, please check it.\n";
			exit 11;
		} else {
			$g{'P'} .= uc $_;
		}
	}
	close IN;
}

print "Processing genome files ...\n";
my $bp_per_len = 50;
open CG2TG, ">$ARGV[2]/CG2TG.fa"  or die( "$!" );
open C2T,   ">$ARGV[2]/C2T.fa"    or die( "$!" );

open SIZE,  ">$ARGV[2]/chr.info"  or die( "$!" );
#open ORIW,  ">$ARGV[2]/watson.fa" or die( "$!" );
#open ORIC,  ">$ARGV[2]/crick.fa"  or die( "$!" );

my ($cnt1, $cnt2, $cnt3, $cnt4) = ( 0, 0, 0, 0 );
foreach $chr ( sort keys %g ) {
	my $watson = $g{$chr};
	my $crick  = reverse $watson;
	$crick =~ tr/ACGT/TGCA/;
	print SIZE "$chr\t", length($watson), "\n";

	## record the raw sequence and chr size information
	my $i = 0;
	open ORIW, ">$ARGV[2]/fasta/w$chr.fa" or die( "$!" );
	print ORIW ">chr$chr\n";
	while( 1 ) {
		my $s = substr( $watson, $i, $bp_per_len );
		my $len = length($s);
		print ORIW "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}
	close ORIW;

	$i = 0;
	open ORIC, ">$ARGV[2]/fasta/c$chr.fa" or die( "$!" );
	print ORIC ">rhr$chr\n";
	while( 1 ) {
		my $s = substr( $crick, $i, $bp_per_len );
		my $len = length($s);
		print ORIC "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}
	close ORIC;

	## CG->TG
	my $seq = $watson;
	$cnt1 += ($seq =~ s/CG/TG/g);
	print CG2TG ">chr$chr\n";
	$i = 0;
	while( 1 ) {
		my $s = substr( $seq, $i, $bp_per_len );
		my $len = length($s);
		print CG2TG "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	$seq = $crick;
	$cnt2 += ($seq =~ s/CG/TG/g);
	print CG2TG ">rhr$chr\n";
	$i = 0;
	while( 1 ) {
		my $s = substr( $seq, $i, $bp_per_len );
		my $len = length($s);
		print CG2TG "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## C->T
	$seq = $watson;
	$cnt3 += ($seq =~ s/C/T/g);
	print C2T ">chr$chr\n";
	$i = 0;
	while( 1 ) {
		my $s = substr( $seq, $i, $bp_per_len );
		my $len = length($s);
		print C2T "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## G->A
	$seq = $crick;
	$cnt4 += ($seq =~ s/C/T/g);
	print C2T ">rhr$chr\n";
	$i = 0;
	while( 1 ) {
		my $s = substr( $seq, $i, $bp_per_len );
		my $len = length($s);
		print C2T "$s\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}
}

print "Done. Conversions:\n",
	  "CG->TG: $cnt1\n",
	  "CG->CA: $cnt2\n",
	  "C -> T: $cnt3\n",
	  "G -> A: $cnt4\n";

close CG2TG;
close C2T;
close SIZE;

