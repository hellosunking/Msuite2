#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : Aug 2021

package MsuiteUtil;

use strict;
use warnings;
use Exporter 'import';
our @EXPORT = qw/$version $ver usage makefile_perchr makefile_methcall mk_samheader detect_cycle check_dependency check_index printRed printGrn printYlw makefile_perchr_v2/;

# long version
our $version = 'v2.0.1 (Dec 2021)';
# short version
our $ver = 'v2.0.1';

sub mk_samheader {
	my $chrinfo   = shift;
	my $genome    = shift;
	my $protocol  = shift;
	my $alignmode = shift;
	my $reads     = shift;
	my $samheader = shift;

	open OUT, ">$samheader" or die( "$!" );

	print OUT "\@HD\tVN:1.0\tSO:coordinate\n";
	open IN, "$chrinfo" or die( "$!" );
	while( <IN> ) {
		chomp;
		my ($chr, $size) = split /\t/;	## chr size
		$chr = "chr$chr";# unless $chr=~/^NC/;
		## do not add chr-prefix for NCBI-style chromosome?
		## need to be compatible to the SAM records
		print OUT "\@SQ\tSN:$chr\tLN:$size\n";
	}
	close IN;

	print OUT "\@PG\tID:Msuite2\tPN:Msuite2\tVN:$ver\tDS:genome=$genome;protocol=$protocol;mode=$alignmode;reads=$reads\n";
	close OUT;
}

sub makefile_perchr {
	my $MsuiteBin = shift;
	my $samtools  = shift;
	my $chrinfo   = shift;
	my $samheader = shift;
	my $seqMode   = shift;	## se or pe
	my $makefile  = shift;
	my $maxins    = shift || 1000;
	my $THREAD    = shift || 1;

	my $job = "";
	my $mkf = "";

	open IN, "$chrinfo" or die( "$!" );
	while( <IN> ) {
		chomp;
		my ($C, $size) = split /\t/;	## chr size
		my $chr = "chr$C";
		$job .= " $chr.srt.bam";
		$mkf .= "$chr.srt.bam: chr$C.sam rhr$C.sam\n";
		$mkf .= "\t\@$MsuiteBin/rmdup.w.$seqMode $maxins chr$C.sam chr$C >chr$C.rmdup.log\n";
		$mkf .= "\t\@$MsuiteBin/rmdup.c.$seqMode $size $maxins rhr$C.sam rhr$C >rhr$C.rmdup.log\n";
		$mkf .= "\t\@cat $samheader chr$C.rmdup.sam rhr$C.c2w.sam | samtools view -bS - | samtools sort -o chr$C.srt.bam -\n\n";
	}
	close IN;

	open MK, ">$makefile" or die( "$!" );

	print MK
"Msuite2.bai.OK: Msuite2.bam.OK
	$samtools index -\@ $THREAD ../Msuite2.final.bam
	\@cat *rmdup.log > ../Msuite2.rmdup.log
	\@touch Msuite2.bai.OK

Msuite2.bam.OK:$job
	$samtools merge -\@ $THREAD ../Msuite2.final.bam chr*.srt.bam
	\@touch Msuite2.bam.OK

$mkf";

	close MK;
}


sub makefile_methcall {
	my $MsuiteBin = shift;
	my $chrinfo   = shift;
	my $fastaDIR  = shift;
	my $seqMode   = shift;	## se or pe
	my $protocol  = shift;
	my $cycle     = shift;
	my $makefile  = shift;
	my $target    = shift || 'CpG';

	my $job = "";
	my $mkf = "";

	open IN, "$chrinfo" or die( "$!" );
	while( <IN> ) {
		chomp;
		my ($C, $size) = split /\t/;	## chr size
		my $chr = "chr$C";
		$job .= " $chr.meth.log";
		$mkf .= "$chr.meth.log: chr$C.rmdup.sam rhr$C.rmdup.sam\n";
		$mkf .= "\t\@$MsuiteBin/meth.caller.$target $seqMode $fastaDIR/w$C.fa chr$C.rmdup.sam $cycle chr$C\n";
		$mkf .= "\t\@$MsuiteBin/meth.caller.$target $seqMode $fastaDIR/c$C.fa rhr$C.rmdup.sam $cycle rhr$C\n";
		$mkf .= "\t\@$MsuiteBin/pair.$target chr$C $fastaDIR/w$C.fa $protocol chr$C.$target.call rhr$C.$target.call >chr$C.$target.meth.log\n\n";
	}
	close IN;

	open MK, ">$makefile" or die( "$!" );

	print MK
"Msuite2.$target.meth.call.OK:$job
	cat $MsuiteBin/meth.header chr*.$target.meth > ../Msuite2.$target.meth.call
	\@cat chr*.$target.meth.bedgraph > ../Msuite2.$target.meth.bedgraph
	\@cat chr*.$target.meth.log > ../Msuite2.$target.meth.log
	\@touch Msuite2.$target.meth.call.OK

$mkf";

	close MK;
}

sub printRed {
	my $info = shift;
	print STDERR "\n\033[1;31m$info\033[0m\n\n";
}

sub printGrn {
	my $info = shift;
	print STDERR "\033[1;32m$info\033[0m\n";
}

sub printYlw {
	my $info = shift;
	print STDERR "\n\033[1;33m$info\033[0m\n\n";
}


sub usage {
	print <<END_OF_USAGE;

########## Msuite2: Multi-mode DNA methylation data analysis suite ##########

Author : Kun Sun (sunkun\@szbl.ac.cn)
Version: $version

\033[1;34mUsage: msuite [options] -x index -1/-U Read1.fq [ -2 Read2.fq ] -o out.dir\033[0m

Compulsory parameters:

  -1/-U Read1.fq   Specify the path to the files containing read 1
                   If your data is Paired-end, specify read 2 files using '-2' option
                   Note that if -U is used, '-2' will be ignored no matter it's set or not

                   If you have multiple files, use ',' to separate them or use '*' syntax
                   (Note that single quotation marks are required for '*' syntax)

  -x index         Specify the genome index
                   Please refer to README file on how to build index for Msuite

  -o out.dir       Specify the output directory
                   Note that your specified directory will be created if it does not exist
                   otherwise the files under that directory could get over-written


Optional parameters:

  -2 Read2.fq      Specify the path to the file containing read 2
                   Use this parameter if your data is generated in paired-end mode

                   If you have multiple files, use ',' to separate them or use '*' syntax
                   Note that all files must be correctly paired in '-1' and '-2'

  -3               Use 3-letter alignment (default)
  -4               Use 4-letter alignment
                   Note that the above two options are mutually exclusive

  -m BS/TAPS       Specify the library protocol (default: BS)
                   Note that only 'TAPS' and 'BS' are acceptable

  -c cycle         Specify the seqeuencing cycles of the data (default: auto-detect)

  -k kit           Specify the library preparation kit (default: illumina)
                   Note that the current version supports 'illumina', 'nextera' and 'bgi'

  --phred33        Read cycle quality scores are in Phred33 format (default)
  --phred64        Read cycle quality scores are in Phred64 format
                   Note that the above two options are mutually exclusive

  -q score         The minimum quality score to keep the cycle (default: 20)
                   Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred

                   Sometimes quality scores start from 35 ('#') in the FASTQ files,
                   in this case you could adjust '-q' option, e.g., '--phred33 -q 22'

  --minsize size   Minimum read size to be kept for alignment (default: 36)

  --cut-r1-head N  Cut the head N cycles in read1 (default: 0)
  --cut-r1-tail N  Cut the tail N cycles in read1 (default: 0)
  --cut-r2-head N  Cut the head N cycles in read2 (default: 0)
  --cut-r2-tail N  Cut the tail N cycles in read2 (default: 0)

                   Note that the total cut basepairs in read 1 and 2 must be the same
                   (i.e., cut-r1-head + cut-r1-tail must equal to cut-r2-head + cut-r2-tail)

                   In BS-seq, read2 starts from the 3'-end which frequently suffer from
                   DNA damage issues, and the DNA repair step in library prepraration usually
                   uses un-methylated Cytosines which lead to bias in DNA methylation calling
                   as commonly seen in the M-bias plot.

  --minins MIN     Minimum insert size (default: 0)
  --maxins MAX     Maximum insert size (default: 1000)
                   Note that the above two options will be ignored for Single-End data

  --align-only     Stop after alignment (i.e., do not perform DNA methylation call and
                   visualization around TSS; default: not set)

  --CpH            Set this flag to call methylation status of CpH sites (default: not set)

  -p threads       Specify how many threads should be used (default: use all threads)

  -h/--help        Show this help information and quit
  -v/--version     Show the software version and quit

Please refer to README file for more information.

END_OF_USAGE
}

sub detect_cycle {
	my $file = shift;
	print "INFO: auto-detect read cycles using '$file'.\n";
	unless( -s "$file" ) {
		printRed( "Error: input file '$file' does not exist!" );
		exit(255);
	}
	if( $file =~ /\.gz$/ ) {
		open FQ, "zcat $file |" or die("$!");
	} else {
		open FQ, "$file" or die("$!");
	}
	my $MAX = 10000;	##read the first 1000 reads, and use the largest read size
	my $num = 0;
	my $cycle = 0;
	while( <FQ> ) {
		++ $num;
		last if $num >= $MAX;

		my $seq = <FQ>;
		<FQ>;
		<FQ>;
		chomp( $seq );
		$cycle = length($seq) if length($seq) > $cycle;
	}
	close FQ;

	print "INFO: Sequencing cycle is determined to be $cycle.\n";
	return $cycle;
}

sub check_dependency {
	my $bowtie2 = `which bowtie2 2>/dev/null`;
	chomp( $bowtie2 );
	my $bowtie2_ver;
	if( $bowtie2 ) {
		my $info = `$bowtie2 --version`;
		my @l = split /\n/, $info;
		foreach ( @l ) {
			if( /version (\S+)/ ) {
				$bowtie2_ver = $1;
				last;
			}
		}
		$bowtie2_ver = "unknown" unless $bowtie2_ver;
		print "INFO: bowtie2 (version $bowtie2_ver) found at '$bowtie2'.\n";
	} else {
		printRed( "Fatal error: could not locate 'bowtie2' in your path!" );
	exit 20;
	}

	my $samtools = `which samtools 2>/dev/null`;
	chomp( $samtools );
	if( $samtools ) {
		my $info = `$samtools 2>&1`;
		my $ver;
		my @l = split /\n/, $info;
		foreach ( @l ) {
			$ver = $1 if /^Version:\s(\S+)/;
		}
		$ver = "unknown" unless $ver;
		print "INFO: samtools (version $ver) found at '$samtools'.\n";
	} else {
		printRed( "Fatal error: could not locate 'samtools' in your path!" );
		exit 20;
	}

	my $R = `which R`;
	chomp( $R );
	if( $R ) {
		my $info = `$R --version`;
		my $ver;
		my @l = split /\n/, $info;
		foreach ( @l ) {
			$ver = $1 if /version (\S+)/;
		}
		$ver = "unknown" unless $ver;
		print "INFO: R (version $ver) found at '$R'.\n";

	} else {
		printRed( "Fatal error: could not locate 'R' in your path!" );
		exit 20;
	}

	return ($bowtie2, $bowtie2_ver, $samtools, $R);
}

sub check_index {
	my $Msuite2 = shift;
	my $indexID = shift;

	$Msuite2 = "$Msuite2/index";

	unless( -s "$Msuite2/$indexID/chr.info" && -d "$Msuite2/$indexID/indices" && -d "$Msuite2/$indexID/fasta" ) {
        printRed( "Fatal error: invalid index ($indexID)!\nPlease refer to README file for how to add indice to Msuite2." );
		exit 11;
	}
	unless( -s "$Msuite2/$indexID/tss.ext.bed" ) {
		printYlw( "WARNING: TSS annotation file for $indexID is missing/empty!" );
		return 0;
	}

	return 1;
}

sub makefile_perchr_v2 {
	my $MsuiteBin = shift;
	my $samtools  = shift;
	my $chrinfo   = shift;
	my $samheader = shift;
	my $seqMode   = shift;	## se or pe
	my $makefile  = shift;
	my $maxins    = shift || 1000;
	my $THREAD    = shift || 1;

	my $job = "";
	my $mkf = "";

	open IN, "$chrinfo" or die( "$!" );
	while( <IN> ) {
		chomp;
		my ($C, $size) = split /\t/;	## chr size
		my $chr = "chr$C";
		$chr = "Lambda" if $C eq 'L';

		$job .= " $chr.srt.bam";
		$mkf .= "$chr.srt.bam: chr$C.sam rhr$C.sam\n";
		$mkf .= "\t\@$MsuiteBin/rmdup.w.$seqMode $maxins chr$C.sam chr$C >chr$C.rmdup.log\n";
		$mkf .= "\t\@$MsuiteBin/rmdup.c.$seqMode $size $maxins rhr$C.sam rhr$C >rhr$C.rmdup.log\n";
		$mkf .= "\t\@cat $samheader chr$C.rmdup.sam rhr$C.c2w.sam | samtools view -bS - | samtools sort -o $chr.srt.bam -\n\n";
	}
	close IN;

	open MK, ">$makefile" or die( "$!" );
	print MK
"../Msuite2.final.bam.bai: ../Msuite2.final.bam
	$samtools index -\@ $THREAD ../Msuite2.final.bam
	\@cat *rmdup.log > ../Msuite2.rmdup.log
	$samtools index -\@ $THREAD Lambda.srt.bam
	mv Lambda.srt.bam* ../

../Msuite2.final.bam:$job
	$samtools merge -\@ $THREAD ../Msuite2.final.bam chr*.srt.bam

$mkf";
}

1;

