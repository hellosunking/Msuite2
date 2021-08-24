#!/usr/bin/perl
#
# Author: Kun Sun (sunkun@szbl.ac.cn)
# This program is part of Msuite2
# Date: Jul 2021
#

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use MsuiteUtil qw($version $ver);

#print STDERR "Usage: $0 [data.dir=.]\n";

my $dir = $ARGV[0] || '.';
my @color = ( '#E8DAEF', '#D6EAF8' );

open OUT, ">$dir/Msuite2.report/index.html" or die("$!");
select OUT;

## load configuration file
## this file should be prepared by Msuite program while NOT the user
my $absPath=`readlink -m $0`;
chomp( $absPath );
my $Msuite = dirname($absPath);
$Msuite =~ s/bin\/?$/Msuite/;

open CONF, "$dir/Msuite2.conf" or die( "$!" );
my %conf;
#my $parameter = "<tr bgcolor=\"$color[0]\"><td>Msuite path</td><td>$Msuite</td></tr>\n";
my $parameter = "";
my $i = 1;
while( <CONF> ){
	chomp;
	next if /^#/;
	next unless /\S/;
	my @l = split /\t/;
	$conf{ $l[0] } = $l[1];
	$l[1] =~ s/:/<br \/>/;	## process R1:R2 files
	$parameter .= "<tr bgcolor=\"$color[$i]\"><td>$l[0]</td><td>$l[1]</td></tr>\n";
	$i = 1 - $i;
}
close CONF;

my $pe   = ( $conf{"Sequencing mode"}  =~ /^P/i   ) ? 1 : 0;
my $TAPS = ( $conf{"Library protocol"} =~ /TAPS/i ) ? 1 : 0;
my $alignonly = ( $conf{"Align-only mode"} =~ /on/i ) ? 1 : 0;
my $CpH       = ( $conf{"Call CpH"} =~ /yes/i ) ? 1 : 0;

print <<HTMLHEADER;
<html>
<head>
<title>Msuite2 Analysis Report</title>
<style type="text/css">
td {
text-align: left;
padding-left: 10px;
}
#fqstat td {
text-align: center;
font-weight:bold;
}
</style>
</head>
<body>
<h1>Msuite2 Analysis Report</h1>

HTMLHEADER

############################################
print "<h2>Alignment statistics</h2>\n";

## load trim log
open LOG, "$dir/Msuite2.trim.log" or die( "$!" );
my $line = <LOG>;
$line =~ /(\d+)/;
my $total = $1;
$line = <LOG>;	##Dropped : 0
$line =~ /(\d+)/;
my $dropped = $1;
my $trim = $total - $dropped;
close LOG;

## load aligner log
open LOG, "$dir/Msuite2.rmdup.log" or die( "$!" );
my ($waligned, $wdiscard, $wduplicate) = (0, 0, 0);
my ($caligned, $cdiscard, $cduplicate) = (0, 0, 0);
my $cntL = 0;
while( <LOG> ) {
	my @l = split /\t/;	## chr total discard dup
	if( $l[0] =~ /^chr/ ) {
		$waligned   += $l[1];
		$wdiscard   += $l[2];
		$wduplicate += $l[3];
	} else {	# rhrXXX
		$caligned   += $l[1];
		$cdiscard   += $l[2];
		$cduplicate += $l[3];
	}

	$cntL += $l[1] if $l[0]=~/^[cr]hrL/;
}
close LOG;
my $aligned = $waligned + $caligned;
my $discard = $wdiscard + $cdiscard;
my $duplicate = $wduplicate + $cduplicate;
my $reported  = $aligned - $discard - $duplicate;

print "<table id=\"alignStat\" width=\"75%\">\n",
		"<tr bgcolor=\"$color[0]\"><td width=\"70%\"><b>Total input reads</b></td>",
		"<td width=\"30%\"><b>", digitalize($total), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td><b>After preprocessing</b></td>",
		"<td><b>", digitalize($trim), sprintf(" (%.2f %%)", $trim/$total*100), "</b></td></tr>\n";

print "<tr bgcolor=\"$color[0]\"><td><b>Total aligned reads</b></td>",
		"<td><b>", digitalize($aligned), sprintf(" (%.2f %%)", $aligned/$trim*100), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td>&nbsp;&nbsp;Forward chain</td>",
		"<td>&nbsp;&nbsp;", digitalize($waligned), sprintf(" (%.2f %%)", $waligned/$trim*100), "</td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td>&nbsp;&nbsp;Reverse chain</td>",
		"<td>&nbsp;&nbsp;", digitalize($caligned), sprintf(" (%.2f %%)", $caligned/$trim*100), "</td></tr>\n";

print "<tr bgcolor=\"$color[1]\"><td><b>Low-quality alignments</b></td>",
		"<td><b>", digitalize($discard), sprintf(" (%.2f %%)", $discard/$aligned*100), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td><b>PCR duplicates</b></td>",
		"<td><b>", digitalize($duplicate), sprintf(" (%.2f %%)", $duplicate/$aligned*100), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td><b>Reported alignments</b></td>",
		"<td><b>", digitalize($reported), sprintf(" (%.2f %%)", $reported/$aligned*100), "</b></td></tr>\n",
	 "</table>\n\n";

############################################
unless( $alignonly ) {
print "<h2>Methylation statistics</h2>\n";
## load CpG.meth log
open LOG, "$dir/Msuite2.CpG.meth.log" or die( "$!" );
my ($wC, $wT, $cC, $cT) = ( 0, 0, 0, 0 );
my $conversionL = 'NA';
while( <LOG> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	#chr Total.wC Total.wT Total.cC Total.cT
	if( $l[0] ne 'chrL' ) {
		$wC += $l[1];
		$wT += $l[2];
		$cC += $l[3];
		$cT += $l[4];
	} else {	## reads mapped to the lambda genome
		if( $cntL ) {
			$conversionL = sprintf( "%.2f %%", ($l[2]+$l[4])/($l[1]+$l[2]+$l[3]+$l[4])*100 );
		}
	}
}
my ($wm, $cm, $tm);
if( $TAPS ) {
	$wm = sprintf("%.2f", $wT/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cT/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wT+$cT)/($wC+$wT+$cC+$cT)*100);
} else {
	$wm = sprintf("%.2f", $wC/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cC/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wC+$cC)/($wC+$wT+$cC+$cT)*100);
}

print "<table id=\"methStat\" width=\"75%\">\n",
		"<tr bgcolor=\"$color[0]\"><td width=\"70%\"><b>Overall CpG methylation density</b></td>" ,
			"<td width=\"30%\"><b>$tm %</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td>&nbsp;&nbsp;Forward chain</td><td>&nbsp;&nbsp;$wm %</td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td>&nbsp;&nbsp;Reverse chain</td><td>&nbsp;&nbsp;$cm %</td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td><b>Reads mapped to Lambda genome</b></td>",
			"<td><b>", digitalize($cntL), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td><b>C-&gt;T conversion rate for lambda genome</b></td>",
			"<td><b>$conversionL</b></td></tr>\n";

if( $CpH ) {
## load CpH.meth log
open LOG, "$dir/Msuite2.CpH.meth.log" or die( "$!" );
($wC, $wT, $cC, $cT) = ( 0, 0, 0, 0 );
while( <LOG> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	#chr Total.wC Total.wT Total.cC Total.cT
	if( $l[0] ne 'chrL' ) {
		$wC += $l[1];
		$wT += $l[2];
		$cC += $l[3];
		$cT += $l[4];
	}
}
my ($wm, $cm, $tm);
if( $TAPS ) {
	$wm = sprintf("%.2f", $wT/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cT/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wT+$cT)/($wC+$wT+$cC+$cT)*100);
} else {
	$wm = sprintf("%.2f", $wC/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cC/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wC+$cC)/($wC+$wT+$cC+$cT)*100);
}
print "<tr bgcolor=\"$color[1]\"><td width=\"70%\"><b>Overall CpH methylation density</b></td>" ,
			"<td width=\"30%\"><b>$tm %</b></td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td>&nbsp;&nbsp;Forward chain</td><td>&nbsp;&nbsp;$wm %</td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td>&nbsp;&nbsp;Reverse chain</td><td>&nbsp;&nbsp;$cm %</td></tr>\n";
}

print "</table>\n\n";
}
###################################################
print '<h2>Base composition in the sequenced reads</h2>
<table id="fqstat">
	<tr><td>Read 1 raw sequence</td><td>Read 1 trimmed</td></tr>
	<tr>
		<td><img src="R1.fqstat.png" alt="Base composition in read 1 raw"></td>
		<td><img src="R1.trimmed.fqstat.png" alt="Base composition in read 1 trimmed"></td>
	</tr>
';
if( $pe ) {
	print 
'	<tr><td>Read 2 raw sequence</td><td>Read 2 trimmed</td></tr>
	<tr>
		<td><img src="R2.fqstat.png" alt="Base composition in read 2 raw"></td>
		<td><img src="R2.trimmed.fqstat.png" alt="Base composition in read 2 trimmed"></td>
	</tr>
';
}
print "</table>\n\n";

if( $pe ) {
	print 
'<h2>Fragment size distribution</h2>
	<img src="Msuite2.size.png" alt="fragment size distribution"><br />
';
}

unless( $alignonly ) {
	print
'<h2>Methylation level per chromosome</h2>
	<img src="DNAm.per.chr.png" alt="DNAm.per.chr"><br />

<h2>Methylation level around TSS</h2>
	<img src="DNAm.around.TSS.png" alt="Methylation level around TSS"><br />

<h2>M-bias plot</h2>
<ul>
	<li>Read 1</li><br />
	<img src="R1.mbias.png" alt="M-bias in read 1"><br />
';
	if( $pe ) {
		print
'	<li>Read 2</li><br />
	<img src="R2.mbias.png" alt="M-bias in read 2"><br />
';
	}
	print "</ul>\n\n";	
}

print "<h2>Analysis Parameters</h2>\n",
		"<table id=\"para\" width=\"80%\">\n",
		"<tr bgcolor=\"#888888\"><td width=\"30%\"><b>Option</b></td><td width=\"70%\"><b>Value</b></td></tr>\n",
		"$parameter</table>\n";

## HTML tail
my ($sec, $min, $hour, $day, $mon, $year, $weekday, $yeardate, $savinglightday) = localtime();
$year+= 1900;
my @month = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
$hour = "0$hour" if $hour < 10;
$min  = "0$min"  if $min  < 10;
$sec  = "0$sec"  if $sec  < 10;
$mon = $month[$mon];
my $time = "$hour:$min:$sec, $mon-$day-$year";

my $url = 'https://github.com/hellosunking/Msuite/';
print "<HR align=\"left\" width=\"80%\"/>\n",
		"<h4>Generated by <a href=\"$url\" target=\"_blank\">Msuite2</a> (version $ver) on $time.</h4>\n",
		"</body>\n</html>\n";

close OUT;

###############################################################
sub digitalize {
	my $v = shift;

	while($v =~ s/(\d)(\d{3})((:?,\d\d\d)*)$/$1,$2$3/){};
	return $v;
}

