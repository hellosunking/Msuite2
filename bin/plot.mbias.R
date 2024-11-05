#
# Author: Kun Sun (hellosunking@foxmail.com)
# This program is part of Msuite2.
# Date: Jul 2021
#

args = commandArgs( T );
if( length(args) < 4 ) {
	print( "usage: <in.w.mbias> <in.c.mbias> <mode=BS|TAPS> <out.prefix> [maxSize=inf]" );
	q();
}

maxSize = 1e9;
if( length(args) > 4 ) {
	maxSize = as.numeric( args[5] );
}

wbias = read.table( args[1], head=T, comment="%" );
cbias = read.table( args[2], head=T, comment="%" );
#Cycle   C       T
#1       2303    2341
#2       2321    2389

wbias= subset(wbias, Cycle <= maxSize)
cbias= subset(cbias, Cycle <= maxSize)

## calculate DNAm
if( args[3] == "BS" || args[3] == "bs" ) {
	w = wbias$C / ( wbias$C + wbias$T )* 100;
	c = cbias$C / ( cbias$C + cbias$T )* 100;
} else {
	w = wbias$T / ( wbias$C + wbias$T )* 100;
	c = cbias$T / ( cbias$C + cbias$T )* 100;
}

out = paste0( args[4], '.pdf' );
pdf( out, width=6, height=4 );
par( mar=c(5,5,1,1) );
plot(  w ~ wbias$Cycle, type='o', pch=1, col="red", ylim=c(0, 100),
		xlab="Sequencing cycles (bp)", ylab="Methylation level (%)", cex.lab=1.5 );
lines( c ~ cbias$Cycle, type='o', pch=8, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c("red", "blue"), pch=c(1, 8), horiz=T, bty='n', cex=1.5 );
dev.off();

out = paste0( args[4], '.png' );
png( out, width=600, height=400 );
par( mar=c(5,5,1,1) );
plot(  w ~ wbias$Cycle, type='o', pch=1, col="red", ylim=c(0, 100),
		xlab="Sequencing cycles (bp)", ylab="Methylation level (%)", cex.lab=1.5 );
lines( c ~ cbias$Cycle, type='o', pch=8, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c("red", "blue"), pch=c(1, 8), horiz=T, bty='n', cex=1.5 );
dev.off();

