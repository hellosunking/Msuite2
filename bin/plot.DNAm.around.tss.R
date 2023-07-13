#
# Author: Kun Sun (hellosunking@foxmail.com)
# This program is part of Msuite.
# Date: Dec 2019
#

args = commandArgs( T );
if( length(args) != 2 ) {
	print( "usage: <in.DNAm.around.tss.stat> <out.prefix>" );
	q();
}

infile    = args[1];
outprefix = args[2];

dat = read.table( infile, head=T );
#Distance	Watson	Crick
#-4000	79.5963141728828	78.6795626576955
#-3800	79.1227330240405	79.9830004249894

## PDF
out = paste( outprefix, 'pdf', sep='.' );
pdf( out, width=9, height=6 );
par( mar=c(5,5,1,1) );
plot(  dat$Watson ~ dat$Distance, type='l', lwd=4, col="red", ylim=c(0,100),
		xlab="Distance from TSS on forward chain (bp)", ylab="Methylation density (%)", cex.lab=1.5 );
lines( dat$Crick  ~ dat$Distance, lwd=4, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.5 );
dev.off();

## PNG
out = paste( outprefix, 'png', sep='.' );
png( out, width=1080, height=540 );
par( mar=c(5,5,1,1) );
plot( dat$Watson ~ dat$Distance, type='l', lwd=4, col="red", ylim=c(0, 100),
		xlab="Distance from TSS (bp)", ylab="Methylation density (%)", cex.lab=1.5 );
lines( dat$Crick  ~ dat$Distance, lwd=4, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.5 );
dev.off();

