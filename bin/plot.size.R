#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# This program is part of Msuite2.
# Date: Ju1 2021
#

argv = commandArgs(T);
if( length(argv) != 3 ) {
	print( 'usage: R --slave --args <out.pdf> <w.size> <c.size> < plot.R' );
	q();
}

wsize = read.table( argv[2], head=T );
csize = read.table( argv[3], head=T );
if( sum(wsize[,2]) == 0 || sum(csize[,2]) == 0 ) {	## no reads, could be no lambda spike-in
	q();
}
wsize[,2] = wsize[,2] / sum(wsize[,2]) * 100;
csize[,2] = csize[,2] / sum(csize[,2]) * 100;
#Size	Count
#99	1741
#100	1740

if( nrow(wsize) == 0 || nrow(csize) == 0 ) {
	q();
}

## optimize plotting parameters
for( j in 1:nrow(wsize) ) {
	if( wsize[j,2] > 0.01 ) {
		break;
	}
}

if ( wsize[j,1] > 100 ) {
	xmin = 100;
} else if ( wsize[j,1] > 50 ) {
	xmin = 50;
} else {
	xmin = 0;
}

n = wsize[nrow(wsize),1];
for ( j in n:1 ) {
	if( wsize[j,2] > 0.01 ) {
		break;
	}
}
if ( wsize[j,1] > 300 ) {
	xmax = wsize[j,1];
} else {
	xmax = 300;
}

## plotting
outfileName = paste0(argv[1], ".pdf");
pdf( outfileName );
par( mar=c(5,5,1,1) );
plot( wsize[,2] ~ wsize[,1], type='l', lwd=2, col='red', xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
lines(csize[,2] ~ csize[,1], lwd=2, col='blue' );
legend( 'topright', c('Watson', 'Crick'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.5);
dev.off();

outfileName = paste0(argv[1], ".png");
png( outfileName, width=600, height=400 );
par( mar=c(5,5,1,1) );
plot( wsize[,2] ~ wsize[,1], type='l', lwd=2, col='red', xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
lines(csize[,2] ~ csize[,1], lwd=2, col='blue' );
legend( 'topright', c('Watson', 'Crick'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.5);
dev.off();

