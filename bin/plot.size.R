#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# This program is part of Msuite2.
# Date: Ju1 2021
#

args = commandArgs(T);
if( length(args) < 3 ) {
	print( 'usage: R --slave --args <out.pdf> <w.size> <c.size> [cut.size=0] < plot.R' );
	q();
}

cutsize = 0;
if( length(args) > 3 ) {
	cutsize = as.numeric( args[4] );
}

wsize = read.table( args[2], head=T );
csize = read.table( args[3], head=T );
if( sum(wsize$Count) == 0 || sum(csize$Count) == 0 ) {	## no reads, could be no lambda spike-in
	q();
}

wsize$Size = wsize$Size + cutsize;
csize$Size = csize$Size + cutsize;

wsize$Count = wsize$Count / sum(wsize$Count) * 100;
csize$Count = csize$Count / sum(csize$Count) * 100;
#Size	Count
#99	1741
#100	1740

if( nrow(wsize) == 0 || nrow(csize) == 0 ) {
	q();
}

## optimize plotting parameters
for( j in 1:nrow(wsize) ) {
	if( (! is.na(wsize[j,2])) && wsize[j,2] > 0.01 ) {
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
	if( (! is.na(wsize[j,2])) && wsize[j,2] > 0.01 ) {
		break;
	}
}
if ( wsize[j,1] > 300 ) {
	xmax = wsize[j,1];
} else {
	xmax = 300;
}

## plotting
outfileName = paste0(args[1], ".pdf");
pdf( outfileName );
par( mar=c(5,5,1,1) );
plot( wsize$Count ~ wsize$Size, type='l', lwd=2, col='red', xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
lines(csize$Count ~ csize$Size, lwd=2, col='blue' );
legend( 'topright', c('Watson', 'Crick'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.5);
dev.off();

outfileName = paste0(args[1], ".png");
png( outfileName, width=600, height=400 );
par( mar=c(5,5,1,1) );
plot( wsize$Count ~ wsize$Size, type='l', lwd=2, col='red', xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
lines(csize$Count ~ csize$Size, lwd=2, col='blue' );
legend( 'topright', c('Watson', 'Crick'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.5);
dev.off();

