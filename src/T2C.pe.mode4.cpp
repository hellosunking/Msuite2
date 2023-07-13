#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>
#include <unistd.h>
#include <unordered_map>
#include "common.h"
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
 * This program is designed to return the C>T conversions and split the data into chromosomes
 *
 * Input : Msuite2.sam (by bowtie2)
 * Output: OUTDIR/Msuite2.w1.sam OUTDIR/Msuite2.c2.sam, etc.
 *
 * In Ktrim, the conversion log format is placed BEFORE the raw seqName. Here is the rule:
 *   read 1:	## note that LINENUMBER has been removed in Msuite2's conversion step!!!
 *       @S|xx;xx$  => for reads with endC and C>T changes, 'S' is its quality score
 *       @xx;xx$    => for reads without endC
 *       @$         => for reads without endC and conversions
 *   read 2:
 *       @S|xx;xx$  => for reads with frontG and G>A changes, 'S' is its quality score
 *       @xx;xx$    => for reads without frontG
 *       @$         => for reads without frontG and conversions
 *
 * This could fasten the de-convert procedure since the whole seqName is untouched
**/

int main( int argc, char *argv[] ) {
	if( argc < 4 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <Msuite2.PE.sam> <output.directory> [thread=1]\n"
			 << "\nThis program is part of Msuite2, designed to change T back to C in the alignment file."
			 << "\nMulti-thread is supported, 4-8 threads are recommanded.\n\n";
		//cerr << "Rescue mode is ON.\n\n";
		return 2;
	}

	int thread = 1;
	if( argc > 4 ) {
		thread = atoi( argv[4] );
		if( thread <= 0 ) {
			cerr << "INFO: all threads will be used.\n";
			thread = omp_get_max_threads();
		}
	}
	// prepare files
	ifstream finfo( argv[1] );
	if( finfo.fail() ) {
		cerr << "Error: could not open chr.info file '" << argv[1] << "'!\n";
		exit( 1 );
	}

	stringstream ss;
	string line, chr, mchr;
	char outfile[128];
	unordered_map<string, FILE *> updatedSAM;
	unordered_map<string, FILE *> :: iterator it;
	while( true ) {
		getline( finfo, line );
		if( finfo.eof() )break;

		if( line[0] == '#' )continue;
		ss.str( line );
		ss.clear();
		ss >> chr;

		mchr = "chr";
		mchr += chr;
		sprintf( outfile, "%s/%s.sam", argv[3], mchr.c_str() );
		FILE * fp = fopen( outfile, "w" );
		if( fp == NULL ) {
			cerr << "FATAL: Could not open file " << outfile << " for write.\n";
			exit(3);
		}
		updatedSAM.emplace( pair<string, FILE *>(mchr, fp) );

		mchr = "rhr";
		mchr += chr;
		sprintf( outfile, "%s/%s.sam", argv[3], mchr.c_str() );
		fp = fopen( outfile, "w" );
		if( fp == NULL ) {
			cerr << "FATAL: Could not open file " << outfile << " for write.\n";
			exit(3);
		}
		updatedSAM.emplace( pair<string, FILE *>(mchr, fp) );
	}
	finfo.close();

	// input sam file
	FILE *insam = fopen( argv[2], "r" );
	if( insam == NULL ) {
		cerr << "Error: cannot open file " << argv[1] << "!\n";
		exit(11);
	}

	// main job
	char **R1, **R2;
	R1 = new char * [ READS_PER_BATCH ];
	R2 = new char * [ READS_PER_BATCH ];
	for(unsigned int i=0; i!=READS_PER_BATCH; ++i ) {
		R1[i] = new char[ MAX_SAMLINE_SIZE ];
		R2[i] = new char[ MAX_SAMLINE_SIZE ];
	}

	unsigned int total = 0;
	// working loop
	while( true ) {
		// load files
		unsigned int loaded = 0;
		while( true ) {
			if( fgets( R1[ loaded ], MAX_SAMLINE_SIZE, insam ) == NULL ) break;
			fgets( R2[ loaded ], MAX_SAMLINE_SIZE, insam );
			// note that fgets KEEPs the tail '\n'

			++ loaded;
			if( loaded == READS_PER_BATCH )break;
		}
		if( loaded == 0 ) break;
		total += loaded;
//		cerr << "\rProgress: " << total << " lines loaded.";

		omp_set_num_threads( thread );
		#pragma omp parallel
		{
			unsigned int tn    = omp_get_thread_num();
			unsigned int start = loaded * tn / thread;
			unsigned int end   = loaded * (tn+1) / thread;

			samRecord read1sam, read2sam;
			char newCGAR1[ MAX_CIGAR_SIZE ], newCGAR2[ MAX_CIGAR_SIZE ];
			register char *r1cigar, *r2cigar;
			char tail_cigar1[8];
			char tail_cigar2[8];
//			char tmpR1[MAX_SAMLINE_SIZE], tmpR2[MAX_SAMLINE_SIZE];

			register char *psam;
			register bool endC, frontG;		// indicators: "C" at the end, "G" at the front, and WATSON strand
			register char QendC, QfrontG;	// quality score for the "C"

			string curr_chr;
			///////////////////////////////////////////// T->C based on conversion log
			for( int index=start; index!=end; ++index ) {
//				cerr << "working " << index << "\n";
				// deal read 1: it always has a smaller genomic coordinate in Msuite2
				psam = R1[index];
				// split the sam record
//				cerr << " Split R1\n";
				parseSAM( psam, read1sam );	// NOTE HERER!!!

				// check whether it is a primary alignment, discard it if not;
				int flag = 0;
				for(register int j=read1sam.flag; psam[j]!='\t'; ++j) {
					flag *= 10;
					flag += psam[j] - '0';
				}
				if( flag & 256 ) {      // this is a secondary alignment
					continue;
				}

				//// deal seqName
				// check endC marker
				register int i = 0;
				if( psam[1] == KEEP_QUAL_MARKER ) {	// there is a C at the end
					if( psam[0] == REPLACEMENT_CHAR_AT ) {
						QendC = '@';
					} else {
						QendC = psam[0];
					}
					endC  = true;
					i = 2;
//					cerr << "    has endC, qual=" << QendC << "\n";
				} else {	// no C at the end
					endC = false;
				}

				// process the conversion log
//				cerr << " T -> C\n";
				register char * r1seq = psam + read1sam.seq;
				if( psam[i] != CONVERSION_LOG_END ) {	// there are changes in this read
					register unsigned int j = 0;
					for( ; psam[i] != CONVERSION_LOG_END; ++i ) {
						if ( psam[i] == CONVERSION_LOG_SEPARATOR ) {	// one change met
							r1seq[j] = 'C';	//convert back to 'C; read 1 is on WATSON
							j = 0;
						} else {
							j <<= 4;
							if( psam[i]<='9' ) {	//0-9, seqName > '0' is always true
								j += psam[i] - '0';
							} else {	//a-f
								j += psam[i] - 87;	// 'a'-10
							}
						}
					}
					r1seq[j] = 'C';	//convert back to C
				}
				read1sam.seqName = i + 1;	// position in seqName for the read ID
//				cerr << "    " << psam+read1sam.seqName << "\n";

				// deal with pos and CIGAR;
//				cerr << " CIGAR\n";
				r1cigar = psam + read1sam.cigar;
				if( endC ) {
					// add 1M to the end of CIGAR
					// CIGAR: xM[yID]zM[tS]
////				psam[ read1sam.cigar    - 1 ] = '\0';
					//psam[ read1sam.mateflag - 1 ] = '\0';
					unsigned int len = read1sam.mateflag - read1sam.cigar - 1;
					add_1M_to_cigar_end(r1cigar, read1sam.mateflag - read1sam.cigar - 1, tail_cigar1);
//					cerr << "endC: " << r1cigar << ", " << tail_cigar1 << "\n";
//				} else {	// do not need to update CIGAR
				}

				/////////////////////////////////////////////////////////////////////////////////////////////////////
				// deal read 2: should be on the reverse chain
//				cerr << " R2\n";
				psam = R2[index];
				// split the sam record
				parseSAM( psam, read2sam );	// NOTE HERER!!!

				// look for the conversion log start
				// look for frontG marker
				register unsigned int len;	// R1 is rev-comp in SAM
				if( psam[1] == KEEP_QUAL_MARKER ) {
					if( psam[0] == REPLACEMENT_CHAR_AT ) {
						QfrontG = '@';
					} else {
						QfrontG = psam[0];
					}
					frontG  = true;
					i = 2;
					len = read2sam.qual - read2sam.seq - 1;	// length of SEQ
//					cerr << "    has frontG, qual=" << QfrontG << "\n";
				} else {
					frontG = false;
					i = 0;
					len = read2sam.qual - read2sam.seq - 1 - 1;	// length of SEQ - 1
				}

				// process the conversion log
//				cerr << " T -> C\n";
				register char * r2seq = psam + read2sam.seq;
				if( psam[i] != CONVERSION_LOG_END ) {	// there are G>A conversions
					register unsigned int j = 0;
					for( ; psam[i] != CONVERSION_LOG_END; ++i ) {
						if ( psam[i] == CONVERSION_LOG_SEPARATOR ) {	// one change meet
							r2seq[len-j] = 'C';	// convert back to C; Read 2 is on CRICK strand
							j = 0;
						} else {
							j <<= 4;
							if( psam[i] <= '9' ) {	// 0-9
								j += psam[i] - '0';
							} else {	// a-f
								j += psam[i] - 87;	// 'a'-10
							}
						}
					}
					r2seq[len-j] = 'C';	//convert back to G
				}
				read2sam.seqName = i + 1;

				// deal with pos and CIGAR
//				cerr << " CIGAR\n";
				r2cigar = psam + read2sam.cigar;
				if( frontG ) {
					// Read2 has been reverse-complimented in SAM record
					// add a 'C' to the _end_ of sequence, add 1M at the _end_ of CIGAR
					// CIGAR: xM[yID]zM
////				psam[ read2sam.cigar    - 1 ] = '\0';
					psam[ read2sam.matedist - 1 ] = '\0';

					add_1M_to_cigar_end(r2cigar, read2sam.mateflag - read2sam.cigar - 1, tail_cigar2);
//					cerr << "frontG: " << r2cigar << ", " << tail_cigar2 << "\n";
//				} else {	// do not need to update CIGAR
				}
//				cerr << "  CIGAR: " << r1cigar << ", " << r2cigar << "\n";

				// write output according to chrosomes, after each batch
//				cerr << " Output\n";
				curr_chr = psam[ read2sam.chr ];
				i = read2sam.chr + 1;
				while( psam[i] != '\t' ) {
					curr_chr += psam[i];
					++ i;
				}

				it = updatedSAM.find( curr_chr );
				if( it == updatedSAM.end() ) {
					cerr << "ERROR: no such chr for " << psam+read1sam.seqName << "\n" ;
					continue;
				}
				FILE *fp = updatedSAM.find( curr_chr )->second;
				// Because this program will be run in multi-thread mode, R1 and R2 MUST be written in 1 fprintf call
//				cerr << "  Update sam, chr=" << chr;
//				ss.str( "" );
//				ss.clear();

				if( frontG ) {
					//read2: CIGAR, fragSize, seq, qual has changed
					register char *p2 = R2[index];
					p2[ read2sam.qual     - 1 ] = '\0';
					p2[ read2sam.remaining- 1 ] = '\0';

					// update fragment size
					psam = R1[index];
					register int fragsize = atoi( psam + read1sam.matedist );
					++ fragsize;

					// read 1
					if( endC ) {	// with endC, CIGAR & fragSize changed
//						cerr << ", endC+frontG\n";
						psam[ read1sam.matedist - 1 ] = '\0';
						psam[ read1sam.qual     - 1 ] = '\0';
						psam[ read1sam.remaining- 1 ] = '\0';

						fprintf(fp, "%s%s\t%s\t%d\t%sC\t%s%c\t%s%s%s\t%s\t-%d\t%sC\t%s%c\t%s",
								psam+read1sam.seqName, tail_cigar1, psam+read1sam.mateflag, fragsize,
								psam+read1sam.seq, psam+read1sam.qual, QendC, psam+read1sam.remaining,
								p2+read2sam.seqName, tail_cigar2, p2+read2sam.mateflag, fragsize,
								p2+read2sam.seq, p2+read2sam.qual, QfrontG, p2+read2sam.remaining);
					} else {	// no endC, only fragSize changed
//						cerr << ", frontG only\n";
						// NO endC, then CIGAR1, seq/qual/remaining is not affected
						psam[ read1sam.matedist - 1 ] = '\0';
						fprintf(fp, "%s\t%d\t%s%s%s\t%s\t-%d\t%sC\t%s%c\t%s",
								psam+read1sam.seqName, fragsize, psam+read1sam.seq,
								p2+read2sam.seqName, tail_cigar2, p2+read2sam.mateflag, fragsize,
								p2+read2sam.seq, p2+read2sam.qual, QfrontG, p2+read2sam.remaining);
					}
				} else {	// no frontG, then fragSize is NOT affected
					// read1
					psam = R1[index];
					if( endC ) { // endC, seq, qual changed
//						cerr << ", endC\n";
						psam[ read1sam.qual     - 1 ] = '\0';
						psam[ read1sam.remaining- 1 ] = '\0';

						fprintf(fp, "%s%s\t%sC\t%s%c\t%s%s",
								psam+read1sam.seqName, tail_cigar1, psam+read1sam.mateflag,
								psam+read1sam.qual, QendC, psam+read1sam.remaining,
								R2[index]+read2sam.seqName);
					} else {	// no endC and frontG, ALMOST all elements are not changed
//						cerr << ", null\n";
						fprintf(fp, "%s%s", psam+read1sam.seqName, R2[index]+read2sam.seqName);
					}
				}
			}// end for loop
		} // end multi-thread loop

		if( feof(insam) )break;
	}	// end file loop
	cerr << "\rDone. Totally " << total << " lines loaded.\n";
	fclose( insam );

	for( it=updatedSAM.begin(); it!=updatedSAM.end(); ++it ) {
		fclose( it->second );
	}

	for(unsigned int i=0; i!=READS_PER_BATCH; ++i ) {
		delete [] R1[i];
		delete [] R2[i];
	}
	delete [] R1;
	delete [] R2;

	return 0;
}

