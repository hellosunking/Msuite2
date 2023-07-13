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
		cerr << "\nUsage: " << argv[0] << " <chr.info> <Msuite2.sam> <output.directory> [thread=1]\n"
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
	char **Reads = new char * [ READS_PER_BATCH ];
	for(unsigned int i=0; i!=READS_PER_BATCH; ++i ) {
		Reads[i] = new char[ MAX_SAMLINE_SIZE ];
	}

	unsigned int total = 0;
	// working loop
	while( true ) {
		// load files
		unsigned int loaded = 0;
		while( true ) {
			if( fgets( Reads[ loaded ], MAX_SAMLINE_SIZE, insam ) == NULL ) break;
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
//			cerr << "Thread " << tn << '\n';

			samRecord readsam;
			register char *psam;
			register bool endC;		// indicators: "C" at the end, "G" at the front, and WATSON strand
			register char QendC;	// quality score for the "C"
			char tail_added[8];		// store the modified cigar tail

			string curr_chr;
			///////////////////////////////////////////// T->C based on conversion log
			for( int index=start; index!=end; ++index ) {
//				cerr << "working " << index << "\n";
				psam = Reads[ index ];
				// split the sam record
				parseSAM( psam, readsam );	// NOTE HERER!!!

				// check whether it is a primary alignment, discard it if not
				int flag = 0;
				for(register int j=readsam.flag; psam[j]!='\t'; ++j) {
					flag *= 10;
					flag += psam[j] - '0';
				}
				if( flag & 256 ) {      // this is a secondary alignment
					continue;
				}

				// get chr
				curr_chr = psam[ readsam.chr ];
				register int i = readsam.chr + 1;
				while( psam[i] != '\t' ) {
					curr_chr += psam[i];
					++ i;
				}
//				cerr << " " << curr_chr << '\n';

				it = updatedSAM.find( curr_chr );
				if( it == updatedSAM.end() ) {
					cerr << "ERROR: no such chr for " << psam+readsam.seqName << "\n" ;
					continue;
				}
				FILE *fp = updatedSAM.find( curr_chr )->second;

				//// deal seqName
				// check endC marker
				i = 0;
				if( psam[1] == KEEP_QUAL_MARKER ) {	// there is a C at the end
					if( psam[0] == REPLACEMENT_CHAR_AT ) {
						QendC = '@';
					} else {
						QendC = psam[0];
					}
					endC = true;
					i = 2;
//					cerr << " has endC, qual=" << QendC << "\n";
				} else {	// no C at the end
					endC = false;
				}

				// process the conversion log
//				cerr << " T -> C\n";
				register char * rseq = psam + readsam.seq;
				if( psam[i] != CONVERSION_LOG_END ) {	// there are changes in this read
					register unsigned int j = 0;
					for( ; psam[i] != CONVERSION_LOG_END; ++i ) {
						if ( psam[i] == CONVERSION_LOG_SEPARATOR ) {	// one change met
							rseq[j] = 'C';	//convert back to 'C; read 1 is on WATSON
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
					rseq[j] = 'C';	//convert back to C
				}
				readsam.seqName = i + 1;	// position in seqName for the read ID
//				cerr << "    " << psam+read1sam.seqName << "\n";

				// deal with pos and CIGAR;
//				cerr << " CIGAR\n";
				if( endC ) {
					// add 1M to the end of CIGAR, 'C' to seq and QendC to qual
					// CIGAR: xM[yID]zM
//					psam[ read1sam.cigar   - 1 ] = '\0';
					psam[ readsam.qual      -1 ] = '\0';
					psam[ readsam.remaining -1 ] = '\0';

					//psam[ read1sam.mateflag - 1 ] = '\0';
					add_1M_to_cigar_end(psam + readsam.cigar, readsam.mateflag - readsam.cigar - 1, tail_added);

					fprintf(fp, "%s%s\t%sC\t%s%c\t%s",
							psam+readsam.seqName, tail_added, psam+readsam.mateflag,
							psam+readsam.qual, QendC, psam+readsam.remaining);
				} else {	// do not need to update CIGAR
					fprintf(fp, "%s", psam+readsam.seqName);
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
		delete [] Reads[i];
	}
	delete [] Reads;

	return 0;
}

