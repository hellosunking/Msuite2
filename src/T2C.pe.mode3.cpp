#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <string>
#include <omp.h>
#include <vector>
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
 *       @xx;xx#    => for reads without endC
 *       @#         => for reads without endC and conversions
 *   read 2:
 *       @xx;xx#    => for reads without frontG
 *       @#         => for reads without frontG and conversions
 *
 * Note that in mode3, there is NO endC and frontG issue!!!
 * so this program is speed up the analysis for mode 3 as it is the most common working mode.
 *
**/

int main( int argc, char *argv[] ) {
	if( argc < 4 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <Msuite2.PE.sam> <output.directory> [thread=1]\n"
			 << "\nThis program is part of Msuite2, designed to change T back to C in the alignment file (mode 3 ONLY)."
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

	// working loop
	register unsigned int cnt = 0;
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
		cnt += loaded;
		cerr << "\rProgress: " << cnt << " lines loaded.";

		omp_set_num_threads( thread );
		#pragma omp parallel
		{
			unsigned int tn    = omp_get_thread_num();
			unsigned int start = loaded * tn / thread;
			unsigned int end   = loaded * (tn+1) / thread;

			samRecord read1sam, read2sam;
			register char *psam;
			string curr_chr;

			///////////////////////////////////////////// T->C based on conversion log
			for( unsigned int index=start; index!=end; ++index ) {
//				cerr << index << "\n";
				// deal read 1: it always has a smaller genomic coordinate in Msuite2
				psam = R1[index];
				// split the sam record
				parseSAM( psam, read1sam );

				// check whether it is a primary alignment, discard it if not;
				int flag = 0;
				for(register int j=read1sam.flag; psam[j]!='\t'; ++j) {
					flag *= 10;
					flag += psam[j] - '0';
				}
				if( flag & 256 ) {      // this is a secondary alignment
					continue;
				}

				// get chr
				curr_chr = psam[ read1sam.chr ];
				unsigned int i = read1sam.chr + 1;
				while( psam[i] != '\t' ) {
					curr_chr += psam[i];
					++ i;
				}

				// process the conversion log
				i = 0;
				register char * r1seq = psam + read1sam.seq;
				if( psam[i] != CONVERSION_LOG_END ) {	// there are changes in this read
					register unsigned int j = 0;
					for( ; psam[i] != CONVERSION_LOG_END; ++i ) {
						if ( psam[i] == CONVERSION_LOG_SEPARATOR ) {	// one change met
							r1seq[j] = 'C';	//convert back to 'C; read 1 is on WATSON
							j = 0;
						} else {
							j <<= 4;
							if( psam[i] <= '9' ) {	//0-9, seqName > '0' is always true
								j += psam[i] - '0';
							} else {	//a-f
								j += psam[i] - 87;	// 'a'-10
							}
						}
					}
					r1seq[j] = 'C';	//convert back to C
				}
				char *R1offset = psam + i + 1;	// position in seqName for the read ID
//				cerr << "    " << psam+read1sam.seqName << "\n";

				/////////////////////////////////////////////////////////////////////////////////////////////////////
				// deal read 2: should be on the reverse chain
				psam = R2[index];
				// split the sam record
				parseSAM( psam, read2sam );

				i = 0;
				register unsigned int len = read2sam.qual - read2sam.seq - 1 - 1;
				// length of SEQ - 1; R2 is rev-comp in SAM

				// process the conversion log
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
				char *R2offset = psam + i + 1;

				// write updated sam
				fprintf( updatedSAM.find(curr_chr)->second, "%s%s", R1offset, R2offset);
			}// end for loop
		}// end multi-thread loop
		if( feof(insam) )break;
	} // end file loop

	cerr << "\rDone. Totally " << cnt << " lines loaded.\n";
	fclose( insam );

	unordered_map<string, FILE *> :: iterator it;
	for( it=updatedSAM.begin(); it!=updatedSAM.end(); ++it ) {
		fclose( it->second );
	}

	for(unsigned int i=0; i!=READS_PER_BATCH; ++i) {
		delete [] R1[i];
		delete [] R2[i];
	}
	delete [] R1;
	delete [] R2;

	return 0;
}

