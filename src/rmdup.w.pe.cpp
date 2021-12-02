#include <iostream>
#include <fstream>
//#include <string>
#include <sstream>
//#include <map>
//#include <set>
//#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <stdio.h>
//#include <stdlib.h>
//#include <memory.h>
#include "common.h"

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
 *
*/

int main( int argc, char *argv[] ) {
	if( argc != 4 ) {
		cerr << "\nUsage: " << argv[0] << " <max.insertion> <in.w.sam> <out.prefix>\n\n"
			 << "This program is designed to remove the duplicate reads that have the same start and end/strand.\n"
			 << "Minimum score to keep the alignment: " << MIN_ALIGN_SCORE_KEEP << '\n'
			 << "Note that for duplicated reads, a random keep one will be kept.\n\n";

		return 2;
	}

	int maxinsertion = atoi( argv[1] );
	if( maxinsertion == 0 ) {
		cerr << "Error: incorrect insertion size!\n";
		exit( 1 );
	}
	++ maxinsertion;

	// prepare file
	ifstream fin(argv[2]);
	if( fin.fail() ) {
		cerr << "Error: could not read file '" << argv[2] << "'!\n";
		exit( 1 );
	}

	string outfile = argv[3];
	outfile += ".rmdup.sam";
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write SAM file!\n";
		fin.close();
		exit( 1 );
	}

	// load sam file
	unordered_set<uint64_t> samHit;
	register uint64_t key;
	register unsigned int total = 0;
	register unsigned int discard = 0;
	register unsigned int dup = 0;

	stringstream ss1, ss2;
	string read1, read2;
	string name1, name2, flag, chr, score_str, cigar1, cigar2;
	string mateflag, matepos;
	string seq1, qual1, seq2, qual2, addTag1, addTag2;
	int pos1, pos2, fragSize;
	int score;
	string tmp;

	int * size = new int [ maxinsertion ];

//	cerr << "Loading sam file ...\n";
	while( true ) {
		getline( fin, read1 );
		if( fin.eof() )break;
		getline( fin, read2 );
		++ total;

//		if( ! (lineNum & 0x3fffff) ) {
//			cerr << '\r' << lineNum << " lines loaded.";
//		}

//131171	99	c1	145801355	42	100M	=	145801513	258	TATCTCCTAGGAAACTC	HHHHHHHAAA	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:0	YT:Z:CP
//131171	147	c1	145801513	42	100M	=	145801355	-258	AACCTAATTCATTCTGGGT	HHHHHHHHHHHHH	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:0	YT:Z:CP

		ss1.str( read1 );
		ss1.clear();
		ss1 >> name1 >> flag >> chr >> pos1 >> score >> cigar1 >> mateflag >> matepos >> fragSize;
		//note that the reads are always on FORWARD strand in Msuite2
		if( score < MIN_ALIGN_SCORE_KEEP || fragSize >= maxinsertion ) {
			++ discard;
			continue;
		}

		ss2.str( read2 );
		ss2.clear();
		ss2 >> name2 >> flag >> chr >> pos2;

		if( pos1 > pos2 ) {	// problematic alignment, discard
			++ discard;
			continue;
		}

		key = pos1;
		key <<= 32;
		key |= fragSize;

		if( samHit.find( key ) == samHit.end() ) {	// key is not found, this is NOT a duplicate
			samHit.insert( key );
			++ size[ fragSize ];

			// process bowtie2 tags
			// remaining tags by bowtie2: I will keep AS and NM tags
			ss1 >> seq1 >> qual1;
			addTag1.clear();
			bool AStag = false;
			bool NMtag = false;
			while( ss1.rdbuf()->in_avail() ) {
				ss1 >> tmp;
				if( tmp[0]=='A' && tmp[1]=='S' ) {
					addTag1 += '\t';
					addTag1 += tmp;
					AStag = true;
					if( NMtag )break;
				} else if( tmp[0]=='N' && tmp[1]=='M' ) {
					addTag1 += '\t';
					addTag1 += tmp;
					NMtag = true;
					if( AStag )break;
				}
			}
			addTag1 += "\tXG:Z:CT\n";	// to be compatible with Msuite1

			// Read 2
			ss2 >> score_str >> cigar2 >> mateflag >> matepos >> tmp >> seq2 >> qual2;
			addTag2.clear();
			AStag = false;
			NMtag = false;
			while( ss2.rdbuf()->in_avail() ) {
				ss2 >> tmp;
				if( tmp[0]=='A' && tmp[1]=='S' ) {
					addTag2 += '\t';
					addTag2 += tmp;
					AStag = true;
					if( NMtag )break;
				} else if( tmp[0]=='N' && tmp[1]=='M' ) {
					addTag2 += '\t';
					addTag2 += tmp;
					NMtag = true;
					if( AStag )break;
				}
			}
			addTag2 += "\tXG:Z:CT\n";

			// write output
			fout << name1 << "\t99\t" << chr << '\t' << pos1 << '\t' << score_str << '\t' << cigar1
				 << "\t=\t" << pos2 << '\t'  << fragSize << '\t' << seq1 << '\t' << qual1 << addTag1
				 << name2 << "\t147\t" << chr << '\t' << pos2 << '\t' << score_str << '\t' << cigar2
				 << "\t=\t" << pos1 << "\t-" << fragSize << '\t' << seq2 << '\t' << qual2 << addTag2;
		} else {	// this is a duplicate, discard it
			++ dup;
		}
	}
	fin.close();;
	fout.close();

	cout << argv[2] << '\t' << total << '\t' << discard << '\t' << dup << '\n';

	outfile = argv[3];
	outfile += ".size";
	fout.open( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write size file!\n";
		exit( 1 );
	}
	for( int i=1; i!=maxinsertion; ++i ) {
		fout << i << '\t' << size[i] << '\n';
	}
	delete [] size;
	
	return 0;
}

