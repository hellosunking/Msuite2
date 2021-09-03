#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
#include <tr1/unordered_set>
#include <stdio.h>
//#include <stdlib.h>
//#include <memory.h>
#include "common.h"
#include "util.h"

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
 *
 * Additional tags in SAM by bowtie2
AS:i:<N> Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if SAM record is for an aligned read.
XS:i:<N> Alignment score for second-best alignment. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if the SAM record is for an aligned read and more than one alignment was found for the read.
YS:i:<N> Alignment score for opposite mate in the paired-end alignment. Only present if the SAM record is for a read that aligned as part of a paired-end alignment.
XN:i:<N> The number of ambiguous bases in the reference covering this alignment. Only present if SAM record is for an aligned read.
XM:i:<N> The number of mismatches in the alignment. Only present if SAM record is for an aligned read.
XO:i:<N> The number of gap opens, for both read and reference gaps, in the alignment. Only present if SAM record is for an aligned read.
XG:i:<N> The number of gap extensions, for both read and reference gaps, in the alignment. Only present if SAM record is for an aligned read.
NM:i:<N> The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and deletions) needed to transform the read string into the reference string. Only present if SAM record is for an aligned read.
YF:Z:<N> String indicating reason why the read was filtered out. See also: Filtering. Only appears for reads that were filtered out.
MD:Z:<S> A string representation of the mismatched reference bases in the alignment. See SAM format specification for details. Only present if SAM record is for an aligned read.
*/

int main( int argc, char *argv[] ) {
	if( argc != 5 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.size> <max.insertion> <in.c.sam> <out.prefix>\n\n"
			 << "This program is designed to remove the duplicate reads and revert crick to watson chain.\n"
			 << "Minimum score to keep the alignment: " << MIN_ALIGN_SCORE_KEEP << '\n'
			 << "Note that for duplicated reads, 1 random one will be kept.\n\n";
		return 1;
	}

	int chrsize = atoi( argv[1] );
	if( chrsize == 0 ) {
		cerr << "ERROR: incorrect chr size!\n";
		exit(1);
	}
	++ chrsize;	// to ease the reversion step

	int maxinsertion = atoi( argv[2] );
	if( maxinsertion == 0 ) {
		cerr << "ERROR: incorrect insertion size!\n";
		exit(1);
	}
	++ maxinsertion;

	// prepare file
	ifstream fin(argv[3]);
	if( fin.fail() ) {
		cerr << "Error: could not read file '" << argv[2] << "'!\n";
		exit( 10 );
	}

	string outfile = argv[4];
	outfile += ".rmdup.sam";	// this file is for meth-calling
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write sam file!\n";
		fin.close();
		exit( 11 );
	}

	outfile = argv[4];
	outfile += ".c2w.sam";	// this file is for the final alignment
	ofstream fc2w( outfile.c_str() );
	if( fc2w.fail() ) {
		cerr << "Error: could not write sam file!\n";
		fin.close();
		fout.close();
		exit( 12 );
	}

	// load sam file
	unordered_set<uint64_t> samHit;
	register uint64_t key;
	register unsigned int total = 0;
	register unsigned int discard = 0;
	register unsigned int dup = 0;
	int * size = new int [ maxinsertion ];

	string read1, read2;
	// r1 and r2 have the same chr and score
	// flag will be overwritten so do not need to remember it
	stringstream ss1, ss2;
	string name1, name2, flag, chr, score_str, cigar1, cigar2;
	string mateflag, matepos;
	string seq1, qual1, seq2, qual2;
	int pos1, pos2, fragSize;
	int score;
	string tmp;

	unsigned int rev_pos1, rev_pos2;
	vector<int> revhelper;
	string rev_cigar1, rev_cigar2;
	string rev_s1, rev_q1, rev_s2, rev_q2;
	string addTag1, addTag2;	// additional tags
	string :: const_reverse_iterator it;

//	cerr << "Loading sam file ...\n";
	while( true ) {
		getline( fin, read1 );
		if( fin.eof() )break;
		getline( fin, read2 );
		++ total;

//		if( ! (lineNum & 0x3fffff) ) {
//			cerr << '\r' << lineNum << " lines loaded.";
//		}

//131171  99 c1 145801355 42 100M = 113  258 TATCTCCTA HHHHHHHAAA AS:i:0 XN:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 MD:Z:100 YS:i:0 YT:Z:CP
//131171 147 c1 145801513 42 100M = 155 -258 AACCTAATT HHHHHHHHHH AS:i:0 XN:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 MD:Z:100 YS:i:0 YT:Z:CP

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

		if( pos1 > pos2 ) {	// problematic reads, discard
			++ discard;
			continue;
		}

		key = pos1;
		key <<= 32;
		key |= fragSize;

		if( samHit.find( key ) == samHit.end() ){	// key is not found, this is NOT a duplicate
			samHit.insert( key );
			fout << read1 << '\n' << read2 << '\n';
			++ size[ fragSize ];

			//// revert to real-watson chain
			// Read 1
			// reclaculate pos
			pos1 += get_readLen_from_cigar( cigar1 ) - 1;
			rev_pos1 = chrsize - pos1;
			// cigar
			revert_cigar( cigar1, rev_cigar1, revhelper );
			// sequence and quality: make reverse compliment
			ss1 >> seq1 >> qual1;
			rev_s1.clear();
			for(it=seq1.crbegin(); it!=seq1.crend(); ++it) {
				switch( *it ) {
					case 'A': rev_s1+='T'; break;
					case 'C': rev_s1+='G'; break;
					case 'G': rev_s1+='C'; break;
					case 'T': rev_s1+='A'; break;
					default : rev_s1+=*it;
				}
			}
			rev_q1.assign( qual1.crbegin(), qual1.crend() );

			// remaining tags by bowtie2: I will keep AS and NM tags
			addTag1 = "\tXG:Z:GA";	// to be compatible to Msuite1
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
			// if there is NO AS an NM tags, now remaining is NULL
			addTag1 += '\n';

			// Read 2
			ss2 >> score_str >> cigar2 >> mateflag >> matepos >> tmp >> seq2 >> qual2;

			pos2 += get_readLen_from_cigar( cigar2 ) - 1;
			rev_pos2 = chrsize - pos2;
			revert_cigar( cigar2, rev_cigar2, revhelper );

			rev_s2.clear();
			for(it=seq2.crbegin(); it!=seq2.crend(); ++it) {
				switch( *it ) {
					case 'A': rev_s2+='T'; break;
					case 'C': rev_s2+='G'; break;
					case 'G': rev_s2+='C'; break;
					case 'T': rev_s2+='A'; break;
					default : rev_s2+=*it;
				}
			}
			rev_q2.assign( qual2.crbegin(), qual2.crend() );
			// remaining tags by bowtie2: I will keep AS and NM tags
			addTag2 = "\tXG:Z:GA";	// to be compatible to Msuite1
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
			// if there is NO AS an NM tags, now remaining is NULL
			addTag2 += '\n';

			// write output
			chr[0] = 'c';	// I use rhr for reversed chromosomes, now change back to chr
			// read1 flag is ALWAYS 83; read2 flag is always 163; mateflag is always '='
			// I will output R2 first to speed-up sorting
			fc2w << name2 << "\t163\t" << chr << '\t' << rev_pos2 << '\t' << score_str << '\t' << rev_cigar2
				 << "\t=\t" << rev_pos1 << '\t'  << fragSize << '\t' << rev_s2 << '\t' << rev_q2 << addTag2
				 << name1 << "\t83\t" << chr << '\t' << rev_pos1 << '\t' << score_str << '\t' << rev_cigar1
				 << "\t=\t" << rev_pos2 << "\t-" << fragSize << '\t' << rev_s1 << '\t' << rev_q1 << addTag1;
		} else {	// this is a duplicate, discard it
			++ dup;
		}
	}
	fin.close();;
	fout.close();
	fc2w.close();

	cout << argv[3] << '\t' << total << '\t' << discard << '\t' << dup << '\n';

	outfile = argv[4];
	outfile += ".size";
	fout.open( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write SAM file!\n";
		exit( 1 );
	}
	for( int i=1; i!=maxinsertion; ++i ) {
		fout << i << '\t' << size[i] << '\n';
	}
	delete [] size;

	return 0;
}

