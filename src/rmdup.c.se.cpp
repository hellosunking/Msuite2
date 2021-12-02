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
*/

int main( int argc, char *argv[] ) {
	if( argc != 5 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.size> <max.insertion=placeholder> <in.c.sam> <out.prefix>\n\n"
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

	// prepare file
	ifstream fin(argv[3]);
	if( fin.fail() ) {
		cerr << "Error: could not read file '" << argv[2] << "'!\n";
		exit( 10 );
	}

	string outfile = argv[4];
	outfile += ".rmdup.sam";
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write sam file!\n";
		fin.close();
		exit( 11 );
	}

	outfile = argv[4];
	outfile += ".c2w.sam";
	ofstream fc2w( outfile.c_str() );
	if( fc2w.fail() ) {
		cerr << "Error: could not write sam file!\n";
		fin.close();
		fout.close();
		exit( 12 );
	}

	// load sam file
	unordered_set<int> samHit;
	register int key;
	register unsigned int total = 0;
	register unsigned int discard = 0;
	register unsigned int dup = 0;

	string read;
	// r1 and r2 have the same chr and score
	// flag will be overwritten so do not need to remember it
	stringstream ss;
	string name, flag, chr, cigar;
	string mateflag, matepos, matedist;
	string seq, qual;
	int pos, score;
	string tmp;

	unsigned int rev_pos;
	vector<int> revhelper;
	string rev_cigar;
	string rev_s, rev_q;
	string addTag;	// additional tags
	string :: const_reverse_iterator it;

//	cerr << "Loading sam file ...\n";
	while( true ) {
		getline( fin, read );
		if( fin.eof() )break;
		++ total;

//		if( ! (lineNum & 0x3fffff) ) {
//			cerr << '\r' << lineNum << " lines loaded.";
//		}

//131171	99	c1	145801355	42	100M	=	145801513	258	TATCTCCTAGGAAACTC	HHHHHHHAAA	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:0	YT:Z:CP
//131171	147	c1	145801513	42	100M	=	145801355	-258	AACCTAATTCATTCTGGGT	HHHHHHHHHHHHH	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:0	YT:Z:CP

		ss.str( read );
		ss.clear();
		ss >> name >> flag >> chr >> pos >> score;
		//note that the reads are always on FORWARD strand in Msuite2
		if( score < MIN_ALIGN_SCORE_KEEP ) {
			++ discard;
			continue;
		}

		key = pos;
		if( samHit.find( key ) == samHit.end() ){	// key is not found, this is NOT a duplicate
			samHit.insert( key );
			fout << read << '\n';

			//// revert to real-watson chain
			ss >> cigar >> mateflag >> matepos >> matedist >> seq >> qual;
			// reclaculate pos
			pos += get_readLen_from_cigar( cigar ) - 1;
			rev_pos = chrsize - pos;
			// cigar
			revert_cigar( cigar, rev_cigar, revhelper );
			// sequence and quality: make reverse compliment
			rev_s.clear();
			for(it=seq.crbegin(); it!=seq.crend(); ++it) {
				switch( *it ) {
					case 'A': rev_s+='T'; break;
					case 'C': rev_s+='G'; break;
					case 'G': rev_s+='C'; break;
					case 'T': rev_s+='A'; break;
					default : rev_s+=*it;
				}
			}
			rev_q.assign( qual.crbegin(), qual.crend() );

			// remaining tags by bowtie2: I will keep AS and NM tags
			addTag.clear();	// to be compatible to Msuite1
			bool AStag = false;
			bool NMtag = false;
			while( ss.rdbuf()->in_avail() ) {
				ss >> tmp;
				if( tmp[0]=='A' && tmp[1]=='S' ) {
					addTag += '\t';
					addTag += tmp;
					AStag = true;
					if( NMtag )break;
				} else if( tmp[0]=='N' && tmp[1]=='M' ) {
					addTag += '\t';
					addTag += tmp;
					NMtag = true;
					if( AStag )break;
				}
			}
			// if there is NO AS an NM tags, now remaining is NULL
			addTag += "\tXG:Z:GA\n";

			// write output
			chr[0] = 'c';	// I use rhr for reversed chromosomes, now change back to chr
			// read1 flag is ALWAYS 83; read2 flag is always 163; mateflag is always '='
			// I will output R2 first to speed-up sorting
			fc2w << name << "\t16\t" << chr << '\t' << rev_pos << '\t' << score << '\t' << rev_cigar
				 << "\t*\t0\t0\t" << rev_s << '\t' << rev_q << addTag;
		} else {	// this is a duplicate, discard it
			++ dup;
		}
	}
	fin.close();;
	fout.close();
	fc2w.close();

	cout << argv[3] << '\t' << total << '\t' << discard << '\t' << dup << '\n';
	return 0;
}

