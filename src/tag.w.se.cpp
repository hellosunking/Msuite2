#include <iostream>
#include <fstream>
//#include <string>
#include <sstream>
//#include <map>
//#include <set>
//#include <tr1/unordered_map>
//#include <tr1/unordered_set>
#include <stdio.h>
#include <stdlib.h>
//#include <memory.h>
#include "common.h"

using namespace std;
//using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Aug 2022
 *
*/

int main( int argc, char *argv[] ) {
	if( argc != 4 ) {
		cerr << "\nUsage: " << argv[0] << " <max.insertion=placeholder> <in.w.sam> <out.prefix>\n\n"
			 << "This program is designed to fix the tags in SAM (without rmdup).\n"
			 << "Minimum score to keep the alignment: " << MIN_ALIGN_SCORE_KEEP << ".\n\n";

		return 2;
	}

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
	register unsigned int total = 0;
	register unsigned int discard = 0;
	register unsigned int dup = 0;

	stringstream ss;
	string read;
	string name, flag, chr, cigar;
	string mateflag, matepos, matedist;
	string seq, qual, addTag;
	int pos, score;
	string tmp;

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

		// process bowtie2 tags
		// remaining tags by bowtie2: I will keep AS and NM tags
		ss >> cigar >> mateflag >> matepos >> matedist >> seq >> qual;
		addTag.clear();
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
		addTag += "\tXG:Z:CT\n";

		// write output
		fout << name << "\t0\t" << chr << '\t' << pos << '\t' << score << '\t' << cigar
			 << "\t*\t0\t0\t" << seq << '\t' << qual << addTag;
	}
	fin.close();;
	fout.close();

	cout << argv[2] << '\t' << total << '\t' << discard << '\t' << dup << '\n';
	return 0;
}

