#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

using namespace std;

/*
 * Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
 * Date  : May 6, 2020
 *
**/

typedef struct {
	unsigned short C;
	unsigned short T;
}meth;

const unsigned int MAX_ERROR_COUNT = 5;
const unsigned int MAX_ERROR_PROPORTION = 10;

int main( int argc, char *argv[] ) {
	if( argc < 4 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <Msuite2.CpG.meth.call> <query.bed>\n"
			 << "\nThis program is designed to calculate the CpG coverage in the given regions.\n"
			 << "\n3 columns will be added to the input BED file: CpG.covered C.count T.count"
			 << "\nBy default, the result will be written to STDOUT, you may redirect it to a file.\n"
			 << "\nIn addition, CpG sites with too many non-CT calls (>10%) will be discarded (it could contain SNPs).\n\n";
		return 1;
	}

	// load info
//	cerr << "Loading info file " << argv[1] << " ...\n";
	ifstream fin( argv[1] );
	if( fin.fail() ) {
		cerr << "Error open info file.\n";
		return 1;
	}

	map<string, meth*> methCall;
	string line, chr, chr_tmp;
	stringstream ss;
	register unsigned int len;
	meth * p;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		ss.str( line );
		ss.clear();
		ss >> chr_tmp >> len;
		if( chr_tmp[0]=='c' && chr_tmp[1]=='h' && chr_tmp[2]=='r' ) {
			chr = chr_tmp;
		} else {
			chr = "chr";
			chr += chr_tmp;
		}
//		cerr << "  Adding " << chr << " ...\n";
		p = new meth[ len+1 ];
		memset( p, 0, (len+1)*sizeof(meth) );
		methCall.insert( pair<string, meth*>( chr, p ) );
	}
	fin.close();

	// load call
//	cerr << "Loading call file " << argv[2] << " ...\n";
	fin.open( argv[2] );
	if( fin.fail() ) {
		cerr << "Error open call file.\n";
		return 2;
	}
	string context;
	unsigned int pos, total_cover, wc, wt, wo, cc, ct, co, errors;
	map<string, meth*> :: iterator it, end;
	end = methCall.end();
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		ss.str( line );
		ss.clear();
		ss >> chr >> pos >> total_cover >> wc >> wt >> wo >> context >> cc >> ct >> co;

		errors = wo + co;
		if( errors >= MAX_ERROR_COUNT && errors*100/total_cover >= MAX_ERROR_PROPORTION )
			continue;

		it = methCall.find(chr);
		if( it == end )
			continue;		// no chr record
		p = it->second + pos;
		p->C = wc + cc;
		p->T = wt + ct;
	}
	fin.close();

	//load qery bed
	register unsigned int j, s, e;
	register unsigned int c, t, m;
//	cerr << "Loading query file " << argv[3] << "\n";
	fin.open( argv[3] );
	if( fin.fail() ) {
		cerr << "Error open query file: skip.\n";
		fin.close();
	}

	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		ss.str( line );
		ss.clear();
		ss >> chr >> s >> e;
		it = methCall.find( chr );
		if( it == end )
			continue;

		p = it->second;
		c = 0;
		t = 0;
		m = 0;
		++ s;	// fix the 0-base thing
		++ e;	// to facilitate the loop
		for( j=s; j!=e; ++j ) {
			c += p[j].C;
			t += p[j].T;
			if( p[j].C || p[j].T )
				++ m;
		}
		cout << line << '\t' << m << '\t' << c << '\t'<< t << '\n';
	}
	fin.close();

	for( it=methCall.begin(); it!=end; ++it ) {
		delete [] it->second;
	}

	return 0;
}
