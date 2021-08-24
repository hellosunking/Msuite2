#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <memory.h>
#include <map>
#include "common.h"
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
 * In this version, M-bias data is provided
*/

int main( int argc, char *argv[] ) {
	if( argc != 6 ) {
        cerr<< "\nUsage: " << argv[0] << " <chr.label> <chr.fa> <mode=BS|TAPS> <w.call> <c.call>\n"
			<< "\nThis program is a component of Msuite2, designed to merge CpH methylation from watson and crick strands.\n";
		return 2;
	}

	string seq;
	loadchr( argv[2], seq );
	int chrsize = seq.length() - 2 + 1;
	//seq contains 2 additional letters "X", "Y"

	bool mode;	// true is BS, false is TAPS
	if( strcmp(argv[3], "BS")==0 || strcmp(argv[3], "bs")==0 ) {
		mode = true;
	} else if( strcmp(argv[3], "TAPS")==0 || strcmp(argv[3], "taps")==0 ) {
		mode = false;
	} else {
		cerr << "ERROR: Unknown mode! Must be BS or TAPS.\n";
		exit(1);
	}

	map<int, pairedmeth> meth;
	map<int, pairedmeth> :: iterator it;

	stringstream ss;
	string line;
	register unsigned int pos, C, T, Z;
	pairedmeth pm;
	pm.wC=0; pm.wT=0; pm.wZ=0;
	pm.cC=0; pm.cT=0; pm.cZ=0;

	// watson
	ifstream fcall( argv[4] );
	if( fcall.fail() ) {
		cerr << "Error file: cannot open " << argv[4] << " to read!\n";
		exit(10);
	}
	while( true ) {
		getline( fcall, line );
		if( fcall.eof() ) break;

		if( line[0] == '#' ) continue;

		ss.clear();
		ss.str( line );
		ss >> pos >> C >> T >> Z;
		pm.wC = C;
		pm.wT = T;
		pm.wZ = Z;
		meth.insert( pair<int, pairedmeth>(pos, pm) );
	}
	fcall.close();

	// crick, which is NOT paired to watson
	fcall.open( argv[5] );
	if( fcall.fail() ) {
		cerr << "Error file: cannot open " << argv[5] << " to read!\n";
		exit(11);
	}
	pm.wC=0; pm.wT=0; pm.wZ=0;
	while( true ) {
		getline( fcall, line );
		if( fcall.eof() ) break;

		if( line[0] == '#' ) continue;

		ss.clear();
		ss.str( line );
		ss >> pos >> C >> T >> Z;
		pos = chrsize - pos;

		pm.cC = C;
		pm.cT = T;
		pm.cZ = Z;
		meth.insert( pair<int, pairedmeth>(pos, pm) );
	}
	fcall.close();

	//write output
	string output = argv[1];
	output += ".CpH.meth";
	ofstream fout( output.c_str() );
	if( fout.fail() ) {
		cerr << "Error file: cannot open " << output << " to write!\n";
		exit(12);
	}

	output = argv[1];
	output += ".CpH.meth.bedgraph";
	ofstream fbed( output.c_str() );
	if( fbed.fail() ) {
		cerr << "Error file: cannot open " << output << " to write!\n";
		exit(13);
	}

	//#chr	Locus	Total	wC	wT	wOther	Context	cC	cT	cOther
	int wC_total = 0;
	int wT_total = 0;
	int cC_total = 0;
	int cT_total = 0;
	for( it=meth.begin(); it!=meth.end(); ++it ) {
		pos = it->first;
		int total_valid = it->second.wC + it->second.wT + it->second.cC + it->second.cT;
		int total = total_valid + it->second.wZ + it->second.cZ;
		fout << argv[1] << '\t' << pos << '\t' << total << '\t'
			 << it->second.wC << '\t' << it->second.wT << '\t' << it->second.wZ << '\t'
			 << seq[pos-1] << seq[pos] << seq[pos+1] << '\t'
			 << it->second.cC << '\t' << it->second.cT << '\t' << it->second.cZ << '\n';

		wC_total += it->second.wC;
		wT_total += it->second.wT;
		cC_total += it->second.cC;
		cT_total += it->second.cT;

		float meth;
		if( mode ) {
			meth = (it->second.wC+it->second.cC)*100.0/total_valid;
		} else {
			meth = (it->second.wT+it->second.cT)*100.0/total_valid;
		}
		fbed << argv[1] << '\t' << pos-1 << '\t' << pos << '\t' << meth << '\n';
	}
	fout.close();
	fbed.close();

	cout << argv[1] << '\t' << wC_total << '\t' << wT_total << '\t' << cC_total << '\t' << cT_total << '\n';
}

