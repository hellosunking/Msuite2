#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h>
#include <string.h>

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
*/

typedef struct {
	unsigned short C;
	unsigned short T;
}methCall;

typedef struct {
	unsigned int cov_C;
	unsigned int cov_T;
}methCoverage;

void usage( const char * prg ) {
	cerr << "\nUsage: " << prg << " <genome.info> <Msuite.meth.call> <tss.ext.bed> <TAPS|BS>\n"
		 << "\nThis program is a component of TAPSuite, designed to profile methylation signal around TSS.\n\n";
}

int main( int argc, char *argv[] ) {
	if( argc != 5 ) {
		usage( argv[0] );
		return 2;
	}

	bool TAPS;
	string mode = argv[4];
	if( mode=="TAPS" || mode=="taps" ) {
		TAPS = true;
//		cerr << "Use TAPS\n";
	} else if( mode=="BS" || mode=="bs" ) {
		TAPS = false;
//		cerr << "Use BS\n";
	} else {
		cerr << "Error: Unknown protocol! Must be TAPS or BS!\n";
		exit( 3 );
	}

	// load genome information
//	cerr << "Loading INFO file " << argv[1] << '\n';
	ifstream fin( argv[1] );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << argv[1] << " !\n";
		exit( 101 );
	}
	map<string, methCall*> methProfile;
	stringstream ss;
	string line, chr, ignore_str;
	register unsigned int chrsize;
	methCall* p_mc;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		if( line[0] == '#' )continue;
		ss.str( line );
		ss.clear();
		ss >> chr >> chrsize;
		p_mc = new methCall[chrsize];
		memset( p_mc, 0, chrsize*sizeof(methCall) );
		line = "chr";	// the chr.info does not has chr-prefix
		line += chr;
		methProfile.insert( pair<string, methCall*>(line, p_mc) );
	}
	fin.close();

	// load meth call
//	cerr << "Loading METH file " << argv[2] << '\n';
	fin.open( argv[2] );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << argv[2] << " to read!\n";
		exit( 102 );
	}
	map<string, methCall*> :: iterator mit;
	map<string, methCall*> :: iterator no_such_chr = methProfile.end();
	register int locus, unk, wC, wT, cC, cT;
	while( 1 ) {
		getline( fin, line );
		//chr   Locus   Total   wC      wT      wOther  Context cC      cT      cOther
		if( fin.eof() )break;

		if( line[0] == '#' )continue;
		ss.str( line );
		ss.clear();
		ss >> chr >> locus >> ignore_str >> wC >> wT >> ignore_str >> ignore_str >> cC >> cT;
		mit = methProfile.find( chr );
		if( mit == no_such_chr )continue;

		p_mc = mit->second;
		p_mc[locus].C = wC + cC;
		p_mc[locus].T = wT + cT;
	}
	fin.close();

	// load region file
//	cerr << "Loading BED file " << argv[3] << '\n';
	fin.open( argv[3] );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << argv[3] << " to read!\n";
		exit( 103 );
	}

	methCoverage mc;
	map<int, methCoverage> mTSS_w, mTSS_c;
	map<int, methCoverage> :: iterator mTSS_it, mTSS_it2;
	unsigned int start, end;
	char strand;
	while( true ) {
		getline( fin, line );
		//chr7      137613934       137614134       CREB3L2 3800    -
		if( fin.eof() ) break;

		ss.str( line );
		ss.clear();
		ss >> chr >> start >> end >> ignore_str >> locus >> strand;
		mit = methProfile.find( chr );
		if( mit == no_such_chr )continue;

		p_mc = mit->second;
		mc.cov_C = 0;
		mc.cov_T = 0;
		for( register int i=start+1; i<=end; ++i ) {
			mc.cov_C += p_mc[ i ].C;
			mc.cov_T += p_mc[ i ].T;
		}

		if( strand == '+' ) {
			mTSS_it = mTSS_w.find( locus );
			if( mTSS_it == mTSS_w.end() ) {
				mTSS_w.insert( pair<int, methCoverage>(locus, mc) );
				mTSS_it = mTSS_w.begin();
				mTSS_it->second.cov_C = mc.cov_C;
				mTSS_it->second.cov_T = mc.cov_T;
//				cerr << "Add " << locus << " for +\t" << mc.C << "\t" << mc.T << "\n";
			} else {
				mTSS_it->second.cov_C += mc.cov_C;
				mTSS_it->second.cov_T += mc.cov_T;
//				cerr << "Found " << locus << " for +\t" << mTSS_it->second.C << "\t" << mTSS_it->second.T << "\n";
			}
		} else {
			mTSS_it = mTSS_c.find( locus );
			if( mTSS_it == mTSS_c.end() ) {
				mTSS_c.insert( pair<int, methCoverage>(locus, mc) );
				mTSS_it = mTSS_c.begin();
				mTSS_it->second.cov_C = mc.cov_C;
				mTSS_it->second.cov_T = mc.cov_T;
//				cerr << "Add " << locus << " for - " << mc.C << "\t" << mc.T << "\n";
			} else {
//				cerr << "Found " << locus << " for -\t" << mTSS_it->second.C << "\t" << mTSS_it->second.T << "\n";
				mTSS_it->second.cov_C += mc.cov_C;
				mTSS_it->second.cov_T += mc.cov_T;
			}
		}
//		cerr << chr<<'\t'<<start<<'\t'<<end << '\t' << locus << '\t' << strand << '\t'
//				<< mc.cov_C << '\t' << mc.cov_T << '\t'
//				<< mTSS_it->second.cov_C << '\t' << mTSS_it->second.cov_T << '\n';
	}
	fin.close();

	// output
//	cerr << "Output\n";
	cout << "Distance\tWatson\tCrick\n";
	for( mTSS_it=mTSS_w.begin(); mTSS_it!=mTSS_w.end(); ++mTSS_it ) {
		cout << mTSS_it->first << '\t';
		if( TAPS ) {
			unsigned int total = mTSS_it->second.cov_C + mTSS_it->second.cov_T;
			if( total == 0 ) {
				cout << "NA\t";
			} else {
				cout << mTSS_it->second.cov_T * 100.0 / total << "\t";
			}

			mTSS_it2 = mTSS_c.find( mTSS_it->first );
			if( mTSS_it2 == mTSS_c.end() ) {
				cout << "NA\n";
			} else {
				total = mTSS_it2->second.cov_C + mTSS_it2->second.cov_T;
				if( total == 0 ) {
					cout << "NA\n";
				} else {
					cout << mTSS_it2->second.cov_T * 100.0 / total << "\n";
				}
			}
		} else {
			unsigned int total = mTSS_it->second.cov_C + mTSS_it->second.cov_T;
			if( total == 0 ) {
				cout << "NA\t";
			} else {
				cout << mTSS_it->second.cov_C * 100.0 / total << "\t";
			}

			mTSS_it2 = mTSS_c.find( mTSS_it->first );
			if( mTSS_it2 == mTSS_c.end() ) {
				cout << "NA\n";
			} else {
				total = mTSS_it2->second.cov_C + mTSS_it2->second.cov_T;
				if( total == 0 ) {
					cout << "NA\n";
				} else {
					cout << mTSS_it2->second.cov_C * 100.0 / total << "\n";
				}
			}
//			cerr << mTSS_it->first << '\t' << mTSS_it->second.C << '\t' << mTSS_it->second.T << '\t'
//										   << mTSS_it2->second->C << '\t' << mTSS_it2->second.T << '\n';
		}
	}

	// free memory
	for( mit=methProfile.begin(); mit!=no_such_chr; ++mit ) {
		delete [] mit->second;
	}
}


