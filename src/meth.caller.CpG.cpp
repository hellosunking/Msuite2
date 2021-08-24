#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <memory.h>
#include <tr1/unordered_map>
#include "common.h"
#include "util.h"

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
 * In this version, M-bias data is provided
*/

// function declarations, the implementation is at the end of this file
// some functions are implemented in util.cpp
void deal_SE_CpG( const char *gfile, const char *samfile, const int cycle, const char *output );
void deal_PE_CpG( const char *gfile, const char *samfile, const int cycle, const char *output );

void callmeth_CpG_mbias( string &realSEQ, string &realQUAL, int pos,
				string &g, unordered_map<int, meth> &methcall, meth *mb, int cycle, bool rev );
void callmeth_mbias(     string &realSEQ, string &realQUAL, int pos,
				string &g, meth *mb, int cycle, bool rev );
void write_methcall( unordered_map<int, meth> &m, const char *pre, const char *suf );
void write_mbias( meth *m, int cycle, const char *pre, const char *suf );

int main( int argc, char *argv[] ) {
	if( argc != 6 ) {
        cerr<< "\nUsage: " << argv[0] << " <mode=SE|PE> <chr.fa> <chr.sam> <cycle> <output.prefix>\n"
			<< "\nThis program is a component of Msuite2, designed to call CpG methylation and M-bias from SAM file.\n"
			<< "Both SE/PE data are supported; indels are also supported.\n\n";
		return 2;
	}

	int cycle = atoi( argv[4] );
	if( cycle == 0 ) {
		cerr << "Error: Invalid cycle!\n";
		exit( 4 );
	}

	string mode = argv[1];
	if( mode=="SE" || mode=="se" ) {	// for SE data, only need to calculate target 1
		deal_SE_CpG( argv[2], argv[3], cycle, argv[5] );
	} else if( mode=="PE" || mode=="pe" ) {
		deal_PE_CpG( argv[2], argv[3], cycle, argv[5] );
	} else {
		cerr << "Error: Unknown mode! Must be PE or SE!\n";
		exit( 5 );
	}

	return 0;
}

// process SE data
void deal_SE_CpG( const char *gfile, const char *samfile, const int cycle, const char *output) {
	// load genome
	string g;
	loadchr( gfile, g );

	unordered_map<int, meth> methcall;
//	meth *methcall = new meth [ g.size() ];
//	memset( methcall, 0, sizeof(meth) * g.size() );
	meth *mbias = new meth[ cycle ];
	memset( mbias, 0, sizeof(meth) * cycle );

	// open sam file
	ifstream fsam( samfile );
	if( fsam.fail() ) {
		cerr << "Error file: cannot open " << samfile << " to read!\n";
		exit(200);
	}
//	cout << "Loading alignment " << samfile << " in SE mode ...\n";
//	unsigned int count = 0;
	string line, seqName, chr, cigar, seq, qual;
	int flag;
	string mateinfo, matepos, dist;   //fields that are ignored; all the sequence are converted to WATSON chain
	register unsigned int pos, score;
	stringstream ss;
	string realSEQ, realQUAL;   //these are CIGAR-processed seq and qual
	line.resize( MAX_SAMLINE_SIZE );
//	bool strand;	// strand is always TRUE in Msuite2
	// load sam file
	while( true ) {
		getline( fsam, line );
		if( fsam.eof() ) break;
		//14_R1	83	chr9	73301642	42	36M	=	73301399	-279	TCCTTCTCTCCCTC	GHHHHHHHHHH	XG:Z:GA

		ss.clear();
		ss.str( line );
		ss >> seqName >> flag >> chr >> pos >> score >> cigar >> mateinfo >> matepos >> dist >> seq >> qual;

		if( score < MIN_ALIGN_SCORE_METH ) {
			//cerr << "Discard " << seqName << " due to poor alignment score.\n";
			continue;
		}

		// process the CIGAR, handle the indels
		if( ! fix_cigar(cigar, realSEQ, realQUAL, seq, qual) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar << ") at line " << line << "!\n";
			continue;
		}

		// call CpG methylation
		callmeth_CpG_mbias( realSEQ, realQUAL, pos, g, methcall, mbias, cycle, false );

		// report progress for every 4 million reads
//		++ count;
//		if( ! (count & 0x003fffff) )
//			cout << '\r' << count << " lines loaded.";
	}
	fsam.close();
//	cout << '\r' << "Done: " << count << " lines loaded.\n";

	write_methcall( methcall, output, ".CpG.call" );
	write_mbias( mbias, cycle, output, ".R1.mbias" );

	delete [] mbias;
}

/////////////////////////////////////////////////////////////////////////////////////////
void deal_PE_CpG( const char *gfile, const char *samfile, const int cycle, const char *output) {
	// load genome
	string g;
	loadchr( gfile, g );

	unordered_map<int, meth> methcall;

	meth *mb1 = new meth[ cycle ];
	memset( mb1, 0, sizeof(meth) * cycle );
	meth *mb2 = new meth[ cycle ];
	memset( mb2, 0, sizeof(meth) * cycle );
	meth *mb3 = new meth[ cycle ];	// for overlapping reads; not used in the current version
	memset( mb3, 0, sizeof(meth) * cycle );

	// open sam file
	ifstream fsam( samfile );
	if( fsam.fail() ) {
		cerr << "Error file: cannot open " << samfile << " to read!\n";
		exit(200);
	}
//	cerr << "Loading alignment " << samfile << " in PE mode ...\n";

//	unsigned int count = 0;
	string line1, line2, seqName, chr, cigar1, seq1, qual1, cigar2, seq2, qual2, score2;
	register int flag;
	string mateinfo, matepos, dist;   //fields that are ignored; all the sequence are converted to WATSON chain
	register unsigned int pos1, pos2, score;
	stringstream ss;
	string realSEQ1, realQUAL1, realSEQ2, realQUAL2;   //these are CIGAR-processed seq and qual
	string mSEQ, mQUAL; //merged sequence and quality if read1 and read2 has overlap
	line1.resize( MAX_SAMLINE_SIZE );
	line2.resize( MAX_SAMLINE_SIZE );
	mSEQ.resize( MAX_MERGED_SEQ );
	mQUAL.resize( MAX_MERGED_SEQ );
//	bool strand;

	// load sam file
	while( true ) {
		getline( fsam, line1 );
		if( fsam.eof() ) break;
		getline( fsam, line2 );

		//14_R1	83	chr9	73301642	42	36M	=	73301399	-279	TCCTCCTTCTCTCCCTC	HHHHHHHHH	XG:Z:CT
		//14_R2	163	chr9	73301399	42	36M	=	73301642	279	TTTATTTTGATCCTGTA	DDCBA@?>=<;986420.

		ss.clear();
		ss.str( line1 );
		ss >> seqName >> flag >> chr >> pos1 >> score >> cigar1 >> mateinfo >> matepos >> dist >> seq1 >> qual1;
//		cerr << seqName << '\n';

		if( score < MIN_ALIGN_SCORE_METH ) {
			//cerr << "Discard " << seqName << " due to poor alignment score.\n";
			continue;
		}

		ss.clear();
		ss.str( line2 );
		ss >> seqName >> flag >> chr >> pos2 >> score2 >> cigar2 >> mateinfo >> matepos >> dist >> seq2 >> qual2;

		if( pos1 > pos2 ) {	// rare scenario that read2 contains read1!!! Mapping error?
//			cerr << "ERROR: Read2 contains Read1 in " << seqName << ", skip!\n";
			continue;
		}

		// process CIGAR 1, handle the indels
		realSEQ1.clear();
		realQUAL1.clear();
		if( ! fix_cigar( cigar1, realSEQ1, realQUAL1, seq1, qual1 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar1 << ") in " << seqName << "!\n";
			continue;
		}
		// process CIGAR 2, handle the indels
		realSEQ2.clear();
		realQUAL2.clear();
		if( ! fix_cigar( cigar2, realSEQ2, realQUAL2, seq2, qual2 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar2 << ") in " << seqName << "!\n";
			continue;
		}

		if( pos1 + realSEQ1.size() <= pos2 ) { //there is NO overlap
			callmeth_CpG_mbias( realSEQ1, realQUAL1, pos1, g, methcall, mb1, cycle, false );
			callmeth_CpG_mbias( realSEQ2, realQUAL2, pos2, g, methcall, mb2, cycle, true  );
		} else {	// there is overlap in read 1 and read 2
			//cerr << "Found overlap in " << seqName << '\n';
			if( pos2+realSEQ2.size() >= pos1+realSEQ1.size() ) {	// most case
				mSEQ.clear();
				mQUAL.clear();
				int len = pos2+realSEQ2.size() - pos1;
				//cerr << "len\t" << len << "\tp1=" << p1 << "\tp2=" << p2 << '\n';
				register int rs = realSEQ1.size();
				register int offset = pos2 - pos1;
				register unsigned int k;
				for( k=0; k != offset; ++k ) {	// read1 only
					mSEQ  += realSEQ1[k];		// can also use substr
					mQUAL += realQUAL1[k];
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				unsigned int j = 0;
				for( ; k != rs; ++k ) {	// overlapped region, peak the one with higher quality
					// j= k - offset;
					if( realQUAL1[k] >= realQUAL2[j] ) {
						mSEQ  += realSEQ1[k];
						mQUAL += realQUAL1[k];
					} else {
						mSEQ  += realSEQ2[j];
						mQUAL += realQUAL2[j];
					}
					++ j;
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				for( ; k != len; ++k ) {	//read2 only
					mSEQ  += realSEQ2[j];
					mQUAL += realQUAL2[j];
					++ j;
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";
				callmeth_CpG_mbias( mSEQ, mQUAL, pos1, g, methcall, mb3, cycle, false );
				callmeth_mbias( realSEQ1, realQUAL1, pos1, g, mb1, cycle, false );
				callmeth_mbias( realSEQ2, realQUAL2, pos2, g, mb2, cycle, true  );
			} else {	// rare case that R1 completely contains R2 => use R1 directly
				callmeth_CpG_mbias( realSEQ1, realQUAL1, pos1, g, methcall, mb1, cycle, false );
				callmeth_mbias( realSEQ2, realQUAL2, pos2, g, mb2, cycle, true  );
			}
		}
//		++ count;
//		if( ! (count & 0x003fffff) )
//			cout << '\r' << count << " lines loaded.";
	}
//	cerr << '\r' << "Done: " << count << " lines loaded.\n";
	fsam.close();

	// write meth call and  M-bias
//	cerr << "Output\n";
	write_mbias( mb1, cycle, output, ".R1.mbias" );
	write_mbias( mb2, cycle, output, ".R2.mbias" );
//	cerr << "Call\n";
	write_methcall( methcall, output, ".CpG.call" );

//	cerr << "Done.\n";
	delete [] mb1;
	delete [] mb2;
	delete [] mb3;
}

// call meth from sequence
void callmeth_CpG_mbias( string &seq, string &qual, int pos,
			string &g, unordered_map<int, meth> &methcall, meth *mb, int cycle, bool rev ) {
	unsigned int rs = seq.size();
	unsigned int os = rs - 1;	// offset for rev-cmp-ed R2
	char c1, c2;
	unordered_map<int, meth> :: iterator it;
	for( unsigned int i=0, j=pos; i!=rs; ++i, ++j) {
		if( qual[i] < MIN_BASEQUAL_SCORE )
			continue;

		c1 = g[j];
		c2 = g[j+1];

		if( (c1!='C' && c1!='c') || (c2!='G' && c2!='g') )	// not a CpG site
			continue;

		// m-bias
		int k = (rev) ? (os-i) : i;
		if( k < cycle ) {
			if( seq[i] == 'C' ) {
				++ mb[k].C;
			} else if( seq[i] == 'T' ) {
				++ mb[k].T;
			} else {
				++ mb[k].Z;
			}
		}

//		cerr << "Meet CpG on " << j << '\n';
		it = methcall.find( j );
		if( it == methcall.end() ) {
			meth m;
			if( seq[i] == 'C' ) {
				m.C=1; m.T=0; m.Z=0;
			} else if( seq[i] == 'T' ) {
				m.C=0; m.T=1; m.Z=0;
			} else {
				m.C=0; m.T=1; m.Z=1;
			}
			methcall.insert( pair<int, meth>( j, m ) );
		} else {
			if( seq[i] == 'C' ) {
				++ it->second.C;
			} else if( seq[i] == 'T' ) {
				++ it->second.T;
			} else {
				++ it->second.Z;
			}
		}
	}
}

// call mbias ONLY
void callmeth_mbias( string &seq, string &qual, int pos, string &g, meth *mb, int cycle, bool rev ) {
	unsigned int rs = seq.size();
	unsigned int os = rs - 1;
	char c1, c2;
	for( unsigned int i=0, j=pos; i!=rs; ++i, ++j) {
		if( qual[i] < MIN_BASEQUAL_SCORE )
			continue;

		c1 = g[j];
		c2 = g[j+1];

		if( (c1!='C' && c1!='c') || (c2!='G' && c2!='g') )	// not a CpG site
			continue;

		int k = (rev) ? (os-i) : i;
		if( k < cycle ) {
			if( seq[i] == 'C' ) {
				++ mb[k].C;
			} else if( seq[i] == 'T' ) {
				++ mb[k].T;
			} else {
				++ mb[k].Z;
			}
		}
	}
}

// write meth call into file
void write_methcall( unordered_map<int, meth> &m, const char *pre, const char *suf ) {
	string outfile = pre;
	outfile += suf;
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "ERROR: write output file " << outfile << " failed.\n";
		exit(20);
	}

	unordered_map<int, meth> :: iterator it;
	fout << "#Locus\tC\tT\tZ\n";
	for( it=m.begin(); it!=m.end(); ++it ) {
		fout << it->first << '\t' << it->second.C << '\t' << it->second.T << '\t' << it->second.Z << '\n';
	}
	fout.close();
}

void write_mbias( meth* m, int cycle, const char *pre, const char *suf ) {
	string outfile = pre;
	outfile += suf;
	ofstream fout( outfile.c_str() );
	if( fout.fail() ) {
		cerr << "ERROR: write output file " << outfile << " failed.\n";
		exit(20);
	}

	fout << "#Locus\tC\tT\tZ\n";
	for( unsigned int i=0; i!=cycle; ++i ) {
//		if( m[i].C || m[i].T || m[i].Z ) {
		fout << i+1 << '\t' << m[i].C << '\t' << m[i].T << '\t' << m[i].Z << '\n';
//		}
	}
	fout.close();
}

