#include <string>
#include <vector>
#include <iostream>
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package
 * Date: Jul 2021
*/

// load genome from multi-fasta
void loadgenome( const char * file, unordered_map<string, string> & genome ) {
	ifstream fin( file );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << file << " !\n";
		exit(200);
	}
	cout << "Loading genome: " << file << '\n';
	string line, chr, tmp="X";	// X is for position-taking
	register unsigned int i, j;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;
		j = line.length();
		if( line[0] == '>' ) {
			chr.clear();
			for( i=1; i!=j; ++i ) {
				if( line[i]==' ' || line[i]=='\t' )break;
				chr += line[i];
			}
			genome.insert( pair<string, string>(chr, tmp) );
		} else {
			genome[chr] += line;
		}
	}
	fin.close();
}

// load sequence from 1 chr
void loadchr( const char * file, string & seq ) {
	ifstream fin( file );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << file << "!\n";
		exit(200);
	}

	seq = 'X';	// X is for position-taking
	string line;
	register unsigned int i, j;
	getline( fin, line );	// drop 1st line of chrXXX
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		seq += line;
	}
	seq += 'Y';	// end marker
	fin.close();
//	cerr << "Load fastq done: length=" << seq.length()-2 << "\n";
}

// fix cigar
int get_readLen_from_cigar( const string &cigar ) {
	register int i, j;
	register int size = 0;
	register int cs = cigar.size();
	register const char * p = cigar.c_str();

	for(i=0, j=0; i!=cs; ++i) {
		if( p[i] <= '9' ) {   // digital
			j *= 10;
			j += p[i] - '0';
		} else {        // MUST be M, I, or D
			if( p[i] == 'M' ) { // match or mismatch, keep
				size += j;
			} else if ( p[i] == 'I' ) { // insertion, ignore
				// do nothing
			} else if ( p[i] == 'D' ) { // deletion, add place holders
				size += j;
			} else {	// unsupported CIGAR element
				return 0;
			}
			j = 0;
		}
	}

	return size;
}

// update sequence and quality using CIGAR information to support indels
bool fix_cigar( string &cigar, string &realSEQ, string &realQUAL, string &seq, string &qual ) {
	register int i, j, k;
	j = 0;
	register int curr = 0;
	register int cs = cigar.size();
	realSEQ.clear();
	realQUAL.clear();
	for(i=0; i!=cs; ++i) {
		if( cigar[i] <= '9' ) {   // digital
			j *= 10;
			j += cigar[i] - '0';
		} else {	// MUST be M, I, or D
			if( cigar[i] == 'M' ) { // match or mismatch, copy seq and qual
				for(k=0; k!=j; ++k) {
					realSEQ  +=  seq[ curr+k ];
					realQUAL += qual[ curr+k ];
				}
				curr += j;
			} else if ( cigar[i] == 'I' ) { // insertion, discard this part
				curr += j;
			} else if ( cigar[i] == 'D' ) { // deletion, add place holders
				for(k=0; k!=j; ++k) {
					realSEQ  +=  'N';
					realQUAL += '\0';
				}
			} else {	// unsupported CIGAR element
				return false;
			}
			j = 0;
		}
	}
	return true;
}

// revert cigar: usually xxMxIxDxxM
void revert_cigar(string &raw, string &rev, vector<int> &seg ) {
	seg.clear();
	seg.push_back( -1 );

	int len = raw.length();
	for( int i=1; i!=len; ++i ) {
		if( raw[i] > '9' ) {
			seg.push_back( i );
		}
	}

	rev.clear();
	int j = seg.size()-2;
	for( ; j>=0; --j ) {
		for( int k=seg[j]+1; k<=seg[j+1]; ++k ) {
			rev += raw[k];
		}
	}
}

// revert MD tag
void revert_MDtag(string &raw, string &rev, vector<int> &seg ) {
	seg.clear();

	bool flag;	// whether previous letter is numeric
	int i = ( raw[0] == '0' ) ? 1 : 0;
	if( raw[i] > '9' ) {
		switch( raw[i] ) {
			case 'A': raw[i]='T';break;
			case 'C': raw[i]='G';break;
			case 'G': raw[i]='C';break;
			case 'T': raw[i]='A';break;
			default: break;
		}
		flag = false;
	} else {
		flag = true;
	}
//	cerr << "insert " << i << '\n';
	seg.push_back( i );
	++ i;

	int len = raw.length();
	bool curr;
	for( ; i!=len; ++i ) {
		if( raw[i] > '9' ) {
			switch( raw[i] ) {
				case 'A': raw[i]='T';break;
				case 'C': raw[i]='G';break;
				case 'G': raw[i]='C';break;
				case 'T': raw[i]='A';break;
				default: break;
			}
			curr = false;
		} else {
			curr = true;
		}
		if( curr != flag ) {
//			cerr << "insert " << i << '\n';
			seg.push_back( i );
		}
		flag = curr;
	}
	seg.push_back( len );

	rev.clear();
	if( ! flag ) rev += '0';	// the last letter is NOT numeric, it's a mismatch
	int j = seg.size()-2;
	for( ; j>=0; --j ) {
		for( int k=seg[j]; k!=seg[j+1]; ++k ) {
			rev += raw[k];
		}
	}
}

// segment SAM record
void parseSAM_mode3( char * psam, samRecord & read ) {
	register unsigned int i = 1;
	read.seqName = 0;
	while( true ) {	//flag
		if( psam[i]=='\t' ) {++ i;read.flag = i;break;}
		++ i;
	}
	while( true ) {	//chr
		if( psam[i]=='\t' ) {++ i;read.chr = i;break;}
		++ i;
	}
	while( true ) {	//pos
		if( psam[i]=='\t' ) {++ i;read.pos = i;break;}
		++ i;
	}
	while( true ) {	//score
		if( psam[i]=='\t' ) {++ i;read.score = i;break;}
		++ i;
	}
	while( true ) {	//cigar
		if( psam[i]=='\t' ) {++ i;read.cigar = i;break;}
		++ i;
	}
	while( true ) {	//mateflag
		if( psam[i]=='\t' ) {++ i;read.mateflag = i;break;}
		++ i;
	}
	while( true ) {	//matepos
		if( psam[i]=='\t' ) {++ i;read.matepos = i;break;}
		++ i;
	}
	while( true ) {	//matedist
		if( psam[i]=='\t' ) {++ i;read.matedist = i;break;}
		++ i;
	}
	while( true ) {	//seq
		if( psam[i]=='\t' ) {++ i;read.seq = i;break;}
		++ i;
	}
	while( true ) {	//qual
		if( psam[i]=='\t' ) {++ i;read.qual = i;break;}
		++ i;
	}
	while( true ) {	//remaining
		if( psam[i]=='\t' ) {++ i;read.remaining = i;break;}
		++ i;
	}
}

void parseSAM_mode4( char * psam, samRecord & read ) {
	register unsigned int i = 1;
	read.seqName = 0;
	while( true ) {	//flag
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.flag = i;break;}
		++ i;
	}
	while( true ) {	//chr
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.chr = i;break;}
		++ i;
	}
	while( true ) {	//pos
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.pos = i;break;}
		++ i;
	}
	while( true ) {	//score
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.score = i;break;}
		++ i;
	}
	while( true ) {	//cigar
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.cigar = i;break;}
		++ i;
	}
	while( true ) {	//mateflag
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.mateflag = i;break;}
		++ i;
	}
	while( true ) {	//matepos
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.matepos = i;break;}
		++ i;
	}
	while( true ) {	//matedist
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.matedist = i;break;}
		++ i;
	}
	while( true ) {	//seq
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.seq = i;break;}
		++ i;
	}
	while( true ) {	//qual
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.qual = i;break;}
		++ i;
	}
	while( true ) {	//remaining
		if( psam[i]=='\t' ) {psam[i]='\0';++ i;read.remaining = i;break;}
		++ i;
	}
}

