#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

/*
 * Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
 * Date  : May 6, 2020
 *
**/

using namespace std;

int main( int argc, char *argv[] ) {
	if( argc != 7 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <in.bed> <wig.header=y|n> <normalization=y|n> <bin.size> <output.file>\n"
			 << "\nThis program is designed to transfer the bed file into a wig file."
			 << "\nThe chr.info file could be found in the index folder for your genome. Note that chrM will be excluded.\n\n";
		exit( 1 );
	}
	unsigned int bin = atoi( argv[5] );
	if( bin == 0 ) {
		cerr << "Error: parameter bin is not a number!\n\n";
		return 127;
	}
	unsigned char head = 'N';
	if( argv[3][0] == 'y' || argv[3][0] == 'Y' )
		head = 'Y';

	unsigned char norm = 'Y';
	if( argv[4][0] == 'n' || argv[4][0] == 'N' )
		norm = 'N';

	cerr << "Chromosome info: " << argv[1] << '\n';
	cerr << "Input bed file : " << argv[2] << '\n';
	cerr << "Write wig head : " << head << '\n';
	cerr << "Normalize data : " << norm << '\n';
	cerr << "Output bin     : " << bin << '\n';
	cerr << "Output file    : " << argv[6] << "\n\n";

	map< string, unsigned int  > info;
	map< string, unsigned int  > chrSize;
	map< string, unsigned int *>  wig;

	ifstream fin;

	cerr << "Reading info file ...\n";
	fin.open( argv[1], ios::in );
	if( fin.fail() ) {
		cerr << "Open file " << argv[1] << " failed!\n\n";
		return 2;
	}
	string chr;
	unsigned int len;
	string line;
	stringstream ss;
	unsigned int * p;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;
		ss.str( line );
		ss.clear();
		ss >> chr >> len;
		if( chr[3] == 'M' )
			continue;

//		cerr << chr << "\t" << len << "\t" << (len/bin + 1 ) << "\n";
		if( info.find( chr ) != info.end() )continue;
		p = new unsigned int[ len/bin + 1 ];
		if( p == NULL ) {
			cerr << "Malloc failed!\n";
			return 5;
		}

		chrSize.insert( pair<string, unsigned int>( chr, len ) );

		info.insert( pair<string, unsigned int>( chr, len/bin+1 ) );
		wig. insert( pair<string, unsigned int *>( chr, p ) );
	}
	fin.close();

	cerr << "Totally " << info.size() << " chromsomes loaded.\n";

	fin.open( argv[2], ios::in );
	if( fin.fail() ) {
		cerr << "Open bed file " << argv[2] << " failed!\n\n";
		return 10;
	}
	register unsigned int pos, start, end;
	register unsigned int count = 0;
	map<string , unsigned int *> :: iterator it;
	cerr << "Reading data ...\n";
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;

		++ count;
		if( !(count & 0xfffff ) )
			cerr << '\r' << count << " lines loaded.";

		ss.str( line );
		ss.clear();
		ss >> chr >> start >> end;
//		cout << chr << "\t" << pos << endl;
		it = wig.find( chr );
		if( it == wig.end() ) {
			cerr << "Warning: no such chrosome '" << line << "', skip\n";
			continue;
		}
		if( end/bin >= info.find(chr)->second ) {
			// ++drop;
			continue;
		}
		for( pos=start+1; pos<=end; ++pos )
			++ (it->second)[ pos/bin ];
	}
	fin.close();
	cerr << "\rDone: " << count << " lines loaded.\n";

	register double weight;
	if( norm == 'Y' ) {
		weight = 10000000.0 / count;
		cerr << "Totally " << count << " reads loaded, weight=" << weight << '\n';
	} else {
		weight = 1.0;
	}

	ofstream fout;
	fout.open( argv[6] );
	if( fout.fail() ) {
		cerr << "Error: write file " << argv[6] << " failed!\n";
		return 20;
	}

	if( head == 'Y' ) {
		fout<< "track type='wiggle_0' name='" << argv[2] <<"' description='bed2wig, file=" << argv[2]
			<< "Normalization=" << norm << ", bin=" << bin << "'\n";
	}

	map<string , unsigned int> :: iterator infoit;
//	fout << "track type=wiggle_0 name=\"" << argv[2] << "_bin" << bin << "\" description=\"" << argv[2] << "_bin" << bin << "\"\n";
	fout.setf( ios::fixed );
	fout.precision(4);
	for( infoit=info.begin(); infoit!=info.end(); ++infoit ) {
		cerr << "\rProcessing " << infoit->first << " ...";
		fout << "variableStep chrom=" << infoit->first << " span=" << bin << "\n";
		it = wig.find( infoit->first );
		if( it == wig.end() ) {
			cerr << "\rWarning: cannot find data for " << infoit->first << " , skip!";
			continue;
		}
		register unsigned int i, j=1, k;
		k = infoit->second;
		unsigned int chr_size = chrSize.find(infoit->first)->second;
		for( i=0; i!=k; ++i ) {
			if( j > chr_size )
				j = chr_size;

			if( (it->second)[i] )
				fout << j << '\t' << (it->second)[i] * weight << '\n';

			j += bin;
		}
		delete [] it->second;
	}
	fout.close();

	cerr << "\rDone. Please check the output file '" << argv[6] << "'.\n\n";

	return 0;
}


