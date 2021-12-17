#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <zlib.h>
#include "common.h"
#include "util.h"

using namespace std;

// hisat2 supports 256 character long of read id, and does not has --sam-no-qname-trunc option
const unsigned int MAX_CONVERTED_READ_ID = 256;

// changes in v2.1: use 2 threads for file loading; change Phred64 to Phred33 when necessary

/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package, adapted from Ktrim
 * Date: Jul 2021
 *
 * The conversion log's rule: NOTE that LINE_NUMBER is NO longer recorded in Msuite2
 *   read 1:
 *	   @S|xx;xx$ => for reads with endC and C>T changes, 'S' is its quality score
 *	   @xx;xx$	 => for reads without endC
 *	   @$		 => for reads without endC and conversions
 *   read 2:
 *	   @S|xx;xx$ => for reads with frontG and G>A changes, 'S' is its quality score
 *	   @xx;xx$	 => for reads without frontG
 *	   @$		 => for reads without frontG and conversions
**/

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * so use 4 is not good
*/
bool check_mismatch_dynamic_PE( const string & s1, string & s2, unsigned int pos, const adapter_info* ai ) {
	register unsigned int mis1=0, mis2=0;
	register unsigned int i, len;
	len = s1.length() - pos;
	if( len > ai->adapter_len )
		len = ai->adapter_len;

	register unsigned int max_mismatch_dynamic = len >> 3;
	// here use 1/8 of the total length, use roof(len/8)

	if( (max_mismatch_dynamic<<3) != len )
		++ max_mismatch_dynamic;

	// check mismatch for each read
	const char * p = s1.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != ai->adapter_r1[i] ) {
			++ mis1;
			if( mis1 > max_mismatch_dynamic )
				return false;
		}
	}
	p = s2.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != ai->adapter_r2[i] ) {
			++ mis2;
			if( mis2 > max_mismatch_dynamic )
				return false;
		}
	}
	// check total mismatch
	max_mismatch_dynamic = (len+1) >> 2;
	if( mis1 + mis2 > max_mismatch_dynamic )
		return false;

	return true;
}

bool inline is_revcomp( const char a, const char b ) {
	switch( a ) {
		case 'A': return b=='T';
		case 'C': return b=='G';
		case 'G': return b=='C';
		case 'T': return b=='A';
		default : return false;
	}
}

unsigned int load_read_batch_gz( gzFile &gfp, char * gz_buffer, string *id, string *seq, string *qual ) {
	unsigned int loaded = 0;
	for( ; loaded!=READS_PER_BATCH; ++loaded ) {
		if( gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE ) == NULL ) break;
		id [ loaded ] = gz_buffer;
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );
		seq [ loaded ] = gz_buffer;
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );	// this line is useless
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );
		qual[ loaded ] = gz_buffer;

		id [ loaded ].pop_back();	// trim the tail '\n'
		seq[ loaded ].pop_back();
		if( qual[loaded].size() != seq[loaded].size() )	//the last read may not contain '\n' for quality line
			qual[loaded].pop_back();
	}
	return loaded;
}

unsigned int load_read_batch( ifstream &fq, string *id, string *seq, string *qual ) {
	unsigned int loaded = 0;
	for( ; loaded!=READS_PER_BATCH; ++loaded ) {
		getline( fq, id  [ loaded ] );
		if( fq.eof() )break;
		getline( fq, seq [ loaded ] );
		getline( fq, qual[ loaded ] );	// this line is useless
		getline( fq, qual[ loaded ] );
	}
	return loaded;
}

int main( int argc, const char *argv[] ) {
	if( argc < 5 ) {
		cerr << "\nUsage: " << argv[0] << " <r1.fq> <r2.fq> <cycle> <out.prefix> "
			 << "[mode] [thread] [min.length] [min.quality] [library] "
			 << "[cut.head.r1] [cut.tail.r1] [cut.head.r2] [cut.tail.r2]\n\n"

			 << "This program is part of Msuite and is designed to do fastq statistics, quality-trimming,\n"
			 << "adapter-trimming and C->T/G->A conversions for Paired-End reads generated by illumina sequencers.\n\n"

			 << "Run modes:\n"
			 << "  0. Do not perform any conversions. Usable for normal DNA/RNA-seq alignment.\n"
			 << "  3. Do C>T (G>A) for all C (G) sites in read1 (read2). Suitable for general BS-seq data.\n"
			 << "  4. Do C>T (G>A) for C (G) in CpG sites in read1 (read2). Suitable for novel TAPS-like data.\n\n"

			 << "Default parameters:\n"
			 << "  mode: 0 (could be 0,3,4)\n"
			 << "  thread: 1\n"
			 << "  min.length: 36\n"
			 << "  min.quality: 53 (33+20 for phred33('!', or '#') scoring system)\n"
			 << "    Phred64 to Phred33 conversion is automatically ON if min.quality >= 74\n"
			 << "  library: illumina\n"
			 << "  cut.head.r1: 0\n"
			 << "  cut.tail.r1: 0\n"
			 << "  cut.head.r2: 0\n"
			 << "  cut.tail.r2: 0\n\n";

		return 2;
	}

	int mode = 0;
	int thread = 1;
	int min_length = 36;
	char quality = 53;
	const char *libraryKit = "illumina";
	int cut_head_r1=0, cut_tail_r1=0, cut_head_r2=0, cut_tail_r2=0;
	const adapter_info* ai;
	unsigned int cycle = atoi( argv[3] );
	if( cycle == 0 ) {
		cerr << "Error: Unacceptable cycle!\n";
		exit(1);
	}
	if( argc > 5 ) {
		mode = atoi( argv[5] );
		if( argc > 6) {
			thread = atoi( argv[6] );
			if( argc > 7 ) {
				min_length = atoi( argv[7] );
				if( argc > 8 ) {
					quality = (unsigned char) atoi( argv[8] );
					if( argc > 9 ) {
						libraryKit = argv[9];
						if( argc > 10 ) {
							cut_head_r1 = atoi( argv[10] );
							if( argc > 11 ) {
								cut_tail_r1 = atoi( argv[11] );
								if( argc > 12 ) {
									cut_head_r2 = atoi( argv[12] );
									if( argc > 13 ) {
										cut_tail_r2 = atoi( argv[13] );
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if( mode != 0 && mode !=3 && mode != 4 ) {
		cerr << "Error: invalid run mode! Must be 0, 3, or 4!\n";
		return 100;
	}
	if( thread == 0 ) {
		cerr << "Warning: thread is set to 0! I will use all threads instead.\n";
		thread = omp_get_max_threads();
	}
	if( min_length == 0 ) {
		cerr << "Error: invalid min_length! Must be a positive number!\n";
		return 101;
	}
	if( quality == 0 ) {
		cerr << "Error: invalid quality! Must be a positive number!\n";
		return 102;
	}
	if( strcmp(libraryKit, "illumina")==0 || strcmp(libraryKit, "Illumina")==0 ) {
		ai = &illumina_adapter;
	} else if ( strcmp(libraryKit, "nextera")==0 || strcmp(libraryKit, "Nextera")==0 ) {
		ai = &nextera_adapter;
	} else if ( strcmp(libraryKit, "bgi")==0 || strcmp(libraryKit, "BGI")==0 ) {
		ai = &bgi_adapter;
	} else {
		cerr << "Error: invalid library kit! Currently only supports illumina, nextera, and bgi!\n";
		return 103;
	}

	bool changePhred = false;
	if( quality >= 74 ) {	// consider it is Phred 64, otherwise it means Phred>=40 which is impossible
		changePhred = true;
	}

	string *id1   = new string [READS_PER_BATCH];
	string *id2   = new string [READS_PER_BATCH];
	string *seq1  = new string [READS_PER_BATCH];
	string *seq2  = new string [READS_PER_BATCH];
	string *qual1 = new string [READS_PER_BATCH];
	string *qual2 = new string [READS_PER_BATCH];
	string unk;

	int *dropped	  = new int [thread];
	int *real_adapter = new int [thread];
	int *tail_adapter = new int [thread];

	fastqstat * AllR1stat = new fastqstat[ cycle ];
	memset( AllR1stat, 0, cycle*sizeof(fastqstat) );
	fastqstat * AllR2stat = new fastqstat[ cycle ];
	memset( AllR2stat, 0, cycle*sizeof(fastqstat) );

	fastqstat * AllR1stat_trimmed = new fastqstat[ cycle ];
	memset( AllR1stat_trimmed, 0, cycle*sizeof(fastqstat) );
	fastqstat * AllR2stat_trimmed = new fastqstat[ cycle ];
	memset( AllR2stat_trimmed, 0, cycle*sizeof(fastqstat) );

	// buffer for storing the modified reads per thread
	char ** buffer1 = new char * [thread];
	char ** buffer2 = new char * [thread];
	int  * b1stored = new int    [thread];
	int  * b2stored = new int    [thread];
	fastqstat **R1stat = new fastqstat * [thread];
	fastqstat **R2stat = new fastqstat * [thread];
	// new in v2: report the fqstat after trimming
	fastqstat **R1stat_trimmed = new fastqstat * [thread];
	fastqstat **R2stat_trimmed = new fastqstat * [thread];

	for(unsigned int i=0; i!=thread; ++i) {
		buffer1[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];
		buffer2[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];

		R1stat[i]  = new fastqstat [ cycle ];
		R2stat[i]  = new fastqstat [ cycle ];

		R1stat_trimmed[i]  = new fastqstat [ cycle ];
		R2stat_trimmed[i]  = new fastqstat [ cycle ];

		dropped[i] = 0;
		real_adapter[i] = 0;
		tail_adapter[i] = 0;
	}

	cerr << "Loading files ...\n";
	// deal with multiple input files
	vector<string> R1s, R2s;
	string fileName="";
	for(unsigned int i=0; argv[1][i]!='\0'; ++i) {
		if( argv[1][i] == FILE_SEPARATOR ) {
			R1s.push_back( fileName );
			fileName.clear();
		} else {
			fileName += argv[1][i];
		}
	}
	R1s.push_back( fileName );

	fileName.clear();
	for(unsigned int i=0; argv[2][i]!='\0'; ++i) {
		if( argv[2][i] == FILE_SEPARATOR ) {
			R2s.push_back( fileName );
			fileName.clear();
		} else {
			fileName += argv[2][i];
		}
	}
	R2s.push_back( fileName );
	
	if( R1s.size() != R2s.size() ) {
		cerr << "Fatal error: Read1 and Read2 do not contain equal sized files!\n";
		return 10;
	}
	unsigned int totalFiles = R1s.size();
	cout << "INFO: " << totalFiles << " paired fastq files will be loaded.\n";

	string base = argv[4];
	ofstream fout1( (base+".R1.fq").c_str() ), fout2( (base+".R2.fq").c_str() );
	if( fout1.fail() || fout2.fail() ) {
		cout << "Error: write file failed!\n";
		fout1.close();
		fout2.close();
		return 3;
	}
	
	ifstream fq1, fq2;
	gzFile gfp1, gfp2;
	char * gz_buffer_1 = new char [ MAX_SEQNAME_SIZE ];
	char * gz_buffer_2 = new char [ MAX_SEQNAME_SIZE ];
	register int totalReads = 0;
	for( int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		register int i_file = R1s[fileCnt].size() - 3;
		register const char * p_file = R1s[fileCnt].c_str();
		if( p_file[i_file]=='.' && p_file[i_file+1]=='g' && p_file[i_file+2]=='z' ) {	// .gz file
			file_is_gz = true;
			gfp1 = gzopen( p_file, "r" );
			gfp2 = gzopen( R2s[fileCnt].c_str(), "r" );
			if( gfp1==NULL || gfp2==NULL ) {
				cerr << "Error: open fastq file failed!\n";
				fout1.close();
				fout2.close();
				return 11;
			}
		} else {	// plain text
			file_is_gz = false;
			fq1.open( R1s[fileCnt] );
			fq2.open( R2s[fileCnt] );
			if( fq1.fail() || fq2.fail() ) {
				cerr << "Error: open fastq file failed!\n";
				fout1.close();
				fout2.close();
				return 11;
			}
		}
		// load and process reads
		while( true ) {
			unsigned int loaded = 0;
			//TODO: use 2 threads to load data
			unsigned int loaded_2 = 0;

			if( file_is_gz ) {
				if( thread == 1 ) {
					// read 1
					loaded = load_read_batch_gz( gfp1, gz_buffer_1, id1, seq1, qual1 );
					// read 2
					loaded_2 = load_read_batch_gz( gfp2, gz_buffer_2, id2, seq2, qual2 );
				} else {
					omp_set_num_threads( 2 );
					#pragma omp parallel
					{
						unsigned int tn = omp_get_thread_num();
						if( tn == 0 ) {
							loaded = load_read_batch_gz( gfp1, gz_buffer_1, id1, seq1, qual1 );
						} else {
							loaded_2 = load_read_batch_gz( gfp2, gz_buffer_2, id2, seq2, qual2 );
						}
					}
				}
			} else {
				if( thread == 1 ) {
					loaded = load_read_batch( fq1, id1, seq1, qual1 );
					loaded_2 = load_read_batch( fq2, id2, seq2, qual2 );
				} else {
					omp_set_num_threads( 2 );
					#pragma omp parallel
					{
						unsigned int tn = omp_get_thread_num();
						if( tn == 0 ) {
							loaded = load_read_batch( fq1, id1, seq1, qual1 );
						} else {
							loaded_2 = load_read_batch( fq2, id2, seq2, qual2 );
						}
					}
				}
			}

			if( loaded != loaded_2 ) {	// error happens
				cerr << "ERROR in loading file!\n";
				exit(10);
			}

			if( loaded == 0 )	// reach the end of file
				break;

			// start parallalization
			omp_set_num_threads( thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				unsigned int start = loaded * tn / thread;
				unsigned int end   = loaded * (tn+1) / thread;

				// normalization
				b1stored[tn] = 0;
				b2stored[tn] = 0;
				memset( R1stat[tn], 0, cycle*sizeof(fastqstat) );
				memset( R2stat[tn], 0, cycle*sizeof(fastqstat) );
				memset( R1stat_trimmed[tn], 0, cycle*sizeof(fastqstat) );
				memset( R2stat_trimmed[tn], 0, cycle*sizeof(fastqstat) );
			
				string conversionLog1, conversionLog2;
				register int i, j;
				register int last_seed;
				vector<int> seed;
				vector<int> :: iterator it;
				const char *p, *q;
				char *conversion = new char [MAX_CONVERSION];
				char numstr[10]; // enough to hold all numbers up to 99,999,999 plus ':'

				for( register int ii=start; ii!=end; ++ii ) {
					// check R1/R2 cycles
					if( seq2[ii].size() != seq1[ii].size() ) {
						if( seq2[ii].size() > seq1[ii].size() ) {
							seq2[ii].resize(  seq1[ii].size()  );
							qual2[ii].resize( qual1[ii].size() );
						} else {
							seq1[ii].resize(  seq2[ii].size()  );
							qual1[ii].resize( qual2[ii].size() );
						}
					}

					//if the reads are longer than "cycle" paramater, only keep the head "cycle" ones
					if( seq1[ii].size() > cycle ) {
						seq1[ii].resize(  cycle );
						qual1[ii].resize( cycle );
						seq2[ii].resize(  cycle );
						qual2[ii].resize( cycle );
					}

					// raw fqstatistics
					p = seq1[ii].c_str();
					q = seq2[ii].c_str();

					j = seq1[ii].size();
					for( i=0; i!=j; ++i ) {
						switch ( p[i] ) {
							case 'a':
							case 'A': R1stat[tn][i].A ++; break;
							case 'c':
							case 'C': R1stat[tn][i].C ++; break;
							case 'g':
							case 'G': R1stat[tn][i].G ++; break;
							case 't':
							case 'T': R1stat[tn][i].T ++; break;
							default : R1stat[tn][i].N ++; break;
						}
					}
					j = seq2[ii].size();
					for( i=0; i!=j; ++i ) {
						switch ( q[i] ) {
							case 'a':
							case 'A': R2stat[tn][i].A ++; break;
							case 'c':
							case 'C': R2stat[tn][i].C ++; break;
							case 'g':
							case 'G': R2stat[tn][i].G ++; break;
							case 't':
							case 'T': R2stat[tn][i].T ++; break;
							default : R2stat[tn][i].N ++; break;
						}
					}

					// quality control
					p = qual1[ii].c_str();
					q = qual2[ii].c_str();
					i = get_quality_trim_cycle_pe( p, q, qual1[ii].size(), min_length, quality );

					if( i < min_length ) { // not long enough
						++ dropped[ tn ];
						continue;
					}
					seq1[ii].resize(  i );
					seq2[ii].resize(  i );
					qual1[ii].resize( i );
					qual2[ii].resize( i );
					if( changePhred ) {
						for( j=0; j!=i; ++j ) {
							qual1[ii][j] -= 31;
							qual2[ii][j] -= 31;
						}
					}

					// looking for seed target, 1 mismatch is allowed for these 2 seeds
					// which means seq1 and seq2 at least should take 1 perfect seed match
					seed.clear();
					for( i=0; (i=seq1[ii].find(ai->adapter_index, i)) != string::npos; ++i )
						seed.push_back( i );
					for( i=0; (i=seq2[ii].find(ai->adapter_index, i)) != string::npos; ++i )
						seed.push_back( i );

					sort( seed.begin(), seed.end() );

					last_seed = impossible_seed;	// a position which cannot be in seed
					for( it=seed.begin(); it!=seed.end(); ++it ) {
						if( *it != last_seed ) {
						// as there maybe the same value in seq1_seed and seq2_seed,
						// use this to avoid re-calculate that pos
							if( check_mismatch_dynamic_PE( seq1[ii], seq2[ii], *it, ai) )
								break;
							last_seed = *it;
						}
					}

					if( it != seed.end() ) {	// adapter found
						++ real_adapter[tn];
						if( *it >= min_length )	{
							seq1[ii].resize(  *it );
							seq2[ii].resize(  *it );
							qual1[ii].resize( *it );
							qual2[ii].resize( *it );
						} else {	// drop this read as its length is not enough
							++ dropped[tn];
							continue;
						}
					} else {	// seed not found, now check the tail, if perfect match, trim the tail
						i = seq1[ii].length() - 2;
						p = seq1[ii].c_str();
						q = seq2[ii].c_str();
						if( p[i]==ai->adapter_r1[0] && p[i+1]==ai->adapter_r1[1] &&
									q[i]==ai->adapter_r2[0] && q[i+1]==ai->adapter_r2[1] ) {
							// if it is a real adapter, then Read1 and Read2 should be complimentary
							// in real data, the heading 5 bp are usually of poor quality, therefore we test the 6th, 7th
							if( is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) ) {
								if( i < min_length ) {
									++ dropped[tn];
									continue;
								}
								seq1[ii].resize( i );
								seq2[ii].resize( i );
								qual1[ii].resize( i );
								qual2[ii].resize( i );

								++ tail_adapter[tn];
							}
						} else {	// tail 2 is not good, check tail 1
							++ i;
							if( p[i] == ai->adapter_r1[0] && q[i] == ai->adapter_r2[0] ) {
								if(is_revcomp(p[5], q[i-6]) && is_revcomp(q[5], p[i-6]) &&
										is_revcomp(p[6], q[i-7]) && is_revcomp(q[6], p[i-7]) ) {
									if( i < min_length ) {
										++ dropped[tn];
										continue;
									}
									seq1[ii].resize(  i );
									seq2[ii].resize(  i );
									qual1[ii].resize( i );
									qual2[ii].resize( i );

									++ tail_adapter[tn];
								}
							}
						}
					}

					// cut head and tail
					if( cut_head_r1 ) {
						seq1[ii].erase(  0, cut_head_r1 );
						qual1[ii].erase( 0, cut_head_r1 );
					}
					if( cut_tail_r1 ) {
						seq1[ii].resize(  seq1[ii].length()  - cut_tail_r1 );
						qual1[ii].resize( qual1[ii].length() - cut_tail_r1 );
					}

					if( seq1[ii].length() < min_length ) {
						++ dropped[tn];
						continue;
					}

					if( cut_head_r2 ) {
						seq2[ii].erase(  0, cut_head_r2 );
						qual2[ii].erase( 0, cut_head_r2 );
					}
					if( cut_tail_r2 ) {
						seq2[ii].resize(  seq2[ii].length()  - cut_tail_r2 );
						qual2[ii].resize( qual2[ii].length() - cut_tail_r2 );
					}

					//fqstatistics after trimming
					p = seq1[ii].c_str();
					q = seq2[ii].c_str();
					j = seq1[ii].size();
					for( i=0; i!=j; ++i ) {
						switch ( p[i] ) {
							case 'a':
							case 'A': R1stat_trimmed[tn][i].A ++; break;
							case 'c':
							case 'C': R1stat_trimmed[tn][i].C ++; break;
							case 'g':
							case 'G': R1stat_trimmed[tn][i].G ++; break;
							case 't':
							case 'T': R1stat_trimmed[tn][i].T ++; break;
							default : R1stat_trimmed[tn][i].N ++; break;
						}
					}
					j = seq2[ii].size();
					for( i=0; i!=j; ++i ) {
						switch ( q[i] ) {
							case 'a':
							case 'A': R2stat_trimmed[tn][i].A ++; break;
							case 'c':
							case 'C': R2stat_trimmed[tn][i].C ++; break;
							case 'g':
							case 'G': R2stat_trimmed[tn][i].G ++; break;
							case 't':
							case 'T': R2stat_trimmed[tn][i].T ++; break;
							default : R2stat_trimmed[tn][i].N ++; break;
						}
					}

					//check if there is any white space in the IDs; if so, remove all the data after the whitespace
					j = id1[ii].size();
					p = id1[ii].c_str();
					for( i=1; i!=j; ++i ) {
						if( p[i]==' ' || p[i]=='\t' ) {	// white space, then trim ID
							id1[ii].resize( i );
							break;
						}
					}
					j = id2[ii].size();
					q = id2[ii].c_str();
					for( i=0; i!=j; ++i ) {
						if( q[i]==' ' || q[i]=='\t' ) {	// white space, then trim ID
							id2[ii].resize( i );
							break;
						}
					}

					// do C->T and G->A conversion
					if( mode == 3 ) {	// in the current implementation, id1 and id2 are different!!!
						// in mode 3, there is NO endC and frontG issues
						id1[ii][0] = CONVERSION_LOG_END;
						j = seq1[ii].size();	// seq1 and seq2 are of the same size
						conversionLog1 = NORMAL_SEQNAME_START;
						for( i=0; i!=j; ++i ) {
							if( seq1[ii][i] == 'C' ) {
								seq1[ii][i] = 'T';
								sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
								conversionLog1 += numstr;
							}
						}
						if( conversionLog1.back() == CONVERSION_LOG_SEPARATOR )
							conversionLog1.pop_back();

						if( conversionLog1.length() > MAX_CONVERTED_READ_ID ) {
							//cerr << "LONG read ID!\n";
							++ dropped[tn];
							continue;
						}
						/*fout1 << NORMAL_SEQNAME_START << line << conversionLog << id1 << '\n'
								<< seq1 << "\n+\n" << qual1 << '\n';*/

						id2[ii][0] = CONVERSION_LOG_END;
						conversionLog2 = NORMAL_SEQNAME_START;	// read2 does not record line number
						for( i=0; i!=j; ++i ) {
							if( seq2[ii][i] == 'G' ) {
								seq2[ii][i] = 'A';
								sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
								conversionLog2 += numstr;
							}
						}
						if( conversionLog2.back() == CONVERSION_LOG_SEPARATOR )
							conversionLog2.pop_back();

						if( conversionLog2.length() > MAX_CONVERTED_READ_ID ) {
							//cerr << "LONG read ID!\n";
							++ dropped[tn];
							continue;
						}

						b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%s%s\n%s\n+\n%s\n",
												conversionLog1.c_str(), id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );
						b2stored[tn] += sprintf( buffer2[tn]+b2stored[tn], "%s%s\n%s\n+\n%s\n",
												conversionLog2.c_str(), id2[ii].c_str(), seq2[ii].c_str(), qual2[ii].c_str() );
					} else if ( mode == 4 ) {	// this is the major task for EMaligner
						// modify id1 to add line number (to facilitate the removing ambigous step)
						// check seq1 for C>T conversion
						id1[ii][0] = CONVERSION_LOG_END;
						conversionLog1 = NORMAL_SEQNAME_START;
						j = seq1[ii].size()-1;
						if( seq1[ii].back() == 'C' ) { //ther is a 'C' and the end, discard it (but record its Quality score);
							//otherwise it may introduce a mismatch in alignment
							if( qual1[ii].back() == '@' ) {
								conversionLog1 += REPLACEMENT_CHAR_AT;
							} else {
								conversionLog1 += qual1[ii].back();
							}
							conversionLog1 += KEEP_QUAL_MARKER;
							seq1[ii].pop_back();
							qual1[ii].pop_back();
						}
						for( i=0; i!=j; ++i ) {
							if( seq1[ii][i]=='C' && seq1[ii][i+1]=='G' ) {
								seq1[ii][i] = 'T';
								sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
								conversionLog1 += numstr;
							}
						}
						if( conversionLog1.back() == CONVERSION_LOG_SEPARATOR )
							conversionLog1.pop_back();

						if( conversionLog1.length() > MAX_CONVERTED_READ_ID ) {
							//cerr << "LONG read ID!\n";
							++ dropped[tn];
							continue;
						}

						// format for ID1:
						// if there is a C at the end
						//	@ x | C1;C2;C3# raw_seq_name
						//	the | is the marker for the existence of tail 'C' and 'x' is its quality score
						//	if exists, | is ALWAYS two bytes after '+' (use this to test its existence)
						// if there is No C at the end
						//	@C1;C2;C3# raw_seq_name
						//
						// All the numbers in line_number and C1,C2,C3... are HEX

						// check seq2 for G>A conversion
						id2[ii][0] = CONVERSION_LOG_END;
						conversionLog2 = NORMAL_SEQNAME_START;
						if( seq2[ii][0] == 'G' ) { //'G' at the front, discard it (but record its Quality score)
							if( qual2[ii][0] == '@' ) {
								conversionLog2 += REPLACEMENT_CHAR_AT;
							} else {
								conversionLog2 += qual2[ii][0];
							}
							conversionLog2 += KEEP_QUAL_MARKER;
						}
						j = seq2[ii].size();
						for( i=1; i!=j; ++i ) {
							if( seq2[ii][i]=='G' && seq2[ii][i-1]=='C' ) {
								seq2[ii][i] = 'A';
								sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
								conversionLog2 += numstr;
							}
						}
						if( conversionLog2.back() == CONVERSION_LOG_SEPARATOR )
							conversionLog2.pop_back();

						if( conversionLog2.length() > MAX_CONVERTED_READ_ID ) {
							//cerr << "LONG read ID!\n";
							++ dropped[tn];
							continue;
						}

						b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%s%s\n%s\n+\n%s\n",
									conversionLog1.c_str(), id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );

						if( seq2[ii][0] != 'G' ) {
							b2stored[tn] += sprintf( buffer2[tn]+b2stored[tn], "%s%s\n%s\n+\n%s\n",
													conversionLog2.c_str(), id2[ii].c_str(), seq2[ii].c_str(), qual2[ii].c_str() );
						} else {
							p = seq2[ii].c_str();
							q = qual2[ii].c_str();
							b2stored[tn] += sprintf( buffer2[tn]+b2stored[tn], "%s%s\n%s\n+\n%s\n",
													conversionLog2.c_str(), id2[ii].c_str(), p+1, q+1 );
						}
					} else {	// mode 0: no need to do conversion
						b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%s\n%s\n+\n%s\n",
									id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );
						b2stored[tn] += sprintf( buffer2[tn]+b2stored[tn], "%s\n%s\n+\n%s\n",
									id2[ii].c_str(), seq2[ii].c_str(), qual2[ii].c_str() );
					}
				}
			}	// parallel body

			// write output and update fastq statistics
			if( thread > 1 ) {	// multi-thread
				omp_set_num_threads( 2 );
				#pragma omp parallel
				{
					unsigned int tn = omp_get_thread_num();

					if( tn == 0 ) {
						for(register int i=0; i!=thread; ++i ) {
							fout1 << buffer1[i];
						}
					} else {
						for(register int i=0; i!=thread; ++i ) {
							fout2 << buffer2[i];
						}
					}
				}
			} else {	// single-thread
				fout1 << buffer1[0];
				fout2 << buffer2[0];
			}

			for(register int i=0; i!=thread; ++i ) {
				for( register int j=0; j!=cycle; ++j ) {
					AllR1stat[j].A += R1stat[i][j].A;
					AllR1stat[j].C += R1stat[i][j].C;
					AllR1stat[j].G += R1stat[i][j].G;
					AllR1stat[j].T += R1stat[i][j].T;
					AllR1stat[j].N += R1stat[i][j].N;

					AllR2stat[j].A += R2stat[i][j].A;
					AllR2stat[j].C += R2stat[i][j].C;
					AllR2stat[j].G += R2stat[i][j].G;
					AllR2stat[j].T += R2stat[i][j].T;
					AllR2stat[j].N += R2stat[i][j].N;

					AllR1stat_trimmed[j].A += R1stat_trimmed[i][j].A;
					AllR1stat_trimmed[j].C += R1stat_trimmed[i][j].C;
					AllR1stat_trimmed[j].G += R1stat_trimmed[i][j].G;
					AllR1stat_trimmed[j].T += R1stat_trimmed[i][j].T;
					AllR1stat_trimmed[j].N += R1stat_trimmed[i][j].N;

					AllR2stat_trimmed[j].A += R2stat_trimmed[i][j].A;
					AllR2stat_trimmed[j].C += R2stat_trimmed[i][j].C;
					AllR2stat_trimmed[j].G += R2stat_trimmed[i][j].G;
					AllR2stat_trimmed[j].T += R2stat_trimmed[i][j].T;
					AllR2stat_trimmed[j].N += R2stat_trimmed[i][j].N;
				}
			}
			totalReads += loaded;
			cerr << '\r' << totalReads << " reads loaded";

			if( file_is_gz ) {
				if( gzeof( gfp1 ) ) break;
			} else {
				if( fq1.eof() ) break;
			}
		}	//load and process file loop

		if( file_is_gz ) {
			gzclose( gfp1 );
			gzclose( gfp2 );
		} else {
			fq1.close();
			fq2.close();
		}
	}

	fout1.close();
	fout2.close();
	cerr << "\rDone: totally " << totalReads << " lines processed.\n";

	// write trim.log
	ofstream fout( "Msuite2.trim.log" );
	if( fout.fail() ) { 
		cerr << "Error: cannot write log file!\n";
		return 4;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=thread; ++i ) {
		dropped_all += dropped[i];
		real_all += real_adapter[i];
		tail_all += tail_adapter[i];
	}
	fout << "Total\t"		<< totalReads	<< '\n'
		 << "Dropped\t"		<< dropped_all	<< '\n'
		 << "Aadaptor\t"	<< real_all		<< '\n'
		 << "Tail Hit\t"	<< tail_all		<< '\n';
	fout.close();

	// write fqstatistics
	fout.open( "R1.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R1.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << AllR1stat[j].A << '\t' << AllR1stat[j].C << '\t'<< AllR1stat[j].G
				<< '\t' << AllR1stat[j].T << '\t'<< AllR1stat[j].N << '\n';
	}
	fout.close();

	fout.open( "R2.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R2.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << AllR2stat[j].A << '\t' << AllR2stat[j].C << '\t'<< AllR2stat[j].G
				<< '\t' << AllR2stat[j].T << '\t'<< AllR2stat[j].N << '\n';
	}
	fout.close();

	fout.open( "R1.trimmed.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R1.trimmed.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << AllR1stat_trimmed[j].A << '\t' << AllR1stat_trimmed[j].C << '\t'<< AllR1stat_trimmed[j].G
				<< '\t' << AllR1stat_trimmed[j].T << '\t'<< AllR1stat_trimmed[j].N << '\n';
	}
	fout.close();

	fout.open( "R2.trimmed.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R2.trimmed.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << AllR2stat_trimmed[j].A << '\t' << AllR2stat_trimmed[j].C << '\t'<< AllR2stat_trimmed[j].G
				<< '\t' << AllR2stat_trimmed[j].T << '\t'<< AllR2stat_trimmed[j].N << '\n';
	}
	fout.close();

	//free memory
	delete [] gz_buffer_1;
	delete [] gz_buffer_2;
	for(unsigned int i=0; i!=thread; ++i) {
		delete buffer1[i];
		delete buffer2[i];
		delete R1stat[i];
		delete R2stat[i];
	}
	delete [] buffer1;
	delete [] buffer2;
	delete [] R1stat;
	delete [] R2stat;
	delete [] AllR1stat;
	delete [] AllR2stat;
	delete [] AllR1stat_trimmed;
	delete [] AllR2stat_trimmed;

	return 0;
}

