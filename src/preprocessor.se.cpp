#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include <zlib.h>
#include <chrono>
#include <thread>
#include "common.h"
#include "util.h"

using namespace std;

// hisat2 supports 256 character long of read id, and does not has --sam-no-qname-trunc option
const unsigned int MAX_CONVERTED_READ_ID = 256;
const static chrono::microseconds waiting_time_for_writing(100);
char * gz_buffer = NULL;
string unk;

/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite2 package, adapted from Ktrim
 * Date: Nov 2024
 *
 * The conversion log's rule:
 *	   @S|xx;xx$ => for reads with endC and C>T changes, 'S' is its quality score
 *	   @xx;xx$	 => for reads without endC
 *	   @$		 => for reads without endC and conversions
 *
**/

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * so use static values (e.g., 4) is not good
*/
bool check_mismatch_dynamic_SE( const string & s, unsigned int pos, const adapter_info* ai ) {
	register unsigned int mis=0;
	register unsigned int i, len;
	len = s.length() - pos;
	if( len > ai->adapter_len )
	  len = ai->adapter_len;
	register unsigned int max_mismatch_dynamic = len >> 3;
	if( (max_mismatch_dynamic<<3) != len )
	  ++ max_mismatch_dynamic;
	const char * p = s.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != ai->adapter_r1[i] ) {
			++ mis;
			if( mis > max_mismatch_dynamic )
			  return false;
		}
	}

	return true;
}

unsigned int load_text_file(ifstream & fq, string *ld_id, string *ld_seq, string *ld_qual) {
	unsigned int loaded = 0;
	while( true ) {
		getline( fq, ld_id  [ loaded ] );
		if( fq.eof() )break;
		getline( fq, ld_seq [ loaded ] );
		getline( fq, unk );
		getline( fq, ld_qual[ loaded ] );

		++ loaded;
		if( loaded == READS_PER_BATCH )
			break;
	}
	return loaded;
}

unsigned int load_gz_file(gzFile & gfp, string *ld_id, string *ld_seq, string *ld_qual) {
	unsigned int loaded = 0;
	while( true ) {
		if( gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE ) == NULL ) break;
		ld_id  [ loaded ] = gz_buffer;
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );
		ld_seq [ loaded ] = gz_buffer;
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );	// this line is useless
		gzgets( gfp, gz_buffer, MAX_SEQNAME_SIZE );
		ld_qual[ loaded ] = gz_buffer;

		ld_id [ loaded ].pop_back();	// trim the tail '\n'
		ld_seq[ loaded ].pop_back();
		if( ld_qual[loaded].size() != ld_seq[loaded].size() )	//the last read may not contain '\n'
			ld_qual[loaded].pop_back();

		++ loaded;
		if( loaded == READS_PER_BATCH )
			break;
	}

	return loaded;
}

int main( int argc, char *argv[] ) {
	if( argc < 5 ) {
		cerr << "\nUsage: " << argv[0] << " <r1.fq> <r2.fq=placeholder> <cycle> <out.prefix> "
			 << "[mode=0|3|4] [thread=4] [min.length=36] [min.quality=53] [libraryKit=illumina] [cut.head=0] [cut.tail=0]\n\n"

			 << "This program is part of Msuite and is designed to do fastq statistics, quality-trimming,\n"
			 << "adapter-trimming and C->T conversions for Single-End reads.\n\n"

			 << "Run modes:\n"
			 << "  0. Do not perform any conversions. Usable for normal DNA/RNA-seq alignment.\n"
			 << "  3. Do C>T for all C sites. Suitable for 3-letter BS-seq data.\n"
			 << "  4. Do C>T for C in CpG sites. Suitable for 4-letter TAPS data.\n\n"

			 << "Default parameters:\n"
			 << "  mode: 0\n"
			 << "  thread: 4 (must be >=2)\n"
			 << "  min.length: 36\n"
			 << "  min.quality: 53 (33+20 for phred33('!') scoring system)\n"
			 << "    Phred64 to Phred33 conversion is automatically ON if min.quality >= 74\n\n";

		return 2;
	}
	// for SE, the r2.fq will be ignored (could be set to /dev/null ro anything else)
	// but to provide an identical interface to preprocessor.pe, its position is kept

	unsigned int mode = 0;
	unsigned int thread = 4;
	unsigned int min_length = 36;
	unsigned char quality = 53;
	const char *libraryKit = "illumina";
	const adapter_info* ai;
	unsigned int cut_head = 0;
	unsigned int cut_tail = 0;

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
							cut_head = (unsigned char) atoi( argv[10] );
							if( argc > 11 )
								cut_tail = (unsigned char) atoi( argv[11] );
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
	if( thread < 2 ) {
		cerr << "Error: thread must be >=2!\n";
		return 103;
	}
	unsigned int real_wk_thread = thread - 1;
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
	} else if ( strcmp(libraryKit, "bgi")==0 || strcmp(libraryKit, "Nextera")==0 ) {
		ai = &bgi_adapter;
	} else {
		cerr << "Error: invalid library kit! Currently only supports illumina, nextera, and bgi!\n";
		return 103;
	}

	bool changePhred = false;
	if( quality >= 74 ) {   // consider it is Phred 64, otherwise it means Phred>=40 which is impossible
		changePhred = true;
	}

	string base = argv[4];
	base += ".R1.fq";
	ofstream fout( base.c_str() );
	if( fout.fail() ) {
		cerr << "Error: write file failed!\n";
		return 3;
	}

	string *id_buffer1   = new string [READS_PER_BATCH];
	string *seq_buffer1  = new string [READS_PER_BATCH];
	string *qual_buffer1 = new string [READS_PER_BATCH];

	string *id_buffer2   = new string [READS_PER_BATCH];
	string *seq_buffer2  = new string [READS_PER_BATCH];
	string *qual_buffer2 = new string [READS_PER_BATCH];
	string unk;
	string *id, *seq, *qual;	// these are for processing
	string *ld_id, *ld_seq, *ld_qual;	//these are for loading

	register unsigned int line = 1;
	int *dropped	  = new int [thread];
	int *real_adapter = new int [thread];
	int *tail_adapter = new int [thread];

	fastqstat * Allstat = new fastqstat[ cycle ];
	memset( Allstat, 0, cycle*sizeof(fastqstat) );
	fastqstat * Allstat_trimmed = new fastqstat[ cycle ];
	memset( Allstat_trimmed, 0, cycle*sizeof(fastqstat) );

	// buffer for storing the modified reads per thread
	char ** buffer = new char *[thread];
	int  * bstored = new int   [thread];
	fastqstat **Rstat = new fastqstat * [thread];
	fastqstat **Rstat_trimmed = new fastqstat * [thread];
	for(unsigned int i=0; i!=thread; ++i) {
		buffer[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];
		Rstat[i]  = new fastqstat [ cycle ];
		Rstat_trimmed[i]  = new fastqstat [ cycle ];
		dropped[i] = 0;
		real_adapter[i] = 0;
		tail_adapter[i] = 0;
	}

	cerr << "Loading files ...\n";
	vector<string> Rs;
	string fileName="";
	for(unsigned int i=0; argv[1][i]!='\0'; ++i) {
		if( argv[1][i] == FILE_SEPARATOR ) {
			Rs.push_back( fileName );
			fileName.clear();
		} else {
			fileName += argv[1][i];
		}
	}
	Rs.push_back( fileName );
	unsigned int totalFiles = Rs.size();
	cout << "INFO: " << totalFiles << " singled fastq files will be loaded.\n";

	ifstream fq;
	gzFile gfp;
	gz_buffer = new char [ MAX_SEQNAME_SIZE ];
	register unsigned int totalReads = 0;
	for( unsigned int fileCnt=0; fileCnt!=totalFiles; ++ fileCnt ) {
		bool file_is_gz = false;
		register unsigned int i_file = Rs[fileCnt].size() - 3;
		register const char * p_file = Rs[fileCnt].c_str();
		if( p_file[i_file]=='.' && p_file[i_file+1]=='g' && p_file[i_file+2]=='z' ) {	// .gz file
			file_is_gz = true;
			gfp = gzopen( p_file, "r" );
			if( gfp==NULL ) {
				cerr << "Error: open gzipped fastq file failed!\n";
				fout.close();
				return 11;
			}
		} else {	// plain text
			file_is_gz = false;
			fq.open( Rs[fileCnt].c_str() );
			if( fq.fail() ) {
				cerr << "Error: open fastq file failed!\n";
				fout.close();
				return 11;
			}
		}

		// load first batch of data
		unsigned int loaded;
		id = id_buffer1;  seq = seq_buffer1; qual = qual_buffer1;
		ld_id = id_buffer2; ld_seq = seq_buffer2; ld_qual = qual_buffer2;
		if( file_is_gz ) {
			loaded = load_gz_file( gfp, id, seq, qual );
		} else {
			loaded = load_text_file( fq, id, seq, qual );
		}
		if( loaded == 0 ) {
			cerr << "Error: No data loaded!\n";
			exit(1);
//		} else {
//			cerr << "Info: start, " << loaded << " lines loaded.\n";
		}

		// load and process reads, batch by batch
		while( loaded ) {
			unsigned int loaded_batch = 0;
			unsigned int write_thread = 0;

			// start parallalization
			omp_set_num_threads( thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				if( tn == real_wk_thread )	{	// the last thread is for loading data
					if( file_is_gz ) {
						loaded_batch = load_gz_file( gfp, ld_id, ld_seq, ld_qual );
					} else {
						loaded_batch = load_text_file( fq, ld_id, ld_seq, ld_qual );
					}
//					cerr << "Info: " << loaded_batch << " lines loaded.\n";
				} else {
					unsigned int tn = omp_get_thread_num();
					unsigned int start = loaded * tn / real_wk_thread;
					unsigned int end   = loaded * (tn+1) / real_wk_thread;

					// normalization
					bstored[tn] = 0;
					memset( Rstat[tn], 0, cycle*sizeof(fastqstat) );
					memset( Rstat_trimmed[tn], 0, cycle*sizeof(fastqstat) );

					string conversionLog;
					register int i, j;
					register unsigned int last_seed;
					vector<unsigned int> seed;
					vector<unsigned int> :: iterator it;
					const char *p;
					char *conversion = new char [MAX_CONVERSION];
					char numstr[10]; // enough to hold all numbers up to 99,999,999 plus ':'

					for( unsigned int ii=start; ii!=end; ++ii ) {
						//if the reads are longer than "cycle" paramater, only keep the head "cycle" ones
						if( seq[ii].size() > cycle ) {
							seq[ii].resize(  cycle );
							qual[ii].resize( cycle );
						}

						// fqstatistics
						p = seq[ii].c_str();
						j = seq[ii].size();
						for( i=0; i!=j; ++i ) {
							switch ( p[i] ) {
								case 'a':
								case 'A': Rstat[tn][i].A ++; break;
								case 'c':
								case 'C': Rstat[tn][i].C ++; break;
								case 'g':
								case 'G': Rstat[tn][i].G ++; break;
								case 't':
								case 'T': Rstat[tn][i].T ++; break;
								default : Rstat[tn][i].N ++; break;
							}
						}

						// quality control
						p = qual[ii].c_str();
						i = get_quality_trim_cycle_se( p, qual[ii].length(), min_length, quality );
						if( i < min_length ) { // not long enough
							++ dropped[ tn ];
							continue;
						}

						seq[ii].resize(  i );
						qual[ii].resize( i );

						// looking for seed target, 1 mismatch is allowed for these 2 seeds
						// which means seq1 and seq2 at least should take 1 perfect seed match
						seed.clear();
						for( i=0; (i=seq[ii].find(ai->adapter_index, i)) != string::npos; ++i )
							seed.push_back( i );

						last_seed = impossible_seed;	// a position which cannot be in seed
						for( it=seed.begin(); it!=seed.end(); ++it ) {
							if( *it != last_seed ) {
								if( check_mismatch_dynamic_SE(seq[ii], *it, ai) )
									break;
								last_seed = *it;
							}
						}
						if( it != seed.end() ) {	// adapter found
							++ real_adapter[tn];
							if( *it >= min_length )	{
								seq[ii].resize(  *it );
								qual[ii].resize( *it );
							} else {	// drop this read as its length is not enough
								++ dropped[tn];
								continue;
							}
						} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
							i = seq[ii].length() - 2;
							p = seq[ii].c_str();
							if( p[i]==ai->adapter_r1[0] && p[i+1]==ai->adapter_r1[1] ) {
								if( i < min_length ) {
									++ dropped[tn];
									continue;
								}
								seq[ii].resize(  i );
								qual[ii].resize( i );

								++ tail_adapter[tn];
	/* it is not good to check tail-1 due to high false-positive
							} else {	// tail 2 is not good, check tail 1
								++ i;
								if( p[i] == ai->adapter_r1[0] ) {
									if( i < min_length ) {
										++ dropped[tn];
										continue;
									}
									seq1[ii].resize(  i );
									qual1[ii].resize( i );

									++ tail_adapter[tn];
								}
	*/
							}
						}

						// cut head and tail
						if( cut_tail ) { 
							if( seq[ii].length() > cut_tail ) {
								seq[ii].resize(  seq[ii].length()  - cut_tail );
								qual[ii].resize( qual[ii].length() - cut_tail );
							} else {
								++ dropped[tn];
								continue;
							}
						}
						// cut head
						if( cut_head ) {
							if( seq[ii].length() > cut_head ) {
								seq[ii].erase(  0, cut_head );
								qual[ii].erase( 0, cut_head );
							} else {
								++ dropped[tn];
								continue;
							}
						}

						if( seq[ii].length() < min_length ) {
							++ dropped[tn];
							continue;
						}

						j = seq[ii].size();
						// convert Phred64 to Phred33 if necessary
						if( changePhred ) {
							for( i=0; i!=j; ++i ) {
								qual[ii][i] -= 31;
							}
						}

						p = seq[ii].c_str();
						for( i=0; i!=j; ++i ) {
							switch ( p[i] ) {
								case 'a':
								case 'A': Rstat_trimmed[tn][i].A ++; break;
								case 'c':
								case 'C': Rstat_trimmed[tn][i].C ++; break;
								case 'g':
								case 'G': Rstat_trimmed[tn][i].G ++; break;
								case 't':
								case 'T': Rstat_trimmed[tn][i].T ++; break;
								default : Rstat_trimmed[tn][i].N ++; break;
							}
						}

						//check if there is any white space in the IDs
						//if so, remove all the data after the whitespace
						j = id[ii].size();
						p = id[ii].c_str();
						for( i=1; i!=j; ++i ) {
							if( p[i]==' ' || p[i]=='\t' ) {	// white space
								id[ii].resize( i );
								break;
							}
						}

						// do C->T conversion
						if( mode == 3 ) {	// in this implementation, id1 and id2 are different!!!
							// in mode 3, there is NO endC and frontG issues
							id[ii][0] = CONVERSION_LOG_END;
							j = seq[ii].size();	// seq1 and seq2 are of the same size
							conversionLog = NORMAL_SEQNAME_START;
							for( i=0; i!=j; ++i ) {
								if( seq[ii][i] == 'C' ) {
									seq[ii][i] = 'T';
									sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
									conversionLog += numstr;
								}
							}
							if( conversionLog.back() == CONVERSION_LOG_SEPARATOR )
								conversionLog.pop_back();

							if( conversionLog.length() > MAX_CONVERTED_READ_ID ) {
								++ dropped[tn];
								continue;
							}

							bstored[tn] += sprintf( buffer[tn]+bstored[tn], "%s%s\n%s\n+\n%s\n",
												conversionLog.c_str(), id[ii].c_str(), seq[ii].c_str(), qual[ii].c_str() );
						} else if ( mode == 4 ) {
							// check seq1 for C>T conversion
							id[ii][0] = CONVERSION_LOG_END;
							conversionLog = NORMAL_SEQNAME_START;
							j = seq[ii].size()-1;
							if( seq[ii].back() == 'C' ) { //ther is a 'C' and the end, discard it (but record its Quality score);
								//otherwise it may introduce a mismatch in alignment
								if( qual[ii].back() == '@' ) {
									conversionLog += REPLACEMENT_CHAR_AT;
								} else {
									conversionLog += qual[ii].back();
								}
								conversionLog += KEEP_QUAL_MARKER;
								seq[ii].pop_back();
								qual[ii].pop_back();
							}
							for( i=0; i!=j; ++i ) {
								if( seq[ii][i]=='C' && seq[ii][i+1]=='G' ) {
									seq[ii][i] = 'T';
									sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
									conversionLog += numstr;
								}
							}
							if( conversionLog.back() == CONVERSION_LOG_SEPARATOR )
								conversionLog.pop_back();

							if( conversionLog.length() > MAX_CONVERTED_READ_ID ) {
								++ dropped[tn];
								continue;
							}

							bstored[tn] += sprintf( buffer[tn]+bstored[tn], "%s%s\n%s\n+\n%s\n",
										conversionLog.c_str(), id[ii].c_str(), seq[ii].c_str(), qual[ii].c_str() );
						} else {	// no need to do conversion
							bstored[tn] += sprintf( buffer[tn]+bstored[tn], "%s\n%s\n+\n%s\n",
										id[ii].c_str(), seq[ii].c_str(), qual[ii].c_str() );
						}
					}

					// wait for my turn to output and update statistics
					while( true ) {
						if( tn == write_thread ) {// output to stdout for pipe with aligners?
							fout << buffer[tn];
							++ write_thread;

							for( unsigned int j=0; j!=cycle; ++j ) {
								Allstat[j].A += Rstat[tn][j].A;
								Allstat[j].C += Rstat[tn][j].C;
								Allstat[j].G += Rstat[tn][j].G;
								Allstat[j].T += Rstat[tn][j].T;
								Allstat[j].N += Rstat[tn][j].N;

								Allstat_trimmed[j].A += Rstat_trimmed[tn][j].A;
								Allstat_trimmed[j].C += Rstat_trimmed[tn][j].C;
								Allstat_trimmed[j].G += Rstat_trimmed[tn][j].G;
								Allstat_trimmed[j].T += Rstat_trimmed[tn][j].T;
								Allstat_trimmed[j].N += Rstat_trimmed[tn][j].N;
							}
							break;
						} else {
							this_thread::sleep_for( waiting_time_for_writing );
						}
					}
				}	// parallel body for each batch


			}
			totalReads += loaded;
			cerr << '\r' << totalReads << " reads finished";

			loaded = loaded_batch;
			//swap id and ld_id
			string *tmp;
			tmp=id;   id=ld_id;     ld_id=tmp;
			tmp=seq;  seq=ld_seq;   ld_seq=tmp;
			tmp=qual; qual=ld_qual; ld_qual=tmp;
		} // loop for current file
		// close file
		if( file_is_gz ) {
			gzclose( gfp );
		} else {
			fq.close();
		}
	}//loop for all files
	cerr << "\rDone: " << totalReads << " lines processed.\n";
	fout.close();

	// write trim.log
	fout.open( "Msuite2.trim.log" );
	if( fout.fail() ) { 
		cerr << "Error: cannot write log file!\n";
		return 4;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=real_wk_thread; ++i ) {
		dropped_all += dropped[i];
		real_all += real_adapter[i];
		tail_all += tail_adapter[i];
	}
	fout << "Total\t"	 << totalReads	<< '\n'
		 << "Dropped : " << dropped_all << '\n'
		 << "Aadaptor: " << real_all	<< '\n'
		 << "Tail Hit: " << tail_all	<< '\n';
	fout.close();

	// write fqstatistics
	fout.open( "R1.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R1.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << Allstat[j].A << '\t' << Allstat[j].C << '\t'<< Allstat[j].G
				<< '\t' << Allstat[j].T << '\t'<< Allstat[j].N << '\n';
	}
	fout.close();

	fout.open( "R1.trimmed.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R1.trimmed.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << Allstat_trimmed[j].A << '\t' << Allstat_trimmed[j].C
			 << '\t' << Allstat_trimmed[j].G << '\t' << Allstat_trimmed[j].T
			 << '\t'<< Allstat_trimmed[j].N << '\n';
	}
	fout.close();

	//free memory
	for(unsigned int i=0; i!=thread; ++i) {
		delete buffer[i];
		delete Rstat[i];
	}
	delete [] buffer;
	delete [] Rstat;
	delete [] Allstat;
	delete [] Allstat_trimmed;

	delete [] id_buffer1;
	delete [] seq_buffer1;
	delete [] qual_buffer1;
	delete [] id_buffer2;
	delete [] seq_buffer2;
	delete [] qual_buffer2;

	if( gz_buffer != NULL )
		delete [] gz_buffer;

	return 0;
}
