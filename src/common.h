/*
 * common.h
 *
 * This header file records the constants used in Msuite2
 *
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date  : Jan 2022
 *
 *  Adapters
 *
 *  sequence model:
 *
 *                *read1 -->
 *  5' adapter - sequence sequence sequence - 3' adapter
 *                                <-- read2*
 *
 *  so read1 may contains 3' adapter, name it adapter_r1
 *  read2 may contains reversed and ACGT-paired 5' adapter, i.e., adapter_r2
 *
 *  Original adapter pair from illumina:
 *  5': ACACTCTTTCCCTACACGACGCTCTTCCGATCT
 *  3': AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
 *
 * */

#ifndef _MSUITE_COMMON_
#define _MSUITE_COMMON_

// fastq statistics structure
typedef struct {
	unsigned int A;
	unsigned int C;
	unsigned int G;
	unsigned int T;
	unsigned int N;
} fastqstat;

typedef struct {
	const char* adapter_r1;
	const char* adapter_r2;
	const char* adapter_index;
	const unsigned int adapter_len;
} adapter_info;

//illumina sequencing adapters
const char* illumina_adapter_sequence = "AGATCGGAAGAGC";
const unsigned int illumina_adapter_len = 13;	//strlen(illumina_adapter_sequence)
const char* illumina_adapter_index = "AGA";
const adapter_info illumina_adapter = {illumina_adapter_sequence, illumina_adapter_sequence,
										illumina_adapter_index, illumina_adapter_len};

//nextera sequencing adapters
const char* nextera_adapter_sequence = "CTGTCTCTTATACACATCT";
const unsigned int nextera_adapter_len = 19;	//strlen(nextera_adapter_sequence)
const char* nextera_adapter_index = "CTG";
const adapter_info nextera_adapter = {nextera_adapter_sequence, nextera_adapter_sequence,
										nextera_adapter_index, nextera_adapter_len};

//bgi sequencing adapters
const char* bgi_adapter1_sequence = "AAGTCGGAGGCCAAGCGGTC";
const char* bgi_adapter2_sequence = "AAGTCGGATCGTAGCCATGT";
const unsigned int bgi_adapter_len = 19;	//strlen(bgi_adapter_sequence)
const char* bgi_adapter_index = "AAG";
const adapter_info bgi_adapter = {bgi_adapter1_sequence, bgi_adapter2_sequence,
										bgi_adapter_index, bgi_adapter_len};

const char FILE_SEPARATOR = ',';		// separator if multiple files are provided

const int READS_PER_BATCH  = 1 << 20;	// process 1M reads per batch (for parallelization)
const int BUFFER_SIZE_PER_BATCH_READ = 1 << 30;	// 1G buffer for each thread to convert FASTQ

// seed and error configurations
const unsigned int impossible_seed = 10000;
const unsigned int MAX_READ_LENGTH = 1024;
const unsigned int MAX_CONVERSION  = 1024;

// C>T and G>A conversion related
const char NORMAL_SEQNAME_START  = '@';
//const char LINE_NUMBER_SEPARATOR = '+';
//const char CONVERSION_LOG_START  = '~';
const char CONVERSION_LOG_END    = '#';	// original version uses '&'
//const char CONVERSION_LOG_SEPARATOR = ':';	//original version uses ':'
const char CONVERSION_LOG_SEPARATOR = ';';	//':' is replaced because it's commonly used in illumina seqName
const char KEEP_QUAL_MARKER      = '|';
const char REPLACEMENT_CHAR_AT   = '~';		// bowtie2 will fail if the read is like "@@"; ~ is 126, which cannot be a Phred score

const int MAX_CHANGES_CNT  =  256;
const int MAX_SEQ_CYCLE    =  255;
const int MAX_SAMLINE_SIZE = 2048;
const int MAX_SEQNAME_SIZE = 1024;
const int MAX_ITERM_SIZE   =   32;	// MAX sam iterm (e.g., chr, mateinfo) size
const int MAX_CIGAR_SIZE   =  128;
// there could be something like chr12_GA_converted, so set to 32

// quality parameters for keeping the read
const int MIN_ALIGN_SCORE_KEEP = 5;

// quality parameters for methylation calling
const int MIN_ALIGN_SCORE_METH = 20;	// minimum alignment score to call methylation
const int MIN_BASEQUAL_SCORE   = 53;	// minimum phred score for a cycle to be used in meth-calling

#endif

