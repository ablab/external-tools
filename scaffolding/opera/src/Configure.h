#pragma once

#ifndef DEBUG
//#define DEBUG
#endif

#ifndef TIME
#define TIME
#endif

#ifndef LOG
#define LOG
#endif

#ifndef SPLIT
//#define SPLIT
#endif

#ifndef CHECK
//#define CHECK
#endif

#ifndef ORIENTATION
//#define ORIENTATION
#endif

#include <string>
#include <list>
#include <vector>

enum contig_file_format {FASTA, STATISTIC};
enum contig_file_type {VELVET, SOAP, SOAP_CONTIG, NORMAL};
enum map_type {BOWTIE, OPERA};
enum read_ori {IN, OUT, FORWARD, ERROR};
enum common_orientation {PLUS = 0, MINUS = 1, NONE, PLUSMINUS};
enum contig_position {START, END};
enum extend_orientation {RIGHT, LEFT, BOTH, NEITHER};
enum start_type {BORDER, ALL};

using namespace std;

class Configure
{
public:
	Configure(void);
	virtual ~Configure(void);

	//Attributes
	static int FILE_FORMAT;				// fasta or statistic
	static int FILE_TYPE;				// velvet or soap
	static bool FILTER_REPEAT;			// if need to filter the repeat contigs
	static double REPEAT_THRESHOLD;		// repeat threshold, default = 1.5
	static string CONTIG_FILE;			// the contig file
	static int CONTIG_SIZE_THERSHOLD;	// the contig size threshold, default = 500

	static string MAP_FILE;				// the mapping file
	static int MAP_TYPE;				// the mapping type: bowtie
	static bool CALCULATE_LIB;			// if need to calculate the library information
	static int LIB_MEAN;				// the mean length of library, default = 10000
	static int LIB_STD;					// the standard deviation of library, default = 1000
	static bool CALCULATE_ORI;			// if need to calculate the orientation of paired end reads
	static int READ_ORI;				// the orientation of reads(in, out or forward)

	static int STD_TIMES;				// the times of stdandard deviation, default is 6
	static int CLUSTER_THRESHOLD;		// the paired end read cluster threshold
	static string OUTPUT_FOLDER;		// the result directory

	static string DIRECTORY_SIGN;		// the sign of directory, /(linux) or \(windows)

	static string SCAFFOLD_PREFEX;		// name of obtained scaffolds

	static bool ABORT;					// if abort due to long time running

	// heuristic parameters
	static bool HEURISTIC;				// record if use heuristic method to speed up
	static int STEP;					// the max number of step to find unassigned nodes

	// threshold to consider an edge as invalid
	static int MIN_DIS;

	// overlap graph related variables
	static string OVERLAP_GRAPH_FILE;	// the file name of the overlap graph
	static int KMER;					// the size of the kmer

	// multiple libraries information
	static string MULTI_LIB_MEAN;		// the mean of libraries
	static string MULTI_LIB_STD;		// the std of libraries
};
