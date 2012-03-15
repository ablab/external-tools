#include "Configure.h"

Configure::Configure(void)
{
}

Configure::~Configure(void)
{
}

// contig file related parameters
int Configure::FILE_FORMAT = FASTA;
int Configure::FILE_TYPE = VELVET;
bool Configure::FILTER_REPEAT = true;
double Configure::REPEAT_THRESHOLD = 1.5;
string Configure::CONTIG_FILE = "";
int Configure::CONTIG_SIZE_THERSHOLD = 500;

// mapping file related parameters
string Configure::MAP_FILE = "";
int Configure::MAP_TYPE = BOWTIE;
bool Configure::CALCULATE_LIB = true;
int Configure::LIB_MEAN = 10000;
int Configure::LIB_STD = 1000;
bool Configure::CALCULATE_ORI = true;
int Configure::READ_ORI = IN;

// bundle related parameters
int Configure::STD_TIMES = 6;
int Configure::CLUSTER_THRESHOLD = 5;

// result directory
string Configure::OUTPUT_FOLDER = "";
string Configure::DIRECTORY_SIGN = "/";  //"\\";

string Configure::SCAFFOLD_PREFEX = "opera_scaffold_";

bool Configure::ABORT = true;

// heuristic
bool Configure::HEURISTIC = false;
int Configure::STEP = 5;

// minimum valid distance for an edge
int Configure::MIN_DIS = 0;

// overlap graph realated variables
string Configure::OVERLAP_GRAPH_FILE = "";
int Configure::KMER = 39;

// multiple libraries variables
string Configure::MULTI_LIB_MEAN = "";		// the mean of libraries
string Configure::MULTI_LIB_STD = "";		// the std of libraries
