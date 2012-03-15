#include "configureReader.h"


configureReader::configureReader(void)
{
}

configureReader::~configureReader(void)
{
}

// Read the configuration file
// return: -1 failed; 1 succeed
int configureReader::ReadConfigFile( string fileName )
{
	ifstream configReader( fileName.c_str() );

	if( configReader == NULL )
	{
		cout<<"error reading configuration file"<<endl;
		return -1;
	}
		
	string line;					
			
	while( getline( configReader, line ) )
	{
		if( line.size() > 0 && line.at( 0 ) != '#' )
		{
			// analyze the parameters
			if( line == "[LIB]" )
				continue;

			if( AnalyzeParameters( line ) == -1 )
				return -1;
		}
	}

	return 1;
}

// anylize the parameters
// return: -1 if failed
int configureReader::AnalyzeParameters( string line ){
	// split the line according to "="
	int pos = line.find( "=" );
	if( pos == -1 )
		return pos;

	string propertyName = line.substr( 0, pos );
	string value = line.substr( pos + 1, line.length() - pos - 1 );
	if( propertyName == "file_format" ){
		if( value == "fasta" )
			Configure::FILE_FORMAT = FASTA;
		else if( value == "statistic" )
			Configure::FILE_FORMAT = STATISTIC;
		else
			return -1;
	}
	else if( propertyName == "file_type" ){
		if( value == "velvet" )
			Configure::FILE_TYPE = VELVET;
		else if( value == "soap" )
			Configure::FILE_TYPE = SOAP;
		else
			return -1;
	}
	else if( propertyName == "filter_repeat" ){
		if( value == "yes" )
			Configure::FILTER_REPEAT = true;
		else if( value == "no" )
			Configure::FILTER_REPEAT = false;
		else
			return -1;
	}
	else if( propertyName == "contig_file" ){
		Configure::CONTIG_FILE = value;	
	}
	else if( propertyName == "repeat_threshold" ){
		Configure::REPEAT_THRESHOLD = atof( value.c_str() );
	}
	else if( propertyName == "map_type" ){
		if( value == "bowtie" )
			Configure::MAP_TYPE = BOWTIE;
		else if( value == "opera" )
			Configure::MAP_TYPE = OPERA;
		else
			return -1;
	}
	else if( propertyName == "calculate_lib" ){
		if( value == "yes" )
			Configure::CALCULATE_LIB = true;
		else if( value == "no" )
			Configure::CALCULATE_LIB = false;
		else
			return -1;
	}
	else if( propertyName == "lib_mean" ){
		Configure::LIB_MEAN = atoi( value.c_str() );
		if( Configure::MULTI_LIB_MEAN == "" )
			Configure::MULTI_LIB_MEAN = value;
		else
			Configure::MULTI_LIB_MEAN += "," + value;
	}
	else if( propertyName == "lib_std" ){
		Configure::LIB_STD = atoi( value.c_str() );
		if( Configure::MULTI_LIB_STD == "" )
			Configure::MULTI_LIB_STD = value;
		else
			Configure::MULTI_LIB_STD += "," + value;
	}
	else if( propertyName == "calculate_ori" ){
		if( value == "yes" )
			Configure::CALCULATE_ORI = true;
		else if( value == "no" )
			Configure::CALCULATE_ORI = false;
		else
			return -1;
	}
	else if( propertyName == "read_ori" ){
		if( value == "in" )
			Configure::READ_ORI = IN;
		else if( value == "out" )
			Configure::READ_ORI = OUT;
		else if( value == "forward" )
			Configure::READ_ORI = FORWARD;
		else
			return -1;
	}
	else if( propertyName == "map_file" ){
		//Configure::MAP_FILE = value;
		if( Configure::MAP_FILE == "" )
			Configure::MAP_FILE = value;
		else
			Configure::MAP_FILE += "," + value;
	}
	else if( propertyName == "cluster_threshold" ){
		Configure::CLUSTER_THRESHOLD = atoi( value.c_str() );
	}
	else if( propertyName == "contig_size_threshold" ){
		Configure::CONTIG_SIZE_THERSHOLD = atoi( value.c_str() );
	}
	else if( propertyName == "output_folder" ){
		Configure::OUTPUT_FOLDER = value;
		mkdir (Configure::OUTPUT_FOLDER.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		/*if( _access( Configure::OUTPUT_FOLDER.c_str(), 0 ) == -1 ){
			// create folder
			_mkdir( Configure::OUTPUT_FOLDER.c_str() );
		}*/
		if( Configure::OUTPUT_FOLDER.substr( Configure::OUTPUT_FOLDER.length() - 1, 1 ) != Configure::DIRECTORY_SIGN )
			Configure::OUTPUT_FOLDER += Configure::DIRECTORY_SIGN;
	}
	else if( propertyName == "heuristic" ){
		if( value == "on" )
			Configure::HEURISTIC = true;
		else if( value == "off" )
			Configure::HEURISTIC = false;
		else
			return -1;
	}
	else if( propertyName == "step" ){
		Configure::STEP = atoi( value.c_str() );
	}
	else if( propertyName == "scaffold_name" ){
		Configure::SCAFFOLD_PREFEX = value + "_";
	}
	else if( propertyName == "abort" ){
		if( value == "true" )
			Configure::ABORT = true;
		else if( value == "false" )
			Configure::ABORT = false;
		else
			return -1;
	}
	else if( propertyName == "overlap_graph_file" ){
		Configure::OVERLAP_GRAPH_FILE = value;
	}
	else if( propertyName == "kmer" ){
		Configure::KMER = atoi( value.c_str() );
	}

	return 1;
}
