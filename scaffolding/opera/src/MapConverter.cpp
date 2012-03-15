#include "MapConverter.h"

MapConverter::MapConverter(void)
{
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, less_distance>*>;
	m_pairMap = new list<string>;
	m_totalTime = 0;

	m_libString = "";
}

MapConverter::~MapConverter(void)
{
	delete m_singlePetsMap;
	delete m_pairMap;
}


MapConverter::MapConverter( Graph *graph ){
	m_graph = graph;
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, less_distance>*>;
	m_pairMap = new list<string>;
	m_totalTime = 0;
	m_pairTime = 0;
}

// analyze mapping file or opera edge file
int MapConverter::Analyze( string fileName ){
	switch( Configure::MAP_TYPE ){
		case BOWTIE:
			return AnalyzeBowtie( fileName );
		case OPERA:
			return AnalyzeOpera( fileName );
		default:
			return -1;
	}
	return 1;
}

// analyze mapping file or opera edge file
int MapConverter::AnalyzeMultiLib( string fileName, list<PetLibrary*> *libs  ){
	switch( Configure::MAP_TYPE ){
		case BOWTIE:
			return AnalyzeBowtieMultiLib( fileName, libs );
		default:
			return -1;
	}
	return 1;
}

// analyze bowtie file
// return -1 if failed
int MapConverter::AnalyzeBowtie( string fileName ){
	m_numOfPet = 0;

	// calculate lib mean/std and read orientation if needed
	if( Configure::CALCULATE_LIB || Configure::CALCULATE_ORI ){
		m_numOfPet = CalculateLibDis( fileName );
		if( m_numOfPet == -1 )
			return -1;

		if( ConvertBowtieFile() == -1 )
			return -1;
	}
	else if( ConvertBowtieFile( fileName ) == -1 ){
		// convert each paired-end reads into opera's format and output
		return -1;
	}

	return 1;
}

bool LibSort( PetLibrary *&lib1, PetLibrary *&lib2 ){
	return lib1->GetMean() < lib2->GetMean();
}

// analyze bowtie file
int MapConverter::AnalyzeBowtieMultiLib( string fileName, list<PetLibrary*> *libs ){
	// get all mapping files and handle them one by one
	vector<string> *files = new vector<string>;
	Split( fileName, ",", files );
	for( int i = 0; i < files->size(); i++ )
	{
		// first calculate the mean and std for each library
		//cout<<files->at( i )<<endl;
		if( CalculateLibDisMultiLib( files->at( i ), libs ) == -1 ){
			return -1;
		}

		if( ConvertBowtieFileMultiLib( files->at( i ), libs ) == -1 )
			return -1;
	}

	// output library information
	ofstream libWriter( (Configure::OUTPUT_FOLDER + "lib.txt").c_str() );
	if( libWriter == NULL ){
		cout<<"error writing lib.txt file!"<<endl;
		return -1;
	}
	libWriter.write( m_libString.c_str(), m_libString.length() );
	libWriter.close();
	
	// sort according to mean
	libs->sort( LibSort );

	files->clear();
	delete files;
}

// calculate library mean and std
// return: the number of paired reads mapping on different contigs
// return -1 if failed
int MapConverter::CalculateLibDis( string fileName ){
	ifstream mapReader( fileName.c_str() );

	int oriNum[ 4 ];
	for( int i = 0; i < 4; i++ )
		oriNum[ i ] = 0;

	if( mapReader == NULL )
	{
		cout<<"error reading mapping file: "<<fileName<<endl;
		return -1;
	}

	int numOfPet = 0;
	
	string preLine, nextLine;
	getline( mapReader, preLine );

	double sum = 0;
	double num = 0;
	vector<double> distance;

	vector<string> *preColumn = new vector<string>;
	vector<string> *nextColumn = new vector<string>;
	vector<string> *preNames = new vector<string>;
	vector<string> *nextNames = new vector<string>;

	char preName[ 200 ], nextName[ 200 ];
	while( getline( mapReader, nextLine ) ){
		sscanf( preLine.c_str(), "%s", preName );
		sscanf( nextLine.c_str(), "%s", nextName );
		string preNameString = string( preName );
		string nextNameString = string( nextName );
		

		if( IsPair( preNameString, nextNameString ) ){
			// is a pair
			Split( preLine, "\t", preColumn );
			Split( nextLine, "\t", nextColumn );

			// for the name, only select the word before the first space
			(*preColumn)[ 0 ] = preNameString;
			(*nextColumn)[ 0 ] = nextNameString;

			if( (*preColumn)[ 2 ] == (*nextColumn)[ 2 ] ){
				// mapped to the same contig
				double dis = CalculateReadDisOnSameContig( preColumn, nextColumn );
				distance.push_back( dis );
				sum += dis;
				num++;

				// check orientation
				string firstOri = (*preColumn)[ 1 ];
				string secondOri = (*nextColumn)[ 1 ];
				int firstPos = atoi( (*preColumn)[ 3 ].c_str() );
				int secondPos = atoi( (*nextColumn)[ 3 ].c_str() );
				
				if( ( firstOri == "+" && secondOri == "+" && firstPos > secondPos ) || ( firstOri == "-" && secondOri == "-" && firstPos < secondPos )  )
					oriNum[ FORWARD ]++;
				else if( ( firstOri == "+" && secondOri == "-" && firstPos < secondPos ) || ( firstOri == "-" && secondOri == "+" && firstPos > secondPos )  )
					oriNum[ IN ]++;
				else if( ( firstOri == "-" && secondOri == "+" && firstPos < secondPos ) || ( firstOri == "+" && secondOri == "-" && firstPos > secondPos )  )
					oriNum[ OUT ]++;
				else
					oriNum[ ERROR ]++;
			}
			else{
				if( m_graph->GetContigIndex( (*preColumn)[ 2 ] ) != -1 
				&& m_graph->GetContigIndex( (*nextColumn)[ 2 ] ) != -1 ){ 
					m_pairMap->push_back( preLine );
					m_pairMap->push_back( nextLine );
					numOfPet++;
				}
			}

			getline( mapReader, preLine );
		}
		else{
			// not a pair
			preLine = nextLine;
		}
	}

	double mean = sum / num;

	double std = 0;
	// calculate std
	for( int i = 0; i < distance.size(); i++ )
		std += pow( distance.at( i ) - mean, 2 );

	std = std / num;
	std = pow( std, 0.5 );

	if( Configure::CALCULATE_LIB ){
		Configure::LIB_MEAN = (int) mean;
		Configure::LIB_STD = (int) std;
	}

	if( Configure::CALCULATE_ORI ){
		int max = 0;
		int maxID = -1;
		for( int i = 0; i < 4; i++ ){
			if( oriNum[ i ] > max ){
				max = oriNum[ i ];
				maxID = i;
			}
		}

		if( maxID == ERROR ){
			cout<<"The orientation of reads is not \"in\", \"out\", neither \"forward\", Opera could not handle such orientation. "
				<<"Please check the mapping files. Please feel free to contact gaosong@nus.edu.sg for further support. \n";
			return -1;
		}

#ifdef ORIENTATION
		if( maxID == IN )
			cout<<"orientation is in\n";
		else if( maxID == OUT )
			cout<<"orientation is out\n";
		else if( maxID == FORWARD )
			cout<<"orientation is forward\n";
#endif
		Configure::READ_ORI = maxID;
	}

	ofstream libWriter( (Configure::OUTPUT_FOLDER + "lib.txt").c_str() );
	if( libWriter == NULL ){
		cout<<"error writing lib.txt file!"<<endl;
		return -1;
	}

	string libString = "Mean length of library is: " + itos( (int)mean ) + "\n";
	libString.append( "Standard deviation of library is: " + itos( (int) std ) + "\n" );
	if( Configure::CALCULATE_ORI ){
		libString.append( "Orientation of paired-reads is: " );
		if( Configure::READ_ORI == IN )
			libString.append( "in ->...<-\n" );
		else if( Configure::READ_ORI == OUT )
			libString.append( "out <-...->\n" );
		else if( Configure::READ_ORI == FORWARD )
			libString.append( "forward ->...->  2nd read...1st read\n" );
	}
	libWriter.write( libString.c_str(), libString.length() );
	libWriter.close();

	mapReader.close();

	preColumn->clear();
	delete preColumn;
	nextColumn->clear();
	delete nextColumn;
	preNames->clear();
	delete preNames;
	nextNames->clear();
	delete nextNames;
	
	return numOfPet;
}

// calculate library mean and std
int MapConverter::CalculateLibDisMultiLib( string fileName, list<PetLibrary*> *libs ){
	PetLibrary *newLib = new PetLibrary( fileName );
	int libID = libs->size();
	libs->push_back( newLib );
	if( Configure::MULTI_LIB_MEAN != "" )
	{
		vector<string> *mean = new vector<string>;
		vector<string> *std = new vector<string>;
		Split( Configure::MULTI_LIB_MEAN, ",", mean );
		Split( Configure::MULTI_LIB_STD, ",", std );
		newLib->SetMean( atoi( (*mean)[ libID ].c_str() ) );
		newLib->SetStd( atoi( (*std)[ libID ].c_str() ) );
		Configure::CALCULATE_LIB = false;

		mean->clear();
		std->clear();
		delete mean;
		delete std;
	}
	else
		Configure::CALCULATE_LIB = true;

	ifstream mapReader( fileName.c_str() );

	int oriNum[ 4 ];
	for( int i = 0; i < 4; i++ )
		oriNum[ i ] = 0;

	if( mapReader == NULL )
	{
		cout<<"error reading mapping file: "<<fileName<<endl;
		return -1;
	}

	int numOfPet = 0;
	
	string preLine, nextLine;
	getline( mapReader, preLine );

	double sum = 0;
	double num = 0;
	vector<double> distance;

	vector<string> *preColumn = new vector<string>;
	vector<string> *nextColumn = new vector<string>;
	vector<string> *preNames = new vector<string>;
	vector<string> *nextNames = new vector<string>;

	char preName[ 200 ], nextName[ 200 ];
	while( getline( mapReader, nextLine ) ){
		sscanf( preLine.c_str(), "%s", preName );
		sscanf( nextLine.c_str(), "%s", nextName );
		string preNameString = string( preName );
		string nextNameString = string( nextName );
		

		if( IsPair( preNameString, nextNameString ) ){
			// is a pair
			Split( preLine, "\t", preColumn );
			Split( nextLine, "\t", nextColumn );

			// for the name, only select the word before the first space
			(*preColumn)[ 0 ] = preNameString;
			(*nextColumn)[ 0 ] = nextNameString;

			if( (*preColumn)[ 2 ] == (*nextColumn)[ 2 ] ){
				// mapped to the same contig
				double dis = CalculateReadDisOnSameContig( preColumn, nextColumn );
				distance.push_back( dis );
				sum += dis;
				num++;

				// check orientation
				string firstOri = (*preColumn)[ 1 ];
				string secondOri = (*nextColumn)[ 1 ];
				int firstPos = atoi( (*preColumn)[ 3 ].c_str() );
				int secondPos = atoi( (*nextColumn)[ 3 ].c_str() );
				
				if( ( firstOri == "+" && secondOri == "+" && firstPos > secondPos ) || ( firstOri == "-" && secondOri == "-" && firstPos < secondPos )  )
					oriNum[ FORWARD ]++;
				else if( ( firstOri == "+" && secondOri == "-" && firstPos < secondPos ) || ( firstOri == "-" && secondOri == "+" && firstPos > secondPos )  )
					oriNum[ IN ]++;
				else if( ( firstOri == "-" && secondOri == "+" && firstPos < secondPos ) || ( firstOri == "+" && secondOri == "-" && firstPos > secondPos )  )
					oriNum[ OUT ]++;
				else
					oriNum[ ERROR ]++;
			}
			else{
				if( m_graph->GetContigIndex( (*preColumn)[ 2 ] ) != -1 
				&& m_graph->GetContigIndex( (*nextColumn)[ 2 ] ) != -1 ){ 
					m_pairMap->push_back( preLine );
					m_pairMap->push_back( nextLine );
					numOfPet++;
				}
			}

			getline( mapReader, preLine );
		}
		else{
			// not a pair
			preLine = nextLine;
		}
	}

	double mean = sum / num;

	double std = 0;
	// calculate std
	for( int i = 0; i < distance.size(); i++ )
		std += pow( distance.at( i ) - mean, 2 );

	std = std / num;
	std = pow( std, 0.5 );

	if( Configure::CALCULATE_LIB ){
		newLib->SetMean( (int) mean );
		newLib->SetStd( (int) std );
		//cout<<"mean: "<<mean<<endl;
		//cout<<"std: "<<std<<endl;
	}

	if( Configure::CALCULATE_ORI ){
		int max = 0;
		int maxID = -1;
		for( int i = 0; i < 4; i++ ){
			if( oriNum[ i ] > max ){
				max = oriNum[ i ];
				maxID = i;
			}
		}

		if( maxID == ERROR ){
			cout<<"The orientation of reads is not \"in\", \"out\", neither \"forward\", Opera could not handle such orientation. "
				<<"Please check the mapping files. Please feel free to contact gaosong@nus.edu.sg for further support. \n";
			return -1;
		}

#ifdef ORIENTATION
		if( maxID == IN )
			cout<<"orientation is in\n";
		else if( maxID == OUT )
			cout<<"orientation is out\n";
		else if( maxID == FORWARD )
			cout<<"orientation is forward\n";
#endif
		newLib->SetOri( maxID );
	}

	// record the info of this library
	m_libString += "For mapping file: " + fileName + "\n";
	m_libString += "Mean length of library is: " + itos( (int)newLib->GetMean() ) + "\n";
	m_libString.append( "Standard deviation of library is: " + itos( (int)newLib->GetStd() ) + "\n" );
	if( Configure::CALCULATE_ORI ){
		m_libString.append( "Orientation of paired-reads is: " );
		if( Configure::READ_ORI == IN )
			m_libString.append( "in ->...<-\n" );
		else if( Configure::READ_ORI == OUT )
			m_libString.append( "out <-...->\n" );
		else if( Configure::READ_ORI == FORWARD )
			m_libString.append( "forward ->...->  2nd read...1st read\n" );
	}
	m_libString.append( "\n" );

	mapReader.close();

	preColumn->clear();
	delete preColumn;
	nextColumn->clear();
	delete nextColumn;
	preNames->clear();
	delete preNames;
	nextNames->clear();
	delete nextNames;
	
	return numOfPet;
}

// check if two lines represent a pair of reads
bool MapConverter::IsPair( string firstRead, string secondRead ){
	string firstName = firstRead.substr( 0, firstRead.length() -2 );
	string secondName = secondRead.substr( 0, secondRead.length() - 2 );

	if( firstName == secondName ){
		return true;
	}
	else
		return false;
}

// calculate the distance of two reads in the same contig
// including the read length
int MapConverter::CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn ){
	int number[ 4 ];
	number[ 0 ] = atoi( firstColumn->at( 3 ).c_str() );
	number[ 1 ] = number[ 0 ] + firstColumn->at( 4 ).length();
	number[ 2 ] = atoi( secondColumn->at( 3 ).c_str() );
	number[ 3 ] = number[ 2 ] + secondColumn->at( 4 ).length();

	int min = number[ 0 ];
	int max = min;

	for( int i = 1; i < 4; i++ ){
		if( number[ i ] < min )
			min = number[ i ];
		if( number[ i ] > max )
			max = number[ i ];
	}

	return max - min;
}

// convert bowtie format into opera format
// return -1 if failed
int MapConverter::ConvertBowtieFile( string fileName ){
	int sumTime = 0;
	ifstream mapReader( fileName.c_str() );

	if( mapReader == NULL )
	{
		cout<<"error reading mapping file: "<<fileName<<endl;
		return -1;
	}

	ofstream mapWriter( (Configure::OUTPUT_FOLDER + "pairedEndReads").c_str(), ios::out );
	if( mapWriter == NULL ){
		cout<<"error opening pairedEndReads file!"<<endl;
		return -1;
	}
	
	string head = "ID\tFirst Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\n";
	mapWriter.write( head.c_str(), head.length() );
	string preLine, nextLine;
	getline( mapReader, preLine );

	int petID = 0;			// record the paired-end reads ID

	vector<string> *preColumn = new vector<string>;
	vector<string> *nextColumn = new vector<string>;
	vector<string> *preNames = new vector<string>;
	vector<string> *nextNames = new vector<string>;

	char preName[ 200 ], nextName[ 200 ];
	while( getline( mapReader, nextLine ) ){
		sscanf( preLine.c_str(), "%s", preName );
		sscanf( nextLine.c_str(), "%s", nextName );
		string preNameString = string( preName );
		string nextNameString = string( nextName );

		if( IsPair( preNameString, nextNameString ) ){
			// is a pair
			Split( preLine, "\t", preColumn );
			Split( nextLine, "\t", nextColumn );

			// for the name, only select the word before the first space
			(*preColumn)[ 0 ] = preNameString;
			(*nextColumn)[ 0 ] = nextNameString;

			if( (*preColumn)[ 2 ] != (*nextColumn)[ 2 ] ){
				// mapped to the different contig
				// check if both contigs exist
				int pos1 = m_graph->GetContigIndex( (*preColumn)[ 2 ] );
				int pos2 = m_graph->GetContigIndex( (*nextColumn)[ 2 ] );
				if( pos1 == -1 || pos2 == -1 ){
					// one of the contigs does not exist
					getline( mapReader, preLine );
					continue;
				}
				
				string newLine = ConvertBowtieFormat( preColumn, nextColumn, petID++ );
				mapWriter.write( newLine.c_str(), newLine.size() );
			}

			getline( mapReader, preLine );
			
		}
		else{
			// not a pair
			preLine = nextLine;
		}
	}

	delete preColumn;
	delete nextColumn;

	mapWriter.close();
	mapReader.close();

	return 1;
}

// convert bowtie format into opera format
int MapConverter::ConvertBowtieFileMultiLib( string fileName, list<PetLibrary*> *libs ){
	int sumTime = 0;
	
	vector<string> *path = new vector<string>;
	Split( fileName, "/", path );
	vector<string> *names = new vector<string>;
	Split( path->at( path->size() - 1 ), ".", names );
	ofstream mapWriter( (Configure::OUTPUT_FOLDER + "pairedEndReads_" + names->at( 0 )).c_str(), ios::out );
	libs->back()->SetNameWithoutPath( names->at( 0 ) );
	if( mapWriter == NULL ){
		cout<<"error opening pairedEndReads("<<names->at( 0 )<<") file!"<<endl;
		names->clear();
		delete names;
		path->clear();
		delete path;
		return -1;
	}
	names->clear();
	delete names;
	path->clear();
	delete path;
	
	string head = "ID\tFirst Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\n";
	mapWriter.write( head.c_str(), head.length() );
	
	// no pairing information
	if( m_pairMap->empty() )
		return 0;

	string preLine, nextLine;
	preLine = m_pairMap->front();
	m_pairMap->pop_front();

	int petID = 0;			// record the paired-end reads ID

	vector<string> *preColumn = new vector<string>;
	vector<string> *nextColumn = new vector<string>;
	vector<string> *preNames = new vector<string>;
	vector<string> *nextNames = new vector<string>;

	char preName[ 200 ], nextName[ 200 ];
	while( !m_pairMap->empty() ){
		nextLine = m_pairMap->front();
		m_pairMap->pop_front();

		sscanf( preLine.c_str(), "%s", preName );
		sscanf( nextLine.c_str(), "%s", nextName );
		string preNameString = string( preName );
		string nextNameString = string( nextName );

		if( IsPair( preNameString, nextNameString ) ){
			// is a pair
			Split( preLine, "\t", preColumn );
			Split( nextLine, "\t", nextColumn );

			// for the name, only select the word before the first space
			(*preColumn)[ 0 ] = preNameString;
			(*nextColumn)[ 0 ] = nextNameString;

			if( (*preColumn)[ 2 ] != (*nextColumn)[ 2 ] ){
				// mapped to the different contig
				// check if both contigs exist
				int pos1 = m_graph->GetContigIndex( (*preColumn)[ 2 ] );
				int pos2 = m_graph->GetContigIndex( (*nextColumn)[ 2 ] );
				if( pos1 == -1 || pos2 == -1 ){
					// one of the contigs does not exist
					preLine = m_pairMap->front();
					m_pairMap->pop_front();
					continue;
				}
				
				string newLine = ConvertBowtieFormatMultiLib( preColumn, nextColumn, petID++, libs );
				mapWriter.write( newLine.c_str(), newLine.size() );
			}

			if( !m_pairMap->empty() ){
				preLine = m_pairMap->front();
				m_pairMap->pop_front();
			}
			else{
				break;
			}
			
		}
		else{
			// not a pair
			preLine = nextLine;
		}
	}

	delete preColumn;
	delete nextColumn;

	mapWriter.close();

	return 1;
}

// convert bowtie format into opera format
// all the paired reads are already saved to m_pairMap. 
// no need to read file again
// return -1 if failed
int MapConverter::ConvertBowtieFile(){
	int sumTime = 0;
	
	ofstream mapWriter( (Configure::OUTPUT_FOLDER + "pairedEndReads").c_str(), ios::out );
	if( mapWriter == NULL ){
		cout<<"error opening pairedEndReads file!"<<endl;
		return -1;
	}
	
	string head = "ID\tFirst Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\n";
	mapWriter.write( head.c_str(), head.length() );
	
	// no pairing information
	if( m_pairMap->empty() )
		return 0;

	string preLine, nextLine;
	preLine = m_pairMap->front();
	m_pairMap->pop_front();

	int petID = 0;			// record the paired-end reads ID

	vector<string> *preColumn = new vector<string>;
	vector<string> *nextColumn = new vector<string>;
	vector<string> *preNames = new vector<string>;
	vector<string> *nextNames = new vector<string>;

	char preName[ 200 ], nextName[ 200 ];
	while( !m_pairMap->empty() ){
		nextLine = m_pairMap->front();
		m_pairMap->pop_front();

		sscanf( preLine.c_str(), "%s", preName );
		sscanf( nextLine.c_str(), "%s", nextName );
		string preNameString = string( preName );
		string nextNameString = string( nextName );

		if( IsPair( preNameString, nextNameString ) ){
			// is a pair
			Split( preLine, "\t", preColumn );
			Split( nextLine, "\t", nextColumn );

			// for the name, only select the word before the first space
			(*preColumn)[ 0 ] = preNameString;
			(*nextColumn)[ 0 ] = nextNameString;

			if( (*preColumn)[ 2 ] != (*nextColumn)[ 2 ] ){
				// mapped to the different contig
				// check if both contigs exist
				int pos1 = m_graph->GetContigIndex( (*preColumn)[ 2 ] );
				int pos2 = m_graph->GetContigIndex( (*nextColumn)[ 2 ] );
				if( pos1 == -1 || pos2 == -1 ){
					// one of the contigs does not exist
					preLine = m_pairMap->front();
					m_pairMap->pop_front();
					continue;
				}
				
				string newLine = ConvertBowtieFormat( preColumn, nextColumn, petID++ );
				mapWriter.write( newLine.c_str(), newLine.size() );
			}

			if( !m_pairMap->empty() ){
				preLine = m_pairMap->front();
				m_pairMap->pop_front();
			}
			else{
				break;
			}
			
		}
		else{
			// not a pair
			preLine = nextLine;
		}
	}

	delete preColumn;
	delete nextColumn;

	mapWriter.close();

	return 1;
}



// convert bowtie format
string MapConverter::ConvertBowtieFormat( vector<string> *firstColumn, vector<string> *secondColumn, int id ){
	string result = itos( id ) + "\t" + firstColumn->at( 2 ) + "\t";

	int pos1 = m_graph->GetContigIndex( firstColumn->at( 2 ) );
	int pos2 = m_graph->GetContigIndex( secondColumn->at( 2 ) );
	Contig *contig1 = m_graph->GetContig( pos1 );
	Contig *contig2 = m_graph->GetContig( pos2 );
	double leftDis, rightDis;
	int dis;
	string firstOri, secondOri;

	if( Configure::READ_ORI == IN ){
		// orientation is ->...<- (+,-)
		if( firstColumn->at( 1 ) == "+" ){
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
			else{
				// read ori is + -; contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
		}
		else{
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
			else{
				// read ori is - -; contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
		}
	} // end dealing with orientation "in"
	else if( Configure::READ_ORI == OUT ){
		// orientation is <-...-> (-,+)
		if( firstColumn->at( 1 ) == "+" ){
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is + -; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
		else{
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is - -; contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
	}
	else if( Configure::READ_ORI == FORWARD ){
		// orientation is ->...-> (+,+), 2nd read...1st read
		// swap read;
		vector<string> *temp = secondColumn;
		secondColumn = firstColumn;
		firstColumn = temp;

		int tempPos = pos1;
		pos1 = pos2;
		pos2 = tempPos;
		Contig *tempContig = contig1;
		contig1 = contig2;
		contig2 = tempContig;

		/*bool iffind = false;
		if( ( firstColumn->at( 2 ) == "scaffold21450" && secondColumn->at( 2 ) == "scaffold270237" )
			|| ( secondColumn->at( 2 ) == "scaffold21450" && firstColumn->at( 2 ) == "scaffold270237" ) ){
				cout<<firstColumn->at( 0 )<<"\t"<<firstColumn->at( 1 )<<"\t"<<firstColumn->at( 2 )<<"\t"<<firstColumn->at( 3 )<<endl;
				cout<<secondColumn->at( 0 )<<"\t"<<secondColumn->at( 1 )<<"\t"<<secondColumn->at( 2 )<<"\t"<<secondColumn->at( 3 )<<endl;
				iffind = true;
		}*/


		if( firstColumn->at( 1 ) == "+" ){
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is + -; contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
		else{
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();

			}
			else{
				// read ori is - -; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
	}

	dis = (int) (Configure::LIB_MEAN - leftDis - rightDis);
	result += itos( dis ) + "\t" + itos( Configure::LIB_STD ) + "\n";
	if( pos1 > pos2 ){
		// change two contigs
		int tempPos = pos1;
		pos1 = pos2;
		pos2 = tempPos;
		string tempOri = firstOri;
		if( secondOri == "+" )
			firstOri = "-";
		else
			firstOri = "+";
		if( tempOri == "+" )
			secondOri = "-";
		else
			secondOri = "+";

	}
	// create a single pet object
	SinglePet *newPet = new SinglePet( id, pos1, firstOri, pos2, secondOri, dis, Configure::LIB_STD );
	//m_singlePets->push_back( newPet );

	// save to map
	int firstID = pos1;
	int secondID = pos2;
	if( firstOri == "-" )
		firstID = -firstID;
	if( secondOri == "-" )
		secondID = -secondID;

	pair<int, int> newPair( firstID, secondID );
	map<pair<int, int>, multiset<SinglePet*, less_distance>*>::iterator mapIter = m_singlePetsMap->find( newPair );
	if( mapIter == m_singlePetsMap->end() ){
		// new pair, create a new set
		multiset<SinglePet*, less_distance> *newSet = new multiset<SinglePet*, less_distance>;
		newSet->insert( newPet );
		m_singlePetsMap->insert( pair<pair<int, int>, multiset<SinglePet*, less_distance>*>( newPair, newSet ) );
	}
	else{
		// insert the single pet to previous set
		mapIter->second->insert( newPet );
	}

	return result;
}

// convert bowtie format
string MapConverter::ConvertBowtieFormatMultiLib( vector<string> *firstColumn, vector<string> *secondColumn, int id, list<PetLibrary*> *libs ){
	// get the current library
	PetLibrary *currentLib = libs->back();
	string result = itos( id ) + "\t" + firstColumn->at( 2 ) + "\t";

	int pos1 = m_graph->GetContigIndex( firstColumn->at( 2 ) );
	int pos2 = m_graph->GetContigIndex( secondColumn->at( 2 ) );
	Contig *contig1 = m_graph->GetContig( pos1 );
	Contig *contig2 = m_graph->GetContig( pos2 );
	double leftDis, rightDis;
	int dis;
	string firstOri, secondOri;

	if( currentLib->GetOri() == IN ){
		// orientation is ->...<- (+,-)
		if( firstColumn->at( 1 ) == "+" ){
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
			else{
				// read ori is + -; contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
		}
		else{
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
			else{
				// read ori is - -; contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
		}
	} // end dealing with orientation "in"
	else if( currentLib->GetOri() == OUT ){
		// orientation is <-...-> (-,+)
		if( firstColumn->at( 1 ) == "+" ){
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is + -; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
		else{
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is - -; contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
	}
	else if( currentLib->GetOri() == FORWARD ){
		// orientation is ->...-> (+,+), 2nd read...1st read
		// swap read;
		vector<string> *temp = secondColumn;
		secondColumn = firstColumn;
		firstColumn = temp;

		int tempPos = pos1;
		pos1 = pos2;
		pos2 = tempPos;
		Contig *tempContig = contig1;
		contig1 = contig2;
		contig2 = tempContig;

		if( firstColumn->at( 1 ) == "+" ){
			leftDis = contig1->GetLength() - atof( firstColumn->at( 3 ).c_str() );
			firstOri = "+";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is + +, contig ori should be + +
				result += "+\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();
			}
			else{
				// read ori is + -; contig ori should be + -
				result += "+\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
		else{
			leftDis = atoi( firstColumn->at( 3 ).c_str() ) + firstColumn->at( 4 ).length();
			firstOri = "-";
			if( secondColumn->at( 1 ) == "+" ){
				// read ori is - +; contig ori should be - +
				result += "-\t" + secondColumn->at( 2 ) + "\t+\t";
				secondOri = "+";
				rightDis = atoi( secondColumn->at( 3 ).c_str() ) + secondColumn->at( 4 ).length();

			}
			else{
				// read ori is - -; contig ori should be - -
				result += "-\t" + secondColumn->at( 2 ) + "\t-\t";
				secondOri = "-";
				rightDis = contig2->GetLength() - atof( secondColumn->at( 3 ).c_str() );
			}
		}
	}

	dis = (int) (currentLib->GetMean() - leftDis - rightDis);
	result += itos( dis ) + "\t" + itos( currentLib->GetStd() ) + "\n";
	if( pos1 > pos2 ){
		// change two contigs
		int tempPos = pos1;
		pos1 = pos2;
		pos2 = tempPos;
		string tempOri = firstOri;
		if( secondOri == "+" )
			firstOri = "-";
		else
			firstOri = "+";
		if( tempOri == "+" )
			secondOri = "-";
		else
			secondOri = "+";

	}
	// create a single pet object
	SinglePet *newPet = new SinglePet( id, pos1, firstOri, pos2, secondOri, dis, currentLib->GetStd() );

	// save to map
	int firstID = pos1;
	int secondID = pos2;
	if( firstOri == "-" )
		firstID = -firstID;
	if( secondOri == "-" )
		secondID = -secondID;

	// insert to pet map
	currentLib->InsertPair( firstID, secondID, newPet );

	return result;
}



bool less_singlePet( SinglePet*& p1, SinglePet*& p2 ){
	if( p1->GetStartID() < p2->GetStartID() ) 
		return true;
	else if( p1->GetStartID() > p2->GetStartID() )
		return false;
	else{
		if( p1->GetEndID() < p2->GetEndID() )
				return true;
			else if( p1->GetEndID() > p2->GetEndID() )
				return false;
		else{
			if( p1->GetStartOri() < p2->GetStartOri() )
				return true;
			else if( p1->GetStartOri() > p2->GetStartOri() )
				return false;
			else{
				if( p1->GetEndOri() < p2->GetEndOri() )
					return true;
				else if( p1->GetEndOri() > p2->GetEndOri() )
					return false;
				else
					return p1->GetDistance() < p2->GetDistance();
			}// finish checking second contig id
		} // finish checking first contig orientation
	}// finish checking first contig id
}


// start bundle
int MapConverter::Bundle(){
	ofstream clusterWriter( (Configure::OUTPUT_FOLDER + "clusters").c_str() );
	if( clusterWriter == NULL ){
		cout<<"error opening clusters file!"<<endl;
		return -1;
	}
	
	string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
	clusterWriter.write( head.c_str(), head.length() );

	// traverse all group
	map<pair<int, int>, multiset<SinglePet*, less_distance>*>::iterator mapIter = m_singlePetsMap->begin();
	while( mapIter != m_singlePetsMap->end() ){
		// deal with current PET set
		string clusterString = BundleCluster( mapIter->second ); 
		clusterWriter.write( clusterString.c_str(), clusterString.length() );

		// proceed to next one
		delete mapIter->second;
		mapIter++;
	}

	m_singlePetsMap->clear();
	clusterWriter.close();
	return 1;
}

// start bundle
int MapConverter::BundleMultiLib( list<PetLibrary*> *libs ){
	for( list<PetLibrary*>::iterator libIter = libs->begin(); libIter != libs->end(); libIter++ )
	{
		PetLibrary *currentLib = *libIter;

		ofstream clusterWriter( (Configure::OUTPUT_FOLDER + "clusters_" + currentLib->GetNameWithoutPath() ).c_str() );
		if( clusterWriter == NULL ){
			cout<<"error opening clusters file!"<<endl;
			return -1;
		}
	
		string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
		clusterWriter.write( head.c_str(), head.length() );

		// traverse all group
		map<pair<int, int>, multiset<SinglePet*, lessDistance>*> *pairs = currentLib->GetSinglePetsMap();
		map<pair<int, int>, multiset<SinglePet*, lessDistance>*>::iterator mapIter = pairs->begin();
		while( mapIter != pairs->end() ){
			// deal with current PET set
			string clusterString = BundleClusterMultiLib( mapIter->second, currentLib ); 
			clusterWriter.write( clusterString.c_str(), clusterString.length() );

			// proceed to next one
			delete mapIter->second;
			mapIter++;
		}

		pairs->clear();
		clusterWriter.close();

		//cout<<"Total number of clusters: "<<currentLib->GetEdges()->size()<<endl;
	}
	
	return 1;
}

// bundle a certain cluster
// return all the cluster string
string MapConverter::BundleCluster( multiset<SinglePet*, less_distance> *group ){
	PET *result;
	string finalCluster = "";

	if( group->empty() )
		return "";

	while( !group->empty() ){
		// find the median single pet
		int median = group->size() / 2 + 1;

		multiset<SinglePet*, less_distance>::iterator iter = group->begin();
		int counter = 1;
		while( counter < median ){
			iter++;
			counter++;
		}
		SinglePet *medianPet = *iter;

		// initialize formula
		double p = 0;
		double q = 0;
		double newMean = 0;
		double newStd = 0;
		int size = 0;
		
		// calculate distance bound
		int lowerbound = medianPet->GetDistance() - Configure::STD_TIMES * medianPet->GetStd();
		int upperbound = medianPet->GetDistance() + Configure::STD_TIMES * medianPet->GetStd();
		int startID = medianPet->GetStartID();
		int startOri = medianPet->GetStartOri();
		int endID = medianPet->GetEndID();
		int endOri = medianPet->GetEndOri();

		// get all edges in same bundle
		iter = group->begin();
		while( iter != group->end() && (*iter)->GetDistance() < upperbound ){
			if( (*iter)->GetDistance() > lowerbound && (*iter)->GetDistance() < upperbound ){
				// in the bundle
				size++;

				// calculate the value
				p += (*iter)->GetDistance() / pow( (double)(*iter)->GetStd(), 2 );
				q += 1 / pow( (double)(*iter)->GetStd(), 2 );

				// remove that list
				delete (*iter);
				group->erase( iter++ );
			}
			else{
				iter++;
			}
		}

		// calculate new mean and std
		newMean = p / q;
		newStd = 1 / pow( q, 0.5 );
		// consider overlap between contigs
		if( newMean < 0 )
			newMean += (-Configure::MIN_DIS);
		

		// add to graph
		result = new PET( m_graph->GetContig( startID ), startOri, m_graph->GetContig( endID ), endOri, (int) newMean, (int) newStd, size );
		finalCluster += result->ToString() + "\n";
		if( size >= Configure::CLUSTER_THRESHOLD
			&& newMean + Configure::STD_TIMES * newStd >= 0
			&& newMean - Configure::STD_TIMES * newStd <= Configure::LIB_MEAN + Configure::LIB_STD * Configure::STD_TIMES )
			m_graph->AddEdge( result );
		else
			delete result;
	}

	return finalCluster;
}

// bundle a certain cluster
string MapConverter::BundleClusterMultiLib( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib ){
	PET *result;
	string finalCluster = "";

	if( group->empty() )
		return "";

	while( !group->empty() ){
		// find the median single pet
		int median = group->size() / 2 + 1;

		multiset<SinglePet*, lessDistance>::iterator iter = group->begin();
		int counter = 1;
		while( counter < median ){
			iter++;
			counter++;
		}
		SinglePet *medianPet = *iter;

		// initialize formula
		double p = 0;
		double q = 0;
		double newMean = 0;
		double newStd = 0;
		int size = 0;
		
		// calculate distance bound
		int lowerbound = medianPet->GetDistance() - Configure::STD_TIMES * medianPet->GetStd();
		int upperbound = medianPet->GetDistance() + Configure::STD_TIMES * medianPet->GetStd();
		int startID = medianPet->GetStartID();
		int startOri = medianPet->GetStartOri();
		int endID = medianPet->GetEndID();
		int endOri = medianPet->GetEndOri();

		// get all edges in same bundle
		iter = group->begin();
		while( iter != group->end() && (*iter)->GetDistance() < upperbound ){
			if( (*iter)->GetDistance() > lowerbound && (*iter)->GetDistance() < upperbound ){
				// in the bundle
				size++;

				// calculate the value
				p += (*iter)->GetDistance() / pow( (double)(*iter)->GetStd(), 2 );
				q += 1 / pow( (double)(*iter)->GetStd(), 2 );

				// remove that list
				delete (*iter);
				group->erase( iter++ );
			}
			else{
				iter++;
			}
		}

		// calculate new mean and std
		newMean = p / q;
		newStd = 1 / pow( q, 0.5 );
		// consider overlap between contigs
		if( newMean < 0 )
			newMean += (-Configure::MIN_DIS);
		

		// add to graph
		/*********/
		result = new PET( m_graph->GetContig( startID ), startOri, m_graph->GetContig( endID ), endOri, (int) newMean, (int) newStd, size );
		finalCluster += result->ToString() + "\n";
		if( size >= Configure::CLUSTER_THRESHOLD
			&& newMean + Configure::STD_TIMES * newStd >= 0
			&& newMean - Configure::STD_TIMES * newStd <= lib->GetMean() + lib->GetStd() * Configure::STD_TIMES ){
				m_graph->AddEdgeMultiLib( result );
				lib->AddCluster( result );
				//cout<<result->ToString()<<endl;
		}
		else
			delete result;
	}

	return finalCluster;
}

// analyze opera edge file
int MapConverter::AnalyzeOpera( string fileName ){
	ifstream edgeReader( fileName.c_str() );

	if( edgeReader == NULL )
	{
		cout<<"error reading opera edge file: "<<fileName<<endl;
		return -1;
	}

	string line;
	vector<string> *currentPet = new vector<string>;
	while( getline( edgeReader, line ) ){
		Split( line, "\t", currentPet );

		if( currentPet->size() == 7 ){
			// add to graph
			Contig *firstContig = m_graph->GetContig( m_graph->GetContigIndex( (*currentPet)[ 0 ] ) );
			Contig *secondContig = m_graph->GetContig( m_graph->GetContigIndex( (*currentPet)[ 2 ] ) );
			int firstOri = PLUS;
			if( (*currentPet)[ 1 ] == "-" )
				firstOri = MINUS;
			int secondOri = PLUS;
			if( (*currentPet)[ 3 ] == "-" )
				secondOri = MINUS;

			PET *result = new PET( firstContig, firstOri, secondContig, secondOri,
				atoi( (*currentPet)[ 4 ].c_str() ), atoi( (*currentPet)[ 5 ].c_str() ),
				atoi( (*currentPet)[ 6 ].c_str() ) );
			if( atoi( (*currentPet)[ 6 ].c_str() ) >= Configure::CLUSTER_THRESHOLD ){
				if( result->GetDis() + Configure::STD_TIMES * result->GetStd() > 0  
					&& result->GetDis() - Configure::STD_TIMES <= Configure::LIB_MEAN + Configure::LIB_STD * Configure::STD_TIMES )
					m_graph->AddEdge( result );
			}
			else{
				delete result;
			}
		}
		currentPet->clear();
	}
	delete currentPet;
	return 1;
}
