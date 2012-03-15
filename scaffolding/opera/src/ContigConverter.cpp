#include "ContigConverter.h"

ContigConverter::ContigConverter(void)
{
	myContigs = new vector<Contig*>;
	m_repeatContigs = new list<Contig*>;
	m_smallContigs = new list<Contig*>;
}

ContigConverter::~ContigConverter(void)
{
	myContigs->clear();
	delete myContigs;
	DeleteContigs( m_repeatContigs );
	delete m_repeatContigs;
	DeleteContigs( m_smallContigs );
	delete m_smallContigs;
}

void ContigConverter::DeleteContigs( list<Contig*> *contigs ){
	list<Contig*>::iterator iter = contigs->begin();
	while( iter != contigs->end() ){
		Contig* temp = *iter;
		iter = contigs->erase( iter );
		delete temp;
	}
}

// Convert contig file
// return: -1 if failed
int ContigConverter::ConvertContigFile( string fileName, Graph *graph ){
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Reading contig file unsuccessfully"<<endl;
		return -1;
	}

	if( Configure::FILE_FORMAT == FASTA ){
		if( Configure::FILE_TYPE == VELVET ){
			// analyze velvet file
			if( AnalyzeVelvet( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == SOAP ){
			// analyze soapdenovo file
			if( AnalyzeSoapDenovo( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == SOAP_CONTIG ){
			// analyze soapdenovo contig format
			if( AnalyzeSoapDenovoContig( &contigReader, myContigs ) == -1 )
				return -1;
		}
		else if( Configure::FILE_TYPE == NORMAL ){
			// analyze normal fasta file
			if( AnalyzeNormalFasta( &contigReader, myContigs ) == -1 )
				return -1;
		}
	}
	else if( Configure::FILE_FORMAT == STATISTIC ){
		// read the contig file directly in opera format
		if( AnalyzeStatistic( &contigReader, myContigs ) == -1 )
			return -1;
	}

	// if need to remove repeat, remove them
	if( Configure::FILTER_REPEAT ){
		FilterRepeat( myContigs, graph );
	}
	else{
		// save the contigs directly and remove the small contigs
		RemoveSmallContigs( myContigs, graph );
	}

	// generate the ID map of all contigs
	graph->GenerateIDMap();

	// print repeat and small contigs
	if( PrintContigs( m_repeatContigs, Configure::OUTPUT_FOLDER + "repeatContigs" ) == -1 )
		return -1;

	if( PrintContigs( m_smallContigs, Configure::OUTPUT_FOLDER + "smallContigs" ) == -1 )
		return -1;

	// output the opera format contig file
	if( Configure::FILE_FORMAT != STATISTIC ){
		if( graph->OutputContigs( Configure::OUTPUT_FOLDER + "contigs" ) == -1 )
			return -1;
	}

	return 1;
}

// analyze velvet file
// return -1 if failed
int ContigConverter::AnalyzeVelvet( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 0;

	while( getline( *contigReader, line ) ){
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				contigs->push_back( newContig );
			}

			// start new contig
			name = line.substr( 1, line.length() - 1 );
			int pos = line.find_last_of( "_" );
			cov = atof( line.substr( pos + 1, line.length() - pos - 1 ).c_str() );
			length = 0;;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	contigs->push_back( newContig );

	return 1;
}

// analyze soapdenovo file
int ContigConverter::AnalyzeSoapDenovo( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line.length() == 0 )
			continue;
			
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				contigs->push_back( newContig );
			}

			// start new contig
			vector<string> *content = new vector<string>;
			Split( line, " ", content );
			name = line.substr( 1, content->at( 0 ).length() - 1 );
			cov = atof( content->at( 1 ).c_str() );
			length = 0;
			content->clear();
			delete content;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	contigs->push_back( newContig );

	return 1;
}


// analyze soapdenovo file
int ContigConverter::AnalyzeSoapDenovoContig( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line.length() == 0 )
			continue;
			
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				contigs->push_back( newContig );
			}

			// start new contig
			vector<string> *content = new vector<string>;
			Split( line, " ", content );
			name = line.substr( 1, content->at( 0 ).length() - 1 );
			Split( line, "_", content );
			cov = atof( content->at( 1 ).c_str() );
			length = 0;
			content->clear();
			delete content;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	contigs->push_back( newContig );

	return 1;
}

// analyze normal fasta file
int  ContigConverter::AnalyzeNormalFasta( ifstream *contigReader, vector<Contig*> *contigs )
{
	string line;
	
	string name;
	double length = 0;
	double cov = 1;

	while( getline( *contigReader, line ) ){
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( length > 0 ){
				Contig *newContig = new Contig( name, length, cov );
				contigs->push_back( newContig );
			}

			// start new contig
			char nameChar[ 500 ];
			sscanf( line.c_str(), "%s", nameChar );
			name = string( nameChar );
			name = name.substr( 1, name.length() );
			length = 0;
		}
		else{
			length += line.length();
		}
	}

	// save last contig
	Contig *newContig = new Contig( name, length, cov );
	contigs->push_back( newContig );
	return 1;
}

// filter repeat
void ContigConverter::FilterRepeat( vector<Contig*> *contigs, Graph *graph ){
	vector<Contig*> usedContigs;		// save the non-repeat long contigs

	double sum = 0;
	double num = 0;
	for( int i = 0; i < contigs->size(); i++ ){
		sum += contigs->at( i )->GetCov() * contigs->at( i )->GetLength();
		num += contigs->at( i )->GetLength();
	}

	double average = sum / num;
	double threshold = average * Configure::REPEAT_THRESHOLD;

#ifdef DEBUG
	cout<<"average coverage is "<<average<<endl;
#endif

	vector<Contig*>::iterator contigIter = contigs->begin();
	while( contigIter != contigs->end() ){
		if( (*contigIter)->GetCov() > threshold ){
			// remove repeat contigs
			m_repeatContigs->push_back( *contigIter );
		}
		else if( (*contigIter)->GetLength() < Configure::CONTIG_SIZE_THERSHOLD ){
			// remove small contigs
			m_smallContigs->push_back( *contigIter );
		}
		else{
			usedContigs.push_back( *contigIter );
		}
		contigIter++;
	}

	// save to graph
	graph->SetContigNum( usedContigs.size() );
	for( int i = 0; i < usedContigs.size(); i++ )
		graph->AddContig( usedContigs.at( i ) );
}

// analyze opera's contig file
int ContigConverter::AnalyzeStatistic( ifstream *contigReader, vector<Contig*> *contigs ){
	string line;

	while( getline( *contigReader, line ) ){
		int pos = line.find( "\t" );
		if( pos == -1 ){
			return -1;
		}

		string name = line.substr( 0, pos );
		string length = line.substr( pos + 1, line.length() - pos - 1 );
		Contig *newContig = new Contig( name, atoi( length.c_str() ), 1 );
		contigs->push_back( newContig );
	}

	return 1;
}

// remove small contigs only
void ContigConverter::RemoveSmallContigs( vector<Contig*> *contigs, Graph *graph ){
	vector<Contig*> usedContigs;		// save the non-repeat long contigs

	vector<Contig*>::iterator contigIter = contigs->begin();
	while( contigIter != contigs->end() ){
		if( (*contigIter)->GetLength() < Configure::CONTIG_SIZE_THERSHOLD ){
			// remove small contigs
			m_smallContigs->push_back( *contigIter );
		}
		else{
			usedContigs.push_back( *contigIter );
		}
		contigIter++;
	}

	// save to graph
	graph->SetContigNum( usedContigs.size() );
	for( int i = 0; i < usedContigs.size(); i++ )
		graph->AddContig( usedContigs.at( i ) );
}

// print list of contigs
int ContigConverter::PrintContigs( list<Contig*> *contigs, string fileName ){
	ofstream contigWriter( fileName.c_str() );

	if( contigWriter == NULL ){
		cout<<"error writing "<<fileName<<endl;
		return -1;
	}
	
	string head = "Contig ID\tLength\tCoverage\n";
	contigWriter.write( head.c_str(), head.length() );

	list<Contig*>::const_iterator iter;
	for( iter = contigs->begin(); iter != contigs->end(); iter++ ){
		string result = (*iter)->GetName(); 
		char buffer[500];  
		//sprintf( buffer, "\t%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() + 0.5 );
		sprintf( buffer, "\t%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() );
		string tempString( buffer );
		result += tempString;
		contigWriter.write( result.c_str(), result.length() );
	}

	contigWriter.close();

	return 1;
}

