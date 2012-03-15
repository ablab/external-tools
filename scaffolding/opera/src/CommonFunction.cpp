#include "CommonFunction.h"


// convert int to string
string itos( int i ){
	string result;
	stringstream buf;
	buf<<i;
	buf>>result;
	buf.clear();
	return result;
}

void Split( string line, string label, vector<string> *column ){
	column->clear();
	char * cstr, *p;
	cstr = new char [line.length()+1];
	strcpy (cstr, line.c_str());

	p=strtok (cstr, label.c_str() );
	while (p!=NULL)
	{
		column->push_back( string( p ) );
		p=strtok(NULL, label.c_str() );
	}

	delete[] cstr;  
}

void Split( string line, string label, list<string> *column ){
	int pos = line.find( label );
	int id = 0;
	while( pos != -1 ){
		column->push_back( line.substr( 0, pos ) );
		line = line.substr( pos + 1, line.length() - pos - 1 );
		pos = line.find( label );
	}
	if( line != "" )
		column->push_back( line );
}


// get the opposite orientation 
int GetOppositeOri( int ori ){
	if( ori == PLUS )
		return MINUS;
	else
		return PLUS;
}

// print thousand comma format of a number
string PrintNumber( double num ){
	string result;
	char buffer [50];
	if( num > 1000000000 ){
		// gb
		//printf( "%.2f Gb", num/1000000000 );
  	sprintf (buffer, "%.2f Gb", num/1000000000 );
	}
	else if( num > 1000000 ){
		// mb
		//printf( "%.2f Mb", num/1000000 );
		sprintf (buffer, "%.2f Mb", num/1000000 );
	}
	else if( num > 1000 ){
		// kb
		//printf( "%.2f Kb", num/1000 );
		sprintf (buffer, "%.2f Kb", num/1000 );
	}
	else{
		// bp
		//printf( "%d bp", num );
		sprintf (buffer, "%d bp", (int)num );
	}
	
	result = string( buffer );
	return result;
}

// print thousand comma format of a number to file
void PrintNumberToFile( FILE *file, double num ){
	if( num > 1000000000 ){
		// gb
		fprintf( file, "%.2f Gb", num/1000000000 );
	}
	else if( num > 1000000 ){
		// mb
		fprintf( file, "%.2f Mb", num/1000000 );
	}
	else if( num > 1000 ){
		// kb
		fprintf( file, "%.2f Kb", num/1000 );
	}
	else{
		// bp
		fprintf( file, "%d bp", (int)num );
	}
}

// check if it is a number (coverage)
bool IsNumber( string s ){
	for( int i = 0; i < s.size(); i++ ){
		char x = s.at( i );
		if( ( x >= '0' && x <= '9' ) || x == '.' )
			continue;
		else
			return false;
	}

	return true;
}

