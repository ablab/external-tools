#include "Contig.h"

Contig::Contig(void)
{
	m_leftEdges = new list<PET*>;
	m_rightEdges = new list<PET*>;
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_scaffold = toScaffoldString();
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_emptyList = new list< pair<Contig*, int> >;
	m_isInAR = false;

	m_name = "xxx";
	m_id = -1;

	m_leftEdgesMultiLib = new list<PET*>;
	m_rightEdgesMultiLib = new list<PET*>;
}

void Contig::Initialize(){
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_isInAR = false;
}


Contig::Contig( string name, double length, double cov ){
	this->m_name = name;
	this->m_length = length;
	this->m_cov = cov;
	this->m_ori = PLUS;

	m_leftEdges = new list<PET*>;
	m_rightEdges = new list<PET*>;
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_scaffold = toScaffoldString();
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_emptyList = new list< pair<Contig*, int> >;
	m_isInAR = false;

	m_leftEdgesMultiLib = new list<PET*>;
	m_rightEdgesMultiLib = new list<PET*>;
}

Contig::~Contig(void)
{
	DeleteEdges( m_leftEdges );
	delete m_leftEdges;
	DeleteEdges( m_rightEdges );
	delete m_rightEdges;
	delete m_emptyList;

	delete m_leftEdgesMultiLib;
	delete m_rightEdgesMultiLib;
}

void Contig::DeleteEdges( list<PET*> *edges ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		if( !(*iter)->IfVisited() ){
			(*iter)->VisitEdge();
			iter = edges->erase( iter );
		}
		else{
			PET* temp = *iter;
			iter = edges->erase( iter );
			delete temp;
		}
	}
}

// set contig name
void Contig::SetName( string name ){
	m_name = name;
}

// get contig name
string Contig::GetName(){
	return m_name;
}

// set contig length
void Contig::SetLength( double length ){
	m_length = length;
}

// get contig length
double Contig::GetLength(){
	return m_length;
}

// set contig coverage
void Contig::SetCov( double cov ){
	m_cov = cov;
}

// get contig coverage
double Contig::GetCov(){
	return m_cov;
}

// set contig id
void Contig::SetID( int id ){
	m_id = id;
}

// get contig ID
int Contig::GetID(){
	return m_id;
}


// add a pet cluster of this contig
void Contig::AddEdge( PET *p ){	
	// check if it is left edge
	if( ( p->GetEndContig() == this && p->GetEndContigOri() == PLUS )
		|| ( p->GetStartContig() == this && p->GetStartContigOri() == MINUS ) ){
			m_leftEdges->push_back( p );
	}
	else{
		m_rightEdges->push_back( p );
	}
}

// add a pet cluster of this contig, from all libraries
void Contig::AddEdgeMultiLib( PET *p ){
	// check if it is left edge
	if( ( p->GetEndContig() == this && p->GetEndContigOri() == PLUS )
		|| ( p->GetStartContig() == this && p->GetStartContigOri() == MINUS ) ){
			m_leftEdgesMultiLib->push_front( p );
			if( p->GetEndContig() == this )
				p->SetEdgeIterMultiLib( m_leftEdgesMultiLib->begin(), LEFT, END );
			else
				p->SetEdgeIterMultiLib( m_leftEdgesMultiLib->begin(), LEFT, START );
	}
	else{
		m_rightEdgesMultiLib->push_front( p );
		if( p->GetEndContig() == this )
			p->SetEdgeIterMultiLib( m_rightEdgesMultiLib->begin(), RIGHT, END );
		else
			p->SetEdgeIterMultiLib( m_rightEdgesMultiLib->begin(), RIGHT, START );
	}
}

// remove a pet cluster of this contig, from all libraries
void Contig::RemoveEdgeMultiLib( PET *p, int contigPos ){
	if( p->GetEdgeTypeMultiLib( contigPos ) == LEFT ){
		m_leftEdgesMultiLib->erase( p->GetEdgeIterMultiLib( contigPos ) );
	}
	else
		m_rightEdgesMultiLib->erase( p->GetEdgeIterMultiLib( contigPos ) );
}

// set the list iterator of current contig
void Contig::SetListPos( list<Contig*>::iterator iter ){
	m_listPos = iter;
}

// set the border list iterator of current contig
void Contig::SetBorderListPos( list<Contig*>::iterator iter ){
	m_borderListPos = iter;
	m_isBorderContig = true;
}

// generate the scaffold string
string Contig::toScaffoldString(){
	char buffer[ 1000 ];
	int n = sprintf( buffer, "%s\tBE\t%.0f\t%d", m_name.c_str(), m_length, m_gapSize );
	string result( buffer );
	return result; 
}

// check if current contig is a singleton
bool Contig::IsSingleton(){
	if( m_leftEdges->empty() && m_rightEdges->empty() )
		return true;
	else
		return false;
}

// get the super contig(scaffold) id for this contig
int Contig::GetScaffoldID(){
	return m_scaffoldID;
}

list<PET*> * Contig::GetLeftEdges(){
	return m_leftEdges;
}

list<PET*> * Contig::GetRightEdges(){
	return m_rightEdges;
}

list<Contig*>::iterator Contig::GetListPos(){
	return m_listPos;
}

list<Contig*>::iterator Contig::GetBorderListPos(){
	return m_borderListPos;
}

// check if current contig has any edge
bool Contig::HasEdge(){
	return (!m_leftEdges->empty() || !m_rightEdges->empty());
}

// check if current contig is a border contig
bool Contig::IsBorderContig(){
	return m_isBorderContig;
}

// get the left distance
double Contig::GetLeftDis(){
	return m_leftDistance;
}

// get the right distance
double Contig::GetRightDis(){
	return m_rightDistance;
}

// set the direction of extension
void Contig::SetExtensionOri( int ori ){
	if( m_extentionDirection == NEITHER || ori == NEITHER )
		m_extentionDirection = ori;
	else if( ori == BOTH )
		m_extentionDirection = BOTH;
	else if( ori != m_extentionDirection )
		m_extentionDirection = BOTH;
}

// get the direction of extension
int Contig::GetExtensionOri(){
	return m_extentionDirection;
}

// label this contig as in subgraph
void Contig::SetInSubgraph(){
	m_ifInSubgraph = true;
}

// check if this contig is in subgraph
bool Contig::isInSubgraph(){
	return m_ifInSubgraph;
}

// check if this contig has left edges
bool Contig::HasLeftEdge(){
	return !m_leftEdges->empty();
}

// check if this contig has right edges
bool Contig::HasRightEdge(){
	return !m_rightEdges->empty();
}

// set orientation in scaffold
void Contig::SetOri( int ori ){
	m_ori = ori;
}
// get orientation in scaffold
int Contig::GetOri(){
	return m_ori;
}

// set start position with and without gap
void Contig::SetStartPosition( double pos, int gap ){
	m_startPosition = pos;
	m_startPositionWithGap = pos + (double)gap;
}

// get start position without gap
double Contig::GetStartPosition(){
	return m_startPosition;
}

// get start position with gap
double Contig::GetStartPositionWithGap(){
	return m_startPositionWithGap;
}

// get the iterator in unassigned node set
list< pair<Contig*, int> >::iterator Contig::GetUnassignedNodeIter( int type ){
	return m_unassignedNodeIter[ type ];
}

// set unassigned node iterator
void Contig::SetUnassignedNodeIter( int type, list< pair<Contig*, int> >::iterator iter ){
	m_unassignedNodeIter[ type ] = iter;
	m_isUnassignedNode[ type ] = true;
}

// check if it is an unassigned node with certain orientation
bool Contig::IsUnassignedNode( int type ){
	return m_isUnassignedNode[ type ];
	/*if( m_unassignedNodeIter[ type ] != m_endOfList )
		return true;
	return false;*/
}

// set iter to end()
void Contig::ClearUnassignedNodeIter( int type ){
	//m_unassignedNodeIter[ type ] = m_endOfList;
	m_isUnassignedNode[ type ] = false;
}

// set if it is in active region
void Contig::SetIfInAR( bool ar ){
	m_isInAR = ar;
}

// check if it is in active region
bool Contig::IfInAR(){
	return m_isInAR;
}		

// check if it has happy dangling edges in right direction
bool Contig::IfHasHappyDE( int ori ){
	list<PET*> *edges;
	if( ori == PLUS )
		edges = m_rightEdges;
	else
		edges = m_leftEdges;

	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++){
		if( (*iter)->IsDE() && !(*iter)->IfUnhappy() && (*iter)->isInSubgraph() )
			return true;
	}

	return false;
}

void Contig::SetEndOfList( list< pair<Contig*, int> >::iterator iter ){
	m_endOfList = iter;
}

// check the number of valid left Edge for certain orientation
// invalid edges include all dangling edges and unhappy edges
int Contig::GetNumOfValidLeftEdges( int ori ){
	int result = 0;
	list<PET*> *edges;
	if( ori == PLUS )
		edges = m_leftEdges;
	else
		edges = m_rightEdges;

	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() && !(*iter)->IsDE() && !(*iter)->IfUnhappy() )
			result++;
	}

	return result;
}

// generate scaffold String with certain orientation
string Contig::GenScaffoldString( int ori ){
	string result = "";

	if( ori == PLUS ){
		list<string> *contigString = new list<string>;
		Split( m_scaffold, "\t", contigString );
		result = "";
		list<string>::reverse_iterator iter = contigString->rbegin();
		iter++;
		while( iter != contigString->rend() ){
			result = (*iter) + "\t" + result;
			iter++;
		}
		delete contigString;
		return result;
	}

	// generate the minus scaffold string
	list<string> *lines = new list<string>;
	Split( m_scaffold, "\n", lines );
	bool firstLine = true;
	for( list<string>::reverse_iterator iter = lines->rbegin(); iter != lines->rend(); iter++ ){
		string line = *iter;
		vector<string> *contigString = new vector<string>;
		Split( line, "\t", contigString );
		string tempString;
		if( firstLine )
			firstLine = false;
		else{
			result += contigString->at( 3 ) + "\n";
		}

		if( contigString->at( 1 ) == "BE" )
			tempString = contigString->at( 0 ) + "\tEB\t" + contigString->at( 2 ) + "\t";
		else
			tempString = contigString->at( 0 ) + "\tBE\t" + contigString->at( 2 ) + "\t";
		result += tempString;
		delete contigString;
	}

	//result += "0\n";

	delete lines;
	return result;
}

void Contig::SetLeftDis( double d ){
	m_leftDistance = d;
}
	
void Contig::SetRightDis( double d ){
	m_rightDistance = d;
}

void Contig::SetScaffoldID( int id ){
	m_scaffoldID = id;
}

// set scaffold string
void Contig::SetScaffoldString( string sca ){
	m_scaffold = sca;
}

void Contig::SetNotInSubgraph(){
	m_ifInSubgraph = false;
}

// get original contig orientation in scaffold
int Contig::GetContigOriInScaffold( string contigName ){
	list<string> lines;
	Split( m_scaffold, "\n", &lines );
	for( list<string>::iterator iter = lines.begin(); iter != lines.end(); iter++ ){
		vector<string> line;
		Split( *iter, "\t", &line );
		if( line.at( 0 ) == contigName ){
			if( line.at( 1 ) == "BE" )
				return PLUS;
			else
				return MINUS;
		}
	}
	return -1;
}

// reverse left and right edges
void Contig::ReverseEdges(){
	list<PET*> *temp = m_leftEdges;
	m_leftEdges = m_rightEdges;
	m_rightEdges = temp;
}

// get orientation of complex contig
string Contig::GetOriOfFirstContig( int ori ){
	list<string> line;
	Split( m_scaffold, "\t", &line );

	list<string>::iterator iter = line.begin();
	iter++;
	string result = *iter;

	if( ori == PLUS ){
		if( result == "BE" )
			return "BE";
		else
			return "EB";
	}
	else{
		if( result == "BE" )
			return "EB";
		else
			return "BE";
	}
}

string Contig::GetNameOfFirstContig(){
	list<string> line;
	Split( m_scaffold, "\t", &line );

	list<string>::iterator iter = line.begin();
	return *iter;
}

string Contig::GetScaffoldString(){
	return m_scaffold;
}

// check if old contig has the same orientation in new contig
bool Contig::CheckOldContigOri( Contig *c ){
	// get old contig information
	list<string> line;
	Split( c->GetScaffoldString(), "\t", &line );

	list<string>::iterator iter = line.begin();
	string oldContigName = *iter;
	iter++;
	string oldContigOri = *iter;

	// get new contig information
	list<string> lines;
	Split( m_scaffold, "\n", &lines );

	iter = lines.begin();
	while( iter != lines.end() ){
		vector<string> content;
		Split( *iter, "\t", &content );
		if( content.at( 0 ) == oldContigName ){
			if( content.at( 1 ) == oldContigOri )
				return true;
			else
				return false;
		}
		iter++;
	}

	return false;
}

// gap related
// set the id in scaffold
void Contig::SetScaffoldPos( int pos ){
	m_scaffoldPos = pos;
}

// get the id in scaffold
int Contig::GetScaffoldPos(){
	return m_scaffoldPos;
}

// set gap size after this contig
void Contig::SetGap( int gap ){
	m_gapSize = gap;
}

int Contig::GetGap(){
	return m_gapSize;
}

// heuristic related
void Contig::SetStep( int s ){
	m_step = s;
}

int Contig::GetStep(){
	return m_step;
}

// check if this contig has multiple libraries edge
bool Contig::HasMultiLibEdge()
{
	if( m_leftEdgesMultiLib->empty() && m_rightEdgesMultiLib->empty() )
		return false;
	else
		return true;
}

// get the left edge of multiple libraries
list<PET*>* Contig::GetLeftEdgeMultiLib()
{
	return m_leftEdgesMultiLib;
}

// get the right edge of multiple libraries
list<PET*>* Contig::GetRightEdgeMultiLib()
{
	return m_rightEdgesMultiLib;
}