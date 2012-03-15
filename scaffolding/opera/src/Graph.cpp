#include "Graph.h"

Graph::Graph(void)
{
	//m_contigNameMap = new map<string, int>;
	m_contigNameHashMap = new hash_map<const char*, int, hash<const char*>, eqName>;
	m_numberOfContigs = 0;
	m_numberOfEdges = 0;

	m_contigsList = new list<Contig*>;		
	m_scaffoldsList = new list<Contig*>;
	m_borderContigsList = new list<Contig*>;

	ifInitializeContig = false;

	m_contigID = new map<int, Contig*>();
}

Graph::~Graph(void)
{
	if( ifInitializeContig )
		delete[] m_contigsArray;
	DeleteContigs( m_contigsList );
	delete m_contigsList;
	DeleteContigs( m_scaffoldsList );
	delete m_scaffoldsList;
	m_borderContigsList->clear();
	delete m_borderContigsList;

	//m_contigNameMap->clear();
	//delete m_contigNameMap;
	m_contigNameHashMap->clear();
	delete m_contigNameHashMap;

	m_contigID->clear();
	delete m_contigID;
}

void Graph::DeleteContigs( list<Contig*> *contigs ){
	list<Contig*>::iterator iter = contigs->begin();
	while( iter != contigs->end() ){
		Contig* temp = *iter;
		iter = contigs->erase( iter );
		delete temp;
	}
}

// set total number of contigs
void Graph::SetContigNum( int n ){
	m_contigsArray = new Contig*[ n ];
	ifInitializeContig = true;
}

// Add non-repeat contigs
void Graph::AddContig( Contig *c ){
	// add to array
	m_contigsArray[ m_numberOfContigs ] = c;
	c->SetID( m_numberOfContigs );
	//m_contigNameMap->insert( pair<string, int>(c->GetName(), c->GetID() ) );
	m_contigNameHashMap->insert( pair<const char*, int>(c->GetName().c_str(), c->GetID() ) );
	m_numberOfContigs++;

	// add to list
	m_contigsList->push_front( c );
	c->SetListPos( m_contigsList->begin() );
}

// initialize contig list for a new library
void Graph::InitializeContig( int libNum ){
	// do not need to initialize for the first library
	if( libNum == 1 )
		return;

	m_contigNameHashMap->clear();

	m_numberOfContigs = 0;
	list<Contig*>::iterator scaIter = this->m_scaffoldsList->begin();
	//cout<<"number of contigs: "<<m_scaffoldsList->size()<<endl;
	while( scaIter != m_scaffoldsList->end() ){
		this->m_contigsList->push_front( *scaIter );
		(*scaIter)->SetListPos( m_contigsList->begin() );
		(*scaIter)->Initialize();
		m_contigsArray[ m_numberOfContigs ] = *scaIter;
		(*scaIter)->SetID( m_numberOfContigs );
		m_contigNameHashMap->insert( pair<const char*, int>((*scaIter)->GetName().c_str(), (*scaIter)->GetID() ) );
		
		m_numberOfContigs++;
		scaIter = m_scaffoldsList->erase( scaIter );
	}
}

// output the contigs into a new file
// return -1 if failed
int Graph::OutputContigs( string fileName ){
	ofstream contigWriter( fileName.c_str() );

	if( contigWriter == NULL )
	{
		cout<<"error writing contig file"<<endl;
		return -1;
	}

	string head = "Contig ID\tLength\tCoverage\n";
	contigWriter.write( head.c_str(), head.length() );
	
	list<Contig*>::const_iterator iter;
	for( iter = m_contigsList->begin(); iter != m_contigsList->end(); iter++ ){
		string line = (*iter)->GetName() + "\t";
		char buffer[30];
		sprintf( buffer, "%.0f\t%.0f\n", (*iter)->GetLength(), (*iter)->GetCov() );
		string tempString( buffer );
		line += tempString;
		//line += ( itos( (*iter)->GetLength() )+ "\n");
		contigWriter.write( line.c_str(), line.length() );
	}

	return 1;
}

// get contig index according to contig name
// return -1 if not find
int Graph::GetContigIndex( string contigName ){
	//map<string, int>::const_iterator pos = m_contigNameMap->find( contigName );
	hash_map<const char*, int, hash<const char*>, eqName>::iterator pos = m_contigNameHashMap->find( contigName.c_str() );
	if( pos == m_contigNameHashMap->end() )
		return -1;
	else
		return pos->second;

	
}

// get contig in specific index
Contig * Graph::GetContig( int pos ){
	return m_contigsArray[ pos ];
}

// add edge to graph
void Graph::AddEdge( PET *p ){
	p->SetID( m_numberOfEdges++ );
	p->GetStartContig()->AddEdge( p );
	p->GetEndContig()->AddEdge( p );
}

// add edges from all libraries to graph
void Graph::AddEdgeMultiLib( PET *p ){
	//p->SetID( m_numberOfEdges++ );
	p->GetStartContig()->AddEdgeMultiLib( p );
	p->GetEndContig()->AddEdgeMultiLib( p );
}

// generate border contig and singletons
void Graph::GenerateBorderAndScaffold(){
	int threshold = Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD;

	list<Contig*>::iterator iter = m_contigsList->begin();
	while( iter != m_contigsList->end() ){
		if( (*iter)->IsSingleton() ){
			// save to scaffold list
			m_scaffoldsList->push_back( *iter );
			iter = m_contigsList->erase( iter );
		}
		else if( (*iter)->GetLength() >= threshold ){
			// save to border contig List
			m_borderContigsList->push_front( *iter );
			(*iter)->SetBorderListPos( m_borderContigsList->begin() );
			iter++;
		}
		else
			iter++;
	}

	//cout<<"there are "<<this->m_contigsList->size()<<" contigs with edge\n";
}

// update the graph using results of subgraph
void Graph::UpdateGraph( list<Contig*> *subgraph, ScaffoldResult **scaffold ){
	// deal with each contigs
	list<Contig*> *tempListOfContig = new list<Contig*>;	// save all contigs needed to be deleted in subgraph

	list<Contig*>::iterator contigIter = subgraph->begin();
	while( contigIter != subgraph->end() ){
		//cout<<"scaffold\n"<<scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString()<<endl;
		//cout<<"contig:\n"<<(*contigIter)->GetScaffoldString()<<endl;

		bool isFirstContig = false;
		bool sameOri = false;
		// check if need to reverse left and right edges
		sameOri = SameOri( *contigIter, scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString() );
		/*if( (*contigIter)->GetName()== "1347")
		{
			cout<<"1347 same? "<<sameOri<<endl;
			cout<<"scaffold\n"<<scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString()<<endl;
			cout<<"contig:\n"<<(*contigIter)->GetScaffoldString()<<endl;
		}*/
		if( scaffold[ (*contigIter)->GetScaffoldID() ]->GetID() == -1 ){
			if( !sameOri )
				(*contigIter)->ReverseEdges();
			//if( IfSameOriInScaffoldAndEdge( *contigIter, MINUS, scaffold[ (*contigIter)->GetScaffoldID() ] ) )
			//	(*contigIter)->ReverseEdges();
			isFirstContig = true;
			scaffold[ (*contigIter)->GetScaffoldID() ]->SetID( (*contigIter)->GetID() );
			(*contigIter)->SetScaffoldString( scaffold[ (*contigIter)->GetScaffoldID() ]->GetScaffoldString() );
			(*contigIter)->SetLength( scaffold[ (*contigIter)->GetScaffoldID() ]->GetLength() );
			(*contigIter)->SetCov( scaffold[ (*contigIter)->GetScaffoldID() ]->GetCov() );
		}

		//cout<<(*contigIter)->GetScaffoldString()<<endl;

		// traverse each edge
		DeleteEdges( *contigIter, (*contigIter)->GetLeftEdges(),
			isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );

		DeleteEdges( *contigIter, (*contigIter)->GetRightEdges(),
			isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );

		// for multiple libraries
		// update the edges in other libraries
		DeleteEdgesMultiLib( *contigIter, (*contigIter)->GetLeftEdgeMultiLib(),
			isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );

		DeleteEdgesMultiLib( *contigIter, (*contigIter)->GetRightEdgeMultiLib(),
			isFirstContig, scaffold[ (*contigIter)->GetScaffoldID() ], sameOri );

		// delete contig
		if( isFirstContig ){
			// update border contig pointer
			if( !(*contigIter)->IsBorderContig() ){
				m_borderContigsList->push_front( *contigIter );
				(*contigIter)->SetBorderListPos( m_borderContigsList->begin() );
			}
			contigIter++;
		}
		else{
			// delect current contig and remove the pointer
			m_contigsList->erase( (*contigIter)->GetListPos() );
			if( (*contigIter)->IsBorderContig() )
				m_borderContigsList->erase( (*contigIter)->GetBorderListPos() );
			m_contigsArray[ (*contigIter)->GetID() ] = NULL;
			//Contig* temp = *contigIter;
			tempListOfContig->push_back( *contigIter );
			contigIter = subgraph->erase( contigIter );
			//delete temp;
		}
	}

	// delete contigs
	contigIter = tempListOfContig->begin();
	while( contigIter != tempListOfContig->end() ){
		Contig *temp = *contigIter;
		delete temp;
		contigIter = tempListOfContig->erase( contigIter );
	}
	delete tempListOfContig;

	// check remaining new contigs if they are singletons
	contigIter = subgraph->begin();
	while( contigIter != subgraph->end() ){
		(*contigIter)->SetNotInSubgraph();
		if( !(*contigIter)->HasEdge() ){
			// remove to singleton list
			m_scaffoldsList->push_back( *contigIter );
			m_contigsArray[ (*contigIter)->GetID() ] = NULL;
			if( (*contigIter)->IsBorderContig() )
				m_borderContigsList->erase( (*contigIter)->GetBorderListPos() );
			m_contigsList->erase( (*contigIter)->GetListPos() );
		}

		contigIter++;
	}
}

// check if contig c has the same orientation as in scaffold
bool Graph::SameOri( Contig *c, string scaffold ){
	// get old contig information
	list<string> line;
	Split( c->GetScaffoldString(), "\t", &line );

	list<string>::iterator iter = line.begin();
	string oldContigName = *iter;
	iter++;
	string oldContigOri = *iter;
		
	/*if( oldContigName == "1593" )
	{
		cout<<"contig:\n"<<c->GetScaffoldString()<<endl;
		cout<<"scaffold:\n"<<scaffold<<endl;
	}*/

	// get new contig information
	list<string> lines;
	Split( scaffold, "\n", &lines );

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

// release the memory of edges
void Graph::DeleteEdges( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		// if this edge is within the subgraph
		if( (*iter)->IfInSubgraph() ){
			if( !(*iter)->IfVisited() ){
				// if it is the first time, update the status of edge
				(*iter)->VisitEdge();
				iter = edges->erase( iter );
			}
			else{
				// the second time, delete
				PET *temp = *iter;
				iter = edges->erase( iter );
				delete temp;
			}
		}
		else if( c->GetScaffoldID() != -1){
			//cout<<"Edge is needed to be updated"<<endl;
			//cout<<"before: "<<(*iter)->ToString()<<endl;
			//cout<<"if same: "<<ifSame<<endl;
			//cout<<c->GetName()<<endl;
			// if this edge connects outside the subgraph, update it
			if( isFirst ){
				// meet the first contig in a scaffold, replace it
				if( scaffold->GetID() == -1 ){
					scaffold->SetID( c->GetID() );
					// update length of contig
					c->SetLength( scaffold->GetLength() );
				}
			}

			// replace edges
			(*iter)->SetNotInSubgraph();
			(*iter)->SetDE( false );
			int pos = (*iter)->GetPositionOfContig( c );
			int ori = (*iter)->GetOrientationOfContig( c );
			int preOri = ori;
			if( !isFirst )
				(*iter)->ReplaceContig( pos, m_contigsArray[ scaffold->GetID() ] );

			if( !ifSame ){
				if( ori == MINUS )
					ori = PLUS;
				else
					ori = MINUS;
				
				(*iter)->SetOri( pos, ori );
			}

			// recalculate the distance
#ifdef CHECK
			int oldDis = (*iter)->GetDis();
#endif
			if( pos == START ){
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int)c->GetLeftDis() );
				}
			}
			else{
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetLeftDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
			}

			//cout<<"after: "<<(*iter)->ToString()<<endl;

#ifdef CHECK
			if( (*iter)->GetDis() > oldDis ){
				//cout<<"new distance > old distance"<<endl;
				//cout<<"edge is "<<(*iter)->GetOriString()<<endl<<endl;
			}
#endif

			if( !isFirst ){
				// add current edge to super contig
				m_contigsArray[ scaffold->GetID() ]->AddEdge( *iter );

				// remove current edge from current contig
				iter = edges->erase( iter );
			}
			else
				iter++;
		}
	}
}

// release the memory of edges
void Graph::DeleteEdgesMultiLib( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		// if this edge is within the subgraph
		if( (*iter)->GetStartContig()->isInSubgraph() && (*iter)->GetEndContig()->isInSubgraph() ){
			if( !(*iter)->IfVisited() ){
				// if it is the first time, update the status of edge
				(*iter)->VisitEdge();
				iter = edges->erase( iter );
			}
			else{
				// the second time, delete
				PET *temp = *iter;
				// delete from library
				(*iter)->DeleteFromLib();

				iter = edges->erase( iter );
	
				delete temp;
			}
		}
		else if( c->GetScaffoldID() != -1){

			//cout<<"scaffold:\n"<<scaffold->GetScaffoldString()<<endl;
			//cout<<"contig:\n"<<c->GetScaffoldString()<<endl;
			//cout<<"Edge is needed to be updated"<<endl;
			//cout<<"before: "<<(*iter)->ToString()<<endl;
			
			// if this edge connects outside the subgraph, update it
			/*if( isFirst ){
				// meet the first contig in a scaffold, replace it
				if( scaffold->GetID() == -1 ){
					scaffold->SetID( c->GetID() );
					// update length of contig
					c->SetLength( scaffold->GetLength() );
				}
			}*/

			// replace edges
			int pos = (*iter)->GetPositionOfContig( c );
			int ori = (*iter)->GetOrientationOfContig( c );
			int preOri = ori;
			if( !isFirst )
				(*iter)->ReplaceContig( pos, m_contigsArray[ scaffold->GetID() ] );

			if( !ifSame ){
				if( ori == MINUS )
					ori = PLUS;
				else
					ori = MINUS;
				
				(*iter)->SetOri( pos, ori );
			}

			// recalculate the distance
#ifdef CHECK
			int oldDis = (*iter)->GetDis();
#endif
			if( pos == START ){
				if( ori == PLUS  ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int)c->GetLeftDis() );
				}
			}
			else{
				if( ori == PLUS ){
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetLeftDis() );
				}
				else{
						(*iter)->SetDis( (*iter)->GetDis() - (int) c->GetRightDis() );
				}
			}

			//cout<<"after: "<<(*iter)->ToString()<<endl<<endl;
			
#ifdef CHECK
			if( (*iter)->GetDis() > oldDis ){
				//cout<<"new distance > old distance"<<endl;
				//cout<<"edge is "<<(*iter)->GetOriString()<<endl<<endl;
			}
#endif

			if( !isFirst ){
				// add current edge to super contig
				m_contigsArray[ scaffold->GetID() ]->AddEdgeMultiLib( *iter );

				// remove current edge from current contig
				iter = edges->erase( iter );
			}
			else
				iter++;
		}
	}
}

// check if current scaffold has the same orientation in scaffold and edge
bool Graph::IfSameOriInScaffoldAndEdge( Contig *c, int ori, ScaffoldResult *s ){
	/*string oriString = "BE";
	if( ori == MINUS )
		oriString = "EB";*/

	string contigName = c->GetNameOfFirstContig();
	string oriString = c->GetOriOfFirstContig( ori );

	vector<string> lines;
	Split( s->GetScaffoldString(), "\n", &lines );

	for( int i = 0; i < lines.size(); i++ ){
		vector<string> line;
		Split( lines.at( i ), "\t", &line );
		if( line.at( 0 ) == contigName ){
			if( line.at( 1 ) == oriString )
				return true;
			else 
				return false;
		}
	}

	return true;
}


// find subgraph
list<Contig*> * Graph::FindSubgraph( int &numOfContig, int &numOfBorderContig ){
	list<Contig*> *subgraph;
		
	Contig *c;

	// check if there is any border contigs
	if( !m_borderContigsList->empty() ){
		c = *( m_borderContigsList->begin() );
		if( c->HasRightEdge() )
			subgraph = TraverseSubgraph( *(m_borderContigsList->begin()), RIGHT, 
					numOfContig, numOfBorderContig );
		else
			subgraph = TraverseSubgraph( *(m_borderContigsList->begin()), LEFT, 
					numOfContig, numOfBorderContig );
	}
	else if( !m_contigsList->empty() ){
		c = *( m_contigsList->begin() );
		if( c->HasRightEdge() )
			subgraph = TraverseSubgraph( *(m_contigsList->begin()), RIGHT, 
						numOfContig, numOfBorderContig );
		else
			subgraph = TraverseSubgraph( *(m_contigsList->begin()), LEFT, 
					numOfContig, numOfBorderContig );
	}

	// check if need to add edges
	for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		list<PET*> *edges = (*iter)->GetLeftEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() )
					(*edgeIter)->SetInSubgraph();
			}
		}

		edges = (*iter)->GetRightEdges();
		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( (*edgeIter)->GetStartContig()->isInSubgraph() && (*edgeIter)->GetEndContig()->isInSubgraph() ){
				if( !(*edgeIter)->isInSubgraph() )
					(*edgeIter)->SetInSubgraph();
			}
		}
	}

	return subgraph;
}

int Graph::GetNumOfContigs(){
	return m_numberOfContigs;
}

// traverse and find subgraph starting from contig c 
// c is the start contig; direction is the extension direction
list<Contig*> * Graph::TraverseSubgraph( Contig *c, int direction, int &numOfContig, int &numOfBorderContig ){
	list<Contig*> *subgraph = new list<Contig*>;

	numOfContig = 0;
	numOfBorderContig = 0;
		
	list<Contig*> possibleContigs;
	possibleContigs.push_back( c );
	c->SetExtensionOri( direction );

	while( !possibleContigs.empty() ){
		// get the first contig
		Contig *currentContig = *possibleContigs.begin();
		possibleContigs.pop_front();

		// label this contig as in subgraph if it is not
		if( !currentContig->isInSubgraph() ){
			currentContig->SetInSubgraph();
			subgraph->push_back( currentContig );
			numOfContig++;
			if( currentContig->IsBorderContig() )
				numOfBorderContig++;
		}

		// traverse edges
		if( currentContig->GetExtensionOri() == RIGHT )
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs );
		else if( currentContig->GetExtensionOri() == LEFT )
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs );
		else if( currentContig->GetExtensionOri() == BOTH ){
			AddEdgesToSubgraph( currentContig, currentContig->GetLeftEdges(), &possibleContigs );
			AddEdgesToSubgraph( currentContig, currentContig->GetRightEdges(), &possibleContigs );
		}
	}

	return subgraph;
}


// add list of edges into subgraph
// add edges to subgraph, and all contigs except c to contigList
void Graph::AddEdgesToSubgraph( Contig *c, list<PET*> *edges, list<Contig*> *contigList ){
	list<PET*>::iterator iter;
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() )
			continue;
		else{
			(*iter)->SetInSubgraph();
			(*iter)->SetDE( false );
			Contig *otherContig = (*iter)->GetOtherContig( c );
			int otherPos = (*iter)->GetPositionOfContig( otherContig );
			int otherOri = (*iter)->GetOrientationOfContig( otherContig );
			//int ext;
			if( otherContig->GetLength() < Configure::LIB_MEAN + Configure::LIB_STD * Configure::STD_TIMES ){
				// short contig, just add and extend to both orientation
				otherContig->SetExtensionOri( BOTH );
			}
			else{
				// border contig, need consider orientation
				if( ( otherPos == START && otherOri == MINUS ) 
				|| ( otherPos == END && otherOri == PLUS ) ){
					// new contig should extend to left
					otherContig->SetExtensionOri( LEFT );
				}
				else
					otherContig->SetExtensionOri( RIGHT );
			}

			contigList->push_back( otherContig );
		}
	}
}

// check if current graph still has edges
bool Graph::HasEdges(){
	return !m_contigsList->empty();
}

// output final scaffolds to files
// return -1 if failed
int Graph::OutputScaffolds( string fileName ){
	ofstream scafWriter( fileName.c_str() );
	if( scafWriter == NULL ){
		cout<<"error writing opera_scaffold.txt file!"<<endl;
		return -1;
	}

	// sort scaffold
	m_scaffoldsList->sort( compare_scaffold );

	int num = 1;
	for( list<Contig*>::iterator iter = m_scaffoldsList->begin(); iter != m_scaffoldsList->end(); iter++ ){
		string scaf = (*iter)->GenScaffoldString( PLUS ) + "0\n";
		char buffer[ 200 ];
		sprintf( buffer, ">%s%d\tlength: %.0f\tcov: %.1f\n", Configure::SCAFFOLD_PREFEX.c_str(), num, (*iter)->GetLength(), (*iter)->GetCov() );
		string head( buffer );
		scafWriter.write( head.c_str(), head.length() );
		scafWriter.write( scaf.c_str(), scaf.length() );
		num++;
	}
	scafWriter.close();
	return 1;
}

// comparison, not case sensitive.
bool compare_scaffold( Contig *first, Contig *second)
{
	if( first->GetLength() > second->GetLength() )
		return true;
	else 
		return false;
}

// generate the id map of all contigs
void Graph::GenerateIDMap()
{
	if( Configure::FILE_TYPE == VELVET ){
		// handle velvet format
		for( int i = 0; i < m_numberOfContigs; i++ ){
			vector<string> *name = new vector<string>;
			Split( m_contigsArray[ i ]->GetName(), "_", name );
			m_contigID->insert( pair<int, Contig*>( atoi( name->at( 1 ).c_str() ), m_contigsArray[ i ] ) );
		}
	}
}

// get the contig using velvet id
Contig* Graph::GetContigUsingID( int id )
{
	return (*m_contigID)[ id ];
}

// clear the edges number
void Graph::ClearEdge(){
	m_numberOfEdges = 0;
}