#include "opera.h"
// optimal scaffolding


opera::opera(void)
{
	m_graph = new Graph();
	m_activeRegion = new list<Contig*>;		
	m_happyDanglingEdges = new list<PET*>;
	m_unhappyDanglingEdges = new list<PET*>;
	m_unassignedNodes = new list< pair<Contig*, int> >;
	m_unhappyEdgeString = new list<string>;
	m_visitedTree = new Tree();

	// for multi libraries
	m_libraries = new list<PetLibrary*>;
}

// initialize for another run
void opera::Initialize(){
	delete m_visitedTree;
	m_visitedTree = new Tree();
	m_activeRegion->clear();;		
	m_happyDanglingEdges->clear();
	m_unhappyDanglingEdges->clear();
	m_unassignedNodes->clear();
}

opera::~opera(void)
{
	delete m_graph;
	delete m_activeRegion;
	delete m_happyDanglingEdges;
	delete m_unhappyDanglingEdges;
	delete m_unassignedNodes;
	delete m_unhappyEdgeString;
	delete m_visitedTree;
	fclose( logFile );

	delete m_libraries;
}

// clear scaffold variables
void opera::Clear(){
	m_activeRegion->clear();		
	m_happyDanglingEdges->clear();
	m_unhappyDanglingEdges->clear();
	m_unassignedNodes->clear();
	m_visitedTree->Clear();
}

// check if the names have confliction
// return true if has confliction
// return false if not
bool opera::CheckNameConfliction(){
	ifstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Reading contig file unsuccessfully"<<endl;
		return true;
	}
		
	string line;					
			
	getline( contigReader, line );
	contigReader.close();
	
	// check the contig name
	if( line.substr( 1, Configure::SCAFFOLD_PREFEX.length() ) == Configure::SCAFFOLD_PREFEX )
		return true;
	else
		return false;
}

// check if all files exist
bool opera::CheckFileExist(){
	// check contig file
	fstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"ERROR: Contig file \""<<Configure::CONTIG_FILE<<"\" does not exist!"<<endl;
		return false;
	}
	
	// check mapping file
	// check contig file
	vector<string> *files = new vector<string>;
	Split( Configure::MAP_FILE, ",", files );
	fstream mapReader( files->at( 0 ).c_str() );
	files->clear();
	delete files;

	if( mapReader == NULL )
	{
		cout<<"ERROR: Mapping file \""<<Configure::MAP_FILE<<"\" does not exist!"<<endl;
		return false;
	}
	return true;
}

// start doing scaffolding
// return -1 if failed
int opera::StartScaffold(){
	int finishedScaffold = 0;
	while( m_graph->HasEdges() ){
		// find subgraph
		int numOfContig = 0;
		int numOfBorderContig = 0;
		list<Contig*> *subgraph = m_graph->FindSubgraph( numOfContig, numOfBorderContig );

		// set all contigs unassigned nodes iterator to null
#ifdef LOG
		fprintf( logFile, "start one subgraph\n" );
#endif
		for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
			(*iter)->ClearUnassignedNodeIter( PLUS );
			(*iter)->ClearUnassignedNodeIter( MINUS );
			(*iter)->SetStartPosition( 1, 0 );

#ifdef LOG
			// print 
			/*fprintf( logFile, "%s\tid is%d\n", (*iter)->GetName().c_str(), (*iter)->GetID() );
			list<PET*>::iterator edgeIter = (*iter)->GetLeftEdges()->begin();
			while( edgeIter != (*iter)->GetLeftEdges()->end() ){
				if( (*edgeIter)->IfInSubgraph() )
					fprintf( logFile, "%s\n", (*edgeIter)->GetOriString().c_str() );
				edgeIter++;
			}
			edgeIter = (*iter)->GetRightEdges()->begin();
			while( edgeIter != (*iter)->GetRightEdges()->end() ){
				if( (*edgeIter)->IfInSubgraph() )
					fprintf( logFile, "%s\n", (*edgeIter)->GetOriString().c_str() );
				edgeIter++;
			}*/
#endif

#ifdef LOG
			//fprintf( logFile, "%d\tlength: %.0f\n", (*iter)->GetID(), (*iter)->GetLength() );
#endif
		}

#ifdef LOG
		fprintf( logFile, "%d contigs in subgraph\n", subgraph->size() );
#endif

		// select starting contig
		Contig *startContig = NULL;
		int startContigOri;
		FindStartingContig( subgraph, numOfContig, numOfBorderContig, startContig, startContigOri );

		// scaffolding
		bool ifSucceed = false;
		int maxOfUnhappyEdges = -1;
		// create a scaffold containing startContig
		Clear();
		PartialScaffold *currentScaffold  = CreatScaffold( startContig, startContigOri );
		PartialScaffold *head = currentScaffold;
		GenUnassigndNodes( currentScaffold );

		int numberOfScaffold = 0;
		struct timeval t_start,t_end;		// record the time for one run
		do{
			gettimeofday( &t_start, NULL );
			// starting from 0
			maxOfUnhappyEdges++;

			//pair<string, int> newPair( currentScaffold->GetScaffoldString(), 
			//	currentScaffold->GetNumOfUnhappyEdges() );
			m_visitedTree->IfExist( currentScaffold->GetARString(),
					currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() );
			//m_visitedScaffolds->insert( newPair );
			// check if current scaffold has more than required unhappy edges
			if( currentScaffold->GetNumOfUnhappyEdges() > maxOfUnhappyEdges ){
				maxOfUnhappyEdges++;
				continue;
			}

			numberOfScaffold = Scaffolding( currentScaffold, maxOfUnhappyEdges );
			Clear();		// clear the map tree
			if( numberOfScaffold != -1 ){
				// succeed
				ifSucceed = true;
			}
			else{
				//Clear();
				// set all contigs unassigned nodes iterator to null
				for( list<Contig*>::iterator iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
					//(*iter)->SetEndOfList( m_unassignedNodes->end() );
					(*iter)->ClearUnassignedNodeIter( PLUS );
					(*iter)->ClearUnassignedNodeIter( MINUS );
					(*iter)->SetStartPosition( 1, 0 );
				}

				DeleteAllScaffolds( currentScaffold, NULL );
				currentScaffold  = CreatScaffold( startContig, startContigOri );
				GenUnassigndNodes( currentScaffold );
			}
			gettimeofday( &t_end, NULL );
			if( !ifSucceed && Configure::ABORT ){
				if( (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec)/1000000.0 > 1800 ){
					// if running time is longer than 30 mins, then remind user if he wants to quit
					cout<<"The time required might be quite long, we suggest you to increase the value of \"cluster_threshold\" in configuration file"<<endl;
					cout<<"Quiting the program...Please run again with new parameters"<<endl;
					//delete currentScaffold;
					DeleteAllScaffolds( currentScaffold, NULL );
					return -1;
				}
			}
		} while( !ifSucceed );

		// update graph
#ifdef LOG
		fprintf( logFile, "%d scaffolds\n", numberOfScaffold );
		fprintf( logFile, "initialize %d ScaffoldResult array\n", numberOfScaffold );
#endif
		ScaffoldResult **result = new ScaffoldResult*[ numberOfScaffold ];   
		GenerateResults( currentScaffold, result, numberOfScaffold );
		// calculate gap sizes
		m_graph->UpdateGraph( subgraph, result );

#ifdef LOG
		fprintf( logFile, "finish one subgraph\n\n" );
#endif

#ifdef DEBUG
		PrintScaffoldResult( result, numberOfScaffold );
		cout<<endl;
#endif

#ifdef DEBUG
		finishedScaffold++;
		cout<<"finish "<<itos( finishedScaffold )<<" subgraph"<<endl;
#endif
		subgraph->clear();
		for( int i = 0; i < numberOfScaffold; i++ )
			delete result[ i ];
		delete[] result;
		// delete scaffold
		DeleteAllScaffolds( currentScaffold, NULL );
		delete subgraph;
	}
	
	return 1;
}

// delete all scaffold
inline void opera::DeleteAllScaffolds( PartialScaffold *p, PartialScaffold *head ){
	while( p != head ){
		PartialScaffold *temp = p;
		p = p->GetParent();
		delete temp;
	}
}

bool SortStartPoints( const StartPoint *s1, const StartPoint *s2 ){
	if( s1->GetNumOfUnhappyEdges() == s2->GetNumOfUnhappyEdges() ){
		return s1->GetContig()->GetLength() > s2->GetContig()->GetLength();
	}
	else
		return s1->GetNumOfUnhappyEdges() < s2->GetNumOfUnhappyEdges();
}

// select the starting contig
void opera::FindStartingContig( list<Contig*> *subgraph, int numOfContig, int numOfBorderContig,
							   Contig *&startContig, int &startContigOri ){
	StartPoint **allStart;
	int num = 0;
	int type;
	if( numOfBorderContig != 0 ){
		allStart = new StartPoint*[ numOfBorderContig * 2 ];
		type = BORDER;
	}
	else{
		allStart = new StartPoint*[ numOfContig * 2 ];
		type = ALL;
	}

	list<Contig*>::const_iterator iter;
	for( iter = subgraph->begin(); iter != subgraph->end(); iter++ ){
		if( ( (type == BORDER ) && ( (*iter)->IsBorderContig() ) ) 
			|| (type == ALL ) ){
				allStart[ num++ ] = new StartPoint( *iter, PLUS );
				allStart[ num++ ] = new StartPoint( *iter, MINUS );
		}
	}

	sort( allStart, allStart + num, SortStartPoints );
	startContig = allStart[ 0 ]->GetContig();
	startContigOri = allStart[ 0 ]->GetOri();

	for( int i = 0; i < num; i++ )
		delete allStart[ i ];

	delete[] allStart;
}	

// create a scaffold containing a start contig
PartialScaffold* opera::CreatScaffold( Contig *c, int ori ){

	PartialScaffold *result = new PartialScaffold;
		
	// add a contig to active region
	c->SetOri( ori );
	c->SetStartPosition( 1, 0 );
	result->AddNode( c );
	m_activeRegion->push_back( c );
	
	// add edges to dangling edges
	if( ori == PLUS ){
		// add all right edges
		AddHappyDEToFirstScaffold( c->GetRightEdges(), result );
		AddUnhappyDEToFirstScaffold( c->GetLeftEdges(), result, c );
	}
	else{
		// add all left edges
		AddHappyDEToFirstScaffold( c->GetLeftEdges(), result );
		AddUnhappyDEToFirstScaffold( c->GetRightEdges(), result, c );
	}

	GenARString( result );
	GenUDEString( result );

	// generate unassigned node
		
	return result;
}

// add happy dangling edges to first scaffolds
void opera::AddHappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p ){
	list<PET*>::const_iterator iter = edges->begin();
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() ){
			p->AddAddedHappyDE( *iter );
			m_happyDanglingEdges->push_front( *iter );
			(*iter)->SetHappyDEIter( m_happyDanglingEdges->begin() );
		}
	}
}

// add unhappy dangling edges to first scaffolds
// c is the first contig
void opera::AddUnhappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p, Contig *c ){
	list<PET*>::const_iterator iter = edges->begin();
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() ){
			p->AddAddedUnhappyDE( *iter, (*iter)->GetOtherContig( c ) );
			//cout<<"unhappy due to first contig:"<<(*iter)->GetOriString()<<endl;
			m_unhappyDanglingEdges->push_front( *iter );
			(*iter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		}
	}
}

bool UDESort( PET *&p1, PET *&p2 ){
	return p1->GetID() < p2->GetID();
}

// generate unhappy dangling edges string
void opera::GenUDEString( PartialScaffold *p ){
	if( m_unhappyDanglingEdges->empty() ){
		p->SetUnhappyDEString( "" );
		return;
	}

	// sort unhappy dangling edges
	m_unhappyDanglingEdges->sort( UDESort );

	string result = "";
	for( list<PET*>::const_iterator iter = m_unhappyDanglingEdges->begin();
		iter != m_unhappyDanglingEdges->end(); iter++ ){
			result += itos((*iter)->GetID()) + "\t";
	}

	p->SetUnhappyDEString( result );
}

// generate active region string
void opera::GenARString( PartialScaffold *p ){
	string result = "";

	for( list<Contig*>::const_iterator iter = m_activeRegion->begin(); 
		iter != m_activeRegion->end(); iter++ ){
			result += itos( (*iter)->GetID() ) + "\t";
			if( (*iter)->GetOri() == PLUS )
				result += "+\t";
			else
				result += "-\t";
	}

	p->SetARString( result );
}

// scaffolding
// return -1 if failed
// otherwise, return the number of scaffold
int opera::Scaffolding( PartialScaffold *&currentScaffold, int maxOfUnhappyEdges ){
	int numberOfScaffold = 1;

	// initialize the tree saving all visited scaffolds
	m_visitedTree->Clear();

#ifdef DEBUG
	if( maxOfUnhappyEdges == 6 )
		cout<<"6"<<endl;
#endif
	// while it is not root partial scaffold or it is not finished
	while( currentScaffold->GetParent() != NULL || !currentScaffold->IfBreaked() ){
		if( currentScaffold->IfBreaked() ){
			// if it is already breaked, trace back
			TraceBack( currentScaffold );

			//if( m_activeRegion->empty() && m_happyDanglingEdges->empty() ){
			if( currentScaffold->IfEmptyScaffold() ){
#ifdef SPLIT
				cout<<"minus 1 for number of scaffolds"<<endl;
				if( currentScaffold->IfMannualEmptyScaffold() )
					cout<<"return from a manual splited one"<<endl;
				else
					cout<<"retrun from a normal splited one"<<endl;
#endif
				numberOfScaffold--;
			}
			continue;
		}

		if( !m_activeRegion->empty() ){
			// select an unassigned node and add to the end
			Contig *contig = NULL;
			int contigOri = 0;
			// active region is not empty
			if( !m_unassignedNodes->empty() ){
				// if there still unused unassigned nodes
				// try to add a contig
				if( !AddContigToScaffold( currentScaffold, maxOfUnhappyEdges, contig, contigOri ) ){
					// if adding is not successful, trace back
					// new scaffold is not valid, need to trace back
					TraceBack( currentScaffold );
					if( currentScaffold->IfEmptyScaffold() ){
#ifdef SPLIT
						cout<<"minus 1 for number of scaffolds"<<endl;
						if( currentScaffold->IfMannualEmptyScaffold() )
							cout<<"return from a manual splited one"<<endl;
						else
							cout<<"retrun from a normal splited one"<<endl;
#endif
						numberOfScaffold--;
					}
				}
				else{
					// adding successful, check if new scaffold has been visited yet or not
					// check if new scaffold exist and insert if not
					bool ifExist = m_visitedTree->IfExist( currentScaffold->GetARString(),
						currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() );
					if( ifExist ){
						// no need to continue, trace back
						TraceBack( currentScaffold );
						if( currentScaffold->IfEmptyScaffold() ){
#ifdef SPLIT
							cout<<"minus 1 for number of scaffolds"<<endl;
							if( currentScaffold->IfMannualEmptyScaffold() )
								cout<<"return from a manual splited one"<<endl;
							else
								cout<<"retrun from a normal splited one"<<endl;
#endif
							numberOfScaffold--;
						}	
					}

#ifdef SPLIT
					string logString = "\nmax unhappy edges: " + itos( maxOfUnhappyEdges ) + "\n";
					logString += "new Scaffold\n" + *(currentScaffold->GetARString()) + "\n" + *(currentScaffold->GetDEString()) + "\n";
					fprintf( logFile, logString.c_str() );
#endif
				}
			}
			else{
				// if there is no unassigned nodes
				// check if current scaffold has been break or not
				if( !currentScaffold->IfBreaked() ){
					// if current scaffold is not breaked yet, break it
					currentScaffold->Break();
					PartialScaffold *prePS = currentScaffold;
					BreakScaffold( currentScaffold, maxOfUnhappyEdges );
					if( currentScaffold != prePS ){
						numberOfScaffold++;
#ifdef SPLIT
						cout<<"mannualy break a scaffold"<<numberOfScaffold<<endl;
						cout<<"last contig in last scaffold is "<<prePS->GetAddedContigInAR()->GetName()<<endl;
#endif
					}
				}
			}
		}
		else{
			// active region is empty
			if( m_unhappyDanglingEdges->empty() ){
				// if there is no unhappy dangling edges, find a solution
#ifdef SPLIT
						fprintf( logFile, "return after adding contig\n" );
#endif
				return numberOfScaffold;
			}
			else{
				// if there is some unhappy dangling edges, start a new scaffold
				currentScaffold->SetIfEmptyScaffold( true );
				currentScaffold->Break();
				numberOfScaffold++;
#ifdef SPLIT
				cout<<"split a new scaffold "<<numberOfScaffold<<endl;
				cout<<"last contig in last scaffold is "<<currentScaffold->GetAddedContigInAR()->GetName()<<endl;
#endif
				// need to start a new scaffold
				StartPoint *sp = FindStartOfNewScaffold();
				sp->GetContig()->SetStartPosition( 1, 0 );
				// create a scaffold
				if( !AddContigToScaffold( currentScaffold, maxOfUnhappyEdges, 
							sp->GetContig(), sp->GetOri() ) ){
						// cannot start new scaffold, trace back
						TraceBack( currentScaffold );	// track back to end of last scaffold
						TraceBack( currentScaffold );	// change an end
						numberOfScaffold--;
#ifdef SPLIT
						cout<<"trace back from AddContigToScaffold"<<endl;
#endif
						if( currentScaffold->IfEmptyScaffold() )
							numberOfScaffold--;
				}
				else{
					// create new scaffold successfully
					if( m_visitedTree->IfExist( currentScaffold->GetARString(),
								currentScaffold->GetDEString(), currentScaffold->GetNumOfUnhappyEdges() ) ){
									// if this new scaffold already exist
							TraceBack( currentScaffold );	// track back to end of last scaffold
							TraceBack( currentScaffold );	// change an end
							numberOfScaffold--;
#ifdef SPLIT
						cout<<"trace back from been visited before"<<endl;
#endif
						if( currentScaffold->IfEmptyScaffold() )
							numberOfScaffold--;
					}
				}
				delete sp;
			}
		}
	}

#ifdef DEBUG
	cout<<"total visited partial scaffolds: "<<itos( m_visitedTree->m_totalNum )<<endl;
#endif

	return -1;
}

// Break scaffold
void opera::BreakScaffold( PartialScaffold *&p, int maxUnhappyEdges ){
	PartialScaffold *emptySP = new PartialScaffold();
	emptySP->SetNumOfUnhappyEdges( p->GetNumOfUnhappyEdges() );
	emptySP->SetParent( p );
	emptySP->Break();
	emptySP->SetIfEmptyScaffold( true );
	emptySP->SetIfMannualEmptyScaffold( true );

	int increasedUnhappyEdges = m_happyDanglingEdges->size();
	// set all happy dangling edges as unhappy dangling edges
	if( increasedUnhappyEdges + p->GetNumOfUnhappyEdges() > maxUnhappyEdges ){
		delete emptySP;
		return;
	}

	list<PET*>::iterator edgeIter = m_happyDanglingEdges->begin();
	while( edgeIter != m_happyDanglingEdges->end() ){
		PET *edge = *edgeIter;
		edgeIter = m_happyDanglingEdges->erase( edgeIter );
		m_unhappyDanglingEdges->push_front( edge );
		edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		emptySP->AddRemovedHappyDE( edge );
		emptySP->AddAddedUnhappyDE( edge, edge->GetOtherContig( edge->FindARContig() ) );
	}

	list<Contig*>::iterator contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *contig = *contigIter;
		contigIter = m_activeRegion->erase( contigIter );
		contig->SetIfInAR( false );
		emptySP->AddRemovedContigInAR( contig );
		
		// check if currentContig connects to active region
		if( !CheckActiveRegionConnectivity( contig ) ){
			TraceBack( emptySP );
			return;
		}
	}

	StartPoint *sp = FindStartOfNewScaffold();
	// create a scaffold
	if( !AddContigToScaffold( emptySP, maxUnhappyEdges, sp->GetContig(), sp->GetOri() ) ){
		// cannot start new scaffold
		//PartialScaffold *parent = emptySP->GetParent(); 
		TraceBack( emptySP );
		TraceBack( emptySP );
		delete sp;
		return;
	}

	delete sp;

	p = emptySP;
}


// find start contig of new scaffold
// return start point
StartPoint* opera::FindStartOfNewScaffold(){
	list<StartPoint*> *borderContigs = new list<StartPoint*>;	// border contigs
	list<StartPoint*> *allContigs = new list<StartPoint*>;		// non border contigs
	StartPoint *result = NULL;

	// traverse
	set<Contig*> *visitedContigs = new set<Contig*>;
	list<Contig*> *possibleContigs = new list<Contig*>;
		// start with unhappy dangling edges
	for( list<PET*>::iterator iter = m_unhappyDanglingEdges->begin();
		iter != m_unhappyDanglingEdges->end(); iter++ ){
			Contig *contig = (*iter)->GetUnusedContig();
			if( visitedContigs->insert( contig ).second ){
				// it is the first time for this contig, put into the list
				AddToList( borderContigs, allContigs, contig, PLUS );
				AddToList( borderContigs, allContigs, contig, MINUS );
				possibleContigs->push_back( contig );
			}
	}

	while( !possibleContigs->empty() ){
		Contig *contig = *(possibleContigs->begin());
		possibleContigs->pop_front();

		CheckEdgesForStartPoint( borderContigs, allContigs, contig, contig->GetLeftEdges(), 
			visitedContigs, possibleContigs );
		CheckEdgesForStartPoint( borderContigs, allContigs, contig, contig->GetRightEdges(), 
			visitedContigs, possibleContigs );
	}

	if( !borderContigs->empty() ){
		borderContigs->sort( SortStartPoints );
		result = *(borderContigs->begin());
		borderContigs->pop_front();
	}
	else{
		allContigs->sort( SortStartPoints );
		result = *(allContigs->begin());
		allContigs->pop_front();
	}

	// delete variables
	list<StartPoint*>::iterator iter = borderContigs->begin();
	while( iter != borderContigs->end() ){
		StartPoint *sp = *iter;
		iter = borderContigs->erase( iter );
		delete sp;
	}
	delete borderContigs;
	iter = allContigs->begin();
	while( iter != allContigs->end() ){
		StartPoint *sp = *iter;
		iter = allContigs->erase( iter );
		delete sp;
	}
	delete allContigs;
	delete visitedContigs;
	delete possibleContigs;

	return result;
}

// check edges to find start point
void opera::CheckEdgesForStartPoint( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, 
									Contig *c, list<PET*> *edges, set<Contig*> *visitedContigs,
									list<Contig*> *possibleContigs ){
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		PET *edge = *iter;
		if( edge->isInSubgraph() && !edge->IsDE() && !edge->IfUnhappy() ){
			Contig *otherContig = edge->GetOtherContig( c );
			if( otherContig->isInSubgraph() && visitedContigs->insert( otherContig ).second ){
				// it is the first time for this contig, put into the list
				AddToList( borderContigs, allContigs, otherContig, PLUS );
				AddToList( borderContigs, allContigs, otherContig, MINUS );
				possibleContigs->push_back( otherContig );
			}
		}
	}
}

// add connected contigs to corresponding list
void opera::AddToList( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, Contig *c, int ori ){
	StartPoint *newSP = new StartPoint( c, ori );
	newSP->SetNumOfUnhappyEdges( c->GetNumOfValidLeftEdges( ori ) );
	if( c->IsBorderContig() )
		borderContigs->push_back( newSP );
	else
		allContigs->push_back( newSP );
}

// generate unassigned contig set
void opera::GenUnassigndNodes( PartialScaffold *p ){
	// clear unassigned nodes possible orientation
	for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
		(*iter).first->SetExtensionOri( NEITHER );
	}

	// check new unassigned nodes
	list< pair<Contig*, int> > *possibleContig = new list< pair<Contig*, int> >;	// int is the extension direction
	set<Contig*> *newUnassignedNodes = new set<Contig*>;
	set<PET*> *visitedEdge = new set<PET*>;

	// deal with all contigs in active region
	for( list<Contig*>::iterator iter = m_activeRegion->begin(); iter != m_activeRegion->end(); iter++ ){
		Contig *contig = *iter;
		list<PET*> *edges;
		if( (*iter)->GetOri() == PLUS )
			edges = (*iter)->GetRightEdges();
		else
			edges = (*iter)->GetLeftEdges();

		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			if( !(*edgeIter)->IfUnhappy() && (*edgeIter)->isInSubgraph() && (*edgeIter)->IsDE() ){
				// only consider happy edge
				PET *edge = *edgeIter;
				visitedEdge->insert( edge );
				int ori = edge->GetOrientationOfContig( contig );
				Contig *otherContig = edge->GetOtherContig( contig );
				int otherOri = edge->GetOrientationOfContig( otherContig );
				otherContig->SetStep( 1 );
				int ext;
				if( !otherContig->IsBorderContig() ){
					// it is not border contig, extend both direction
					otherContig->SetExtensionOri( BOTH );
					ext = BOTH;
				}
				else{
					// it is border contig, only extend to one direction
					if( ( ori == contig->GetOri() && otherOri == PLUS ) ||
						( ori != contig->GetOri() && otherOri == MINUS ) ){
							otherContig->SetExtensionOri( LEFT );
							ext = LEFT;
					}
					else{
						otherContig->SetExtensionOri( RIGHT );
						ext = RIGHT;
					}
				}
				pair<Contig*, int> p( otherContig, ext );
				possibleContig->push_back( p );
				newUnassignedNodes->insert( otherContig );
			}
		}
	}

	// start traversing
	while( !possibleContig->empty() ){
		// get the first contig
		pair<Contig*, int> pair = *possibleContig->begin();
		Contig *currentContig = pair.first;
		int ext = pair.second;
		possibleContig->pop_front();

		if( Configure::HEURISTIC ){
			// only consider contigs within several steps
			if( currentContig->GetStep() >= Configure::STEP )
				continue;
		}

		// traverse edges
		if( ext == RIGHT )
			TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
					 newUnassignedNodes, visitedEdge );
		else if( ext == LEFT )
			TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
					 newUnassignedNodes, visitedEdge );
		else if( ext == BOTH ){
			TraverseEdges( currentContig, currentContig->GetLeftEdges(), possibleContig,
					 newUnassignedNodes, visitedEdge );
			TraverseEdges( currentContig, currentContig->GetRightEdges(), possibleContig,
					 newUnassignedNodes, visitedEdge );
		}
	}

	// find removed unassigned nodes and update corresponding variables
	list< pair<Contig*, int> >::iterator pairIter = m_unassignedNodes->begin();
	while( pairIter != m_unassignedNodes->end() ){
		if( (*pairIter).first->GetExtensionOri() == NEITHER ){
			// it is the removed unassigned nodes
			p->AddRemovedUnassignedNodes( (*pairIter).first, (*pairIter).second );
			(*pairIter).first->ClearUnassignedNodeIter( (*pairIter).second );
			pairIter = m_unassignedNodes->erase( (*pairIter).first->GetUnassignedNodeIter( (*pairIter).second ) );
		}
		else
			pairIter++;
	}

	// find added unassigned nodes
	for( set<Contig*>::iterator iter = newUnassignedNodes->begin(); iter != newUnassignedNodes->end(); iter++ ){
		if( !(*iter)->IsUnassignedNode( PLUS ) ){
			// not an unassigned nodes before
			p->AddAddedUnassignedNodes( *iter, PLUS );
			pair<Contig*, int> pair1( *iter, PLUS );
			m_unassignedNodes->push_front( pair1 );
			(*iter)->SetUnassignedNodeIter( PLUS, m_unassignedNodes->begin() );
		}
		if( !(*iter)->IsUnassignedNode( MINUS ) ){
			p->AddAddedUnassignedNodes( *iter, MINUS );
			pair<Contig*, int> pair2( *iter, MINUS );
			m_unassignedNodes->push_front( pair2 );
			(*iter)->SetUnassignedNodeIter( MINUS, m_unassignedNodes->begin() );
		}
	}

	newUnassignedNodes->clear();   delete newUnassignedNodes;
	possibleContig->clear();  delete possibleContig;
	visitedEdge->clear();  delete visitedEdge;
}

// traverse edge and find unassigned node
void opera::TraverseEdges( Contig *c, list<PET*> *edges, list< pair<Contig*, int> > *&possibleContig,
						   set<Contig*> *&unassignedNodes, set<PET*> *&visitedEdges ){
	list<PET*>::iterator iter;
	for( iter = edges->begin(); iter != edges->end(); iter++ ){
		if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() || (*iter)->IsDE() )
			continue;

		if( visitedEdges->find(*iter) != visitedEdges->end() )
			continue;
		else{
			visitedEdges->insert(*iter);
			Contig *otherContig = (*iter)->GetOtherContig( c );
			int otherPos = (*iter)->GetPositionOfContig( otherContig );
			int otherOri = (*iter)->GetOrientationOfContig( otherContig );
			int ext;
			if( !otherContig->IsBorderContig() ){
				// short contig, just add and extend to both orientation
				otherContig->SetExtensionOri( BOTH );
				ext = BOTH;
			}
			else{
				// border contig, need consider orientation
				if( ( otherPos == START && otherOri == MINUS ) 
				|| ( otherPos == END && otherOri == PLUS ) ){
					// new contig should extend to left
					otherContig->SetExtensionOri( LEFT );
					ext = LEFT;
				}
				else{
					otherContig->SetExtensionOri( RIGHT );
					ext = RIGHT;
				}
			}
			otherContig->SetStep( c->GetStep() + 1 );
			pair<Contig*, int> p(otherContig, ext);
			possibleContig->push_back( p );
			unassignedNodes->insert( otherContig );
		}
	}
}


// add one contig to partial scaffold
// if the new scaffold is not valid, return false
// if newContig is NULL, select one contig from unassigned node list
//   or just add newContig to the end of scaffold
bool opera::AddContigToScaffold( PartialScaffold *&p, int maxUnhappyEdge, Contig *newContig, 
								int newContigOri ){

	// generate a new partial scaffold
	PartialScaffold *newScaffold = new PartialScaffold( p );
	newScaffold->SetParent( p );
	newScaffold->SetNumOfUnhappyEdges( p->GetNumOfUnhappyEdges() );
	p = newScaffold;

	// get the added contig
	if( newContig == NULL ){
		pair<Contig*, int> firstPair = *m_unassignedNodes->begin();

		newContig = firstPair.first;
		newContigOri = firstPair.second;
		newContig->SetOri( newContigOri );
		Contig *endContig = m_activeRegion->back();		// the last contig in active region
		newContig->SetStartPosition( endContig->GetStartPosition() + endContig->GetLength(), 0 );
		m_unassignedNodes->erase( newContig->GetUnassignedNodeIter( newContigOri ) );
		newContig->ClearUnassignedNodeIter( newContigOri );
	}
	else{
		newContig->SetOri( newContigOri );
		newContig->SetStartPosition( 1, 0 );
	}

	if( m_activeRegion->empty() )
		newContig->SetStartPosition( 1, 0 );

	// check all edges of new Contig
	CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetLeftEdges(), newContigOri );
	CheckEdgesOfNewContig( newScaffold, newContig, newContig->GetRightEdges(), newContigOri );
	
	// add v to partial scaffold
	newScaffold->AddNode( newContig );
	m_activeRegion->push_back( newContig );
	if( newScaffold->GetParent() != NULL && !newScaffold->GetParent()->IfEmptyScaffold() ){
		// add to parent removed list
		newScaffold->GetParent()->AddRemovedUnassignedNodes( newContig, newContigOri );
	}

	if( !newScaffold->GetParent()->IfEmptyScaffold() ){
		// check if the other orientation exists, remove and add to removed list
		int otherOri = (newContigOri+1) % 2;
		if( newContig->IsUnassignedNode( otherOri ) ){
			// add to current remove list
			newScaffold->AddRemovedUnassignedNodes( newContig, otherOri );
			// remove from unassigned node list
			m_unassignedNodes->erase( newContig->GetUnassignedNodeIter( otherOri ) );
			newContig->ClearUnassignedNodeIter( otherOri );
		}
	}

	// check connectivity and update dangling edges and active region
	if( !CheckScaffold( newScaffold ) )
		return false;

	if( newScaffold->GetNumOfUnhappyEdges() > maxUnhappyEdge )
		return false;

	// generate active region string and unhappy dangling edge string
	GenUDEString( newScaffold );
	GenARString( newScaffold );

	if( !m_activeRegion->empty() )
		GenUnassigndNodes( newScaffold );
	else{
		// clear unassigned nodes
		for( list< pair<Contig*, int> >::iterator iter = m_unassignedNodes->begin(); iter != m_unassignedNodes->end(); iter++ ){
			newScaffold->AddRemovedUnassignedNodes( (*iter).first, (*iter).second );
			(*iter).first->ClearUnassignedNodeIter( (*iter).second );
		}
		m_unassignedNodes->clear();
	}

	return true;
}


// check edges of newly added contig
// ori is the orientation of contig c
void opera::CheckEdgesOfNewContig( PartialScaffold *p, Contig *c, list<PET*> *edges, int ori ){
	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		PET *edge = *iter;
		if( !(*iter)->isInSubgraph() )
			continue;

		if( edge->IsDE() ){
			// current edge is dangling edge
			if( !edge->IfUnhappy() ){
				// current edge was happy, remove from the list
				p->AddRemovedHappyDE( edge );
				m_happyDanglingEdges->erase( edge->GetHappyDEIter() );

				//if it is not happy any more, add the unhappy edge number
				if( !CheckDanglingEdgeHappiness( edge, c, ori ) ){
					p->AddUnhappyEdgeNumber( 1 );
					p->AddAddedUnhappyEdge( edge );
				}
			}
			else{
				// current edge was not happy before, 
				// remove from unhappy dangling edges
				p->AddRmovedUnhappyDE( edge );
				m_unhappyDanglingEdges->erase( edge->GetUnhappyDEIter() );
			}
		}
		else{
			// it was not a dangling edge, check if orientation is correct
			if( !IsLeftEdge( edge, c, ori ) && edge->GetDis() + edge->GetStd() * Configure::STD_TIMES >= -Configure::KMER ){
				// it is happy, set as happy dangling edge
				p->AddAddedHappyDE( edge );
				m_happyDanglingEdges->push_front( edge );
				edge->SetHappyDEIter( m_happyDanglingEdges->begin() );
			}
			else{
				// it is unhappy, set as unhappy dangling edge
				//cout<<"wrong orientation:"<<edge->GetOriString()<<endl;
				p->AddAddedUnhappyDE( edge, edge->GetOtherContig( c ) );
				m_unhappyDanglingEdges->push_front( edge );
				edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
			}
		}
	}
}

// check dangling edge happiness
// newContig is the contig which is not in active region
// ori is the orientation of newContig
bool opera::CheckDanglingEdgeHappiness( PET *edge, Contig *newContig, int ori ){
	// check if edge is left edge( orientation is happy )
	if( !IsLeftEdge( edge, newContig, ori ) )
		return false;

	// check if distance is happy
		// first, find the node in active region
	Contig *activeNode = edge->GetOtherContig( newContig );
	Contig *endContig = *m_activeRegion->rbegin();		// the last contig in active region
	int dis = (int) (newContig->GetStartPosition() - activeNode->GetStartPosition() - activeNode->GetLength());
	if( dis - Configure::KMER > edge->GetDis() + Configure::STD_TIMES * edge->GetStd() )
		return false;
	return true;
}

// check if it is a left edge of current contig with certain orientation
bool opera::IsLeftEdge( PET *edge, Contig *newContig, int ori ){
	int petOri = edge->GetOrientationOfContig( newContig );
	int petPos = edge->GetPositionOfContig( newContig );
	if( ( petOri == ori && petPos == START ) || ( petOri != ori && petPos == END ) )
		return false;
	else
		return true;
}

// check connectivity and update scaffold active region & dangling edges
// return false if it is not connected
bool opera::CheckScaffold( PartialScaffold *p ){
	Contig *endContig = *m_activeRegion->rbegin();		// the last contig in active region

	// check each dangling edge if the distance could be satisfied
	list<PET*>::iterator iter = m_happyDanglingEdges->begin();
	while( iter != m_happyDanglingEdges->end() ){
			PET *edge = *iter;
#ifdef CHECK
			if( edge->GetOriString().find( "contig_872" )!= string::npos && edge->GetOriString().find( "contig_873" ) != string::npos ){
				cout<<edge->GetOriString()<<endl;
				cout<<"new distance is "<<edge->GetDis()<<"  std is "<<edge->GetStd()<<endl;
				cout<<"end contig in active region is "<<endContig->GetScaffoldString()<<endl;
				cout<<"end contig start position is "<<endContig->GetStartPosition()<<endl;
				cout<<"contig in active region is "<<edge->FindARContig()->GetName()<<endl;
				cout<<"contig in active region end position is "<<edge->FindARContig()->GetStartPosition() + edge->FindARContig()->GetLength()<<endl;
			}
#endif
			// get the contig in active region
			Contig *ARContig = edge->FindARContig();
			int dis = (int) (endContig->GetStartPosition() + endContig->GetLength()
				- ARContig->GetStartPosition() - ARContig->GetLength() - 2);
			if( dis -Configure::KMER >= edge->GetDis() + Configure::STD_TIMES * edge->GetStd() ){
				// change to unhappy dangling edge
				iter = m_happyDanglingEdges->erase( edge->GetHappyDEIter() );
				//cout<<"distance not happy:"<<edge->GetOriString()<<endl;
				m_unhappyDanglingEdges->push_front( edge );
				edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
				p->AddRemovedHappyDE( edge );
				p->AddAddedUnhappyDE( edge, edge->GetOtherContig( ARContig ) );
			}
			else
				iter++;
	}

	// update the active region, remove the head nodes which do not have happy dangling edges
	list<Contig*>::iterator contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() ){
		Contig *currentContig = *contigIter;
		if( !currentContig->IfHasHappyDE( currentContig->GetOri() ) ){
			// remove from active region
			contigIter = m_activeRegion->erase( contigIter );
			currentContig->SetIfInAR( false );
			if( currentContig != p->GetAddedContigInAR() )
				p->AddRemovedContigInAR( currentContig );
			else
				p->RemoveOneFromActiveRegion();		// it is the newly added contig, 
													// only minus 1 from number of active region

			// check if currentContig connects to active region
			if( !CheckActiveRegionConnectivity( currentContig ) )
				return false;
			
			//cout<<"remove node from ar successfully: "<<currentContig->GetName()<<endl;
		}
		else
			break;
	}

	// shrink the active region if it is longer than library threshold
	if( m_activeRegion->empty() ){
		return true;
	}

	list<Contig*>::reverse_iterator rIter = m_activeRegion->rbegin();
	double length = 0;
	while( rIter != m_activeRegion->rend() ){
		Contig *currentContig = *rIter;
		length += currentContig->GetLength();
		if( length > Configure::LIB_MEAN + Configure::STD_TIMES * Configure::LIB_STD ){
			// have to remove all contigs before current contig
			rIter++;
			break;
		}
		rIter++;
	}

	if( rIter == m_activeRegion->rend() )
		return true;
	// return to the starting contig
	rIter--;

	contigIter = m_activeRegion->begin();
	while( contigIter != m_activeRegion->end() 
				&& *contigIter != *rIter ){
		Contig *currentContig = *contigIter;
		// remove from active region
		contigIter = m_activeRegion->erase( contigIter );
		currentContig->SetIfInAR( false );

		// set all happy dangling edges to unhappy dangling edges
		list<PET*>::iterator edgeIter = currentContig->GetLeftEdges()->begin();
		while( edgeIter != currentContig->GetLeftEdges()->end() ){
			if( (*edgeIter)->IsDE() && !(*edgeIter)->IfUnhappy() ){
				// change to unhappy dangling edge
				m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
				m_unhappyDanglingEdges->push_front( *edgeIter );
				(*edgeIter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
				p->AddRemovedHappyDE( *edgeIter );
				p->AddAddedUnhappyDE( *edgeIter, (*edgeIter)->GetOtherContig( currentContig ) );
			}
			edgeIter++;
		}
		edgeIter = currentContig->GetRightEdges()->begin();
		while( edgeIter != currentContig->GetRightEdges()->end() ){
			if( (*edgeIter)->IsDE() && !(*edgeIter)->IfUnhappy() ){
				// change to unhappy dangling edge
				m_happyDanglingEdges->erase( (*edgeIter)->GetHappyDEIter() );
				m_unhappyDanglingEdges->push_front( *edgeIter );
				(*edgeIter)->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
				p->AddRemovedHappyDE( *edgeIter );
				p->AddAddedUnhappyDE( *edgeIter, (*edgeIter)->GetOtherContig( currentContig ) );
			}
			edgeIter++;
		}

		if( currentContig != p->GetAddedContigInAR() )
			p->AddRemovedContigInAR( currentContig );
		else
			p->RemoveOneFromActiveRegion();		// it is the newly added contig, 
													// only minus 1 from number of active region

		// check if currentContig connects to active region
		if( !CheckActiveRegionConnectivity( currentContig ) )
			return false;

		//contigIter++;
	}
	
	return true;
}

// check if contig c could connect to active region
// return false if not connected
bool opera::CheckActiveRegionConnectivity( Contig *c ){
	if( m_activeRegion->empty() ){
		//cout<<"return from empty"<<endl;
		return true;
	}

	set<Contig*> *visitedContigs = new set<Contig*>;
	list<Contig*> *possibleContigs = new list<Contig*>;

	possibleContigs->push_back( c );
	visitedContigs->insert( c );

	while( !possibleContigs->empty() ){
		Contig *currentContig = *possibleContigs->begin();
		possibleContigs->pop_front();

		// check left edges
		list<PET*> *edges = currentContig->GetLeftEdges();
		for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
			if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() )
				continue;

			Contig *otherContig = (*iter)->GetOtherContig( currentContig );
			if( otherContig->IfInAR() ){
				//cout<<"The contig "<<c->GetName()<<" which is connected with in active region: "<<otherContig->GetName()<<endl;
				delete visitedContigs;
				delete possibleContigs;
				return true;
			}

			if( visitedContigs->find( otherContig ) == visitedContigs->end() ){
				visitedContigs->insert( otherContig );
				possibleContigs->push_back( otherContig );
			}
		}

		// check right edges
		edges = currentContig->GetRightEdges();
		for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
			if( !(*iter)->isInSubgraph() || (*iter)->IfUnhappy() )
				continue;

			Contig *otherContig = (*iter)->GetOtherContig( currentContig );
			if( otherContig->IfInAR() ){
				//cout<<"The contig "<<c->GetName()<<" which is connected with in active region: "<<otherContig->GetName()<<endl;
				delete visitedContigs;
				delete possibleContigs;
				return true;
			}

			if( visitedContigs->find( otherContig ) == visitedContigs->end() ){
				visitedContigs->insert( otherContig );
				possibleContigs->push_back( otherContig );
			}
		}
	}

	delete visitedContigs;
	delete possibleContigs;
	return false;
}


// trace back to parent partial scaffold
// change p points to its parent
void opera::TraceBack( PartialScaffold *&p ){
	// recover active region
	Contig *addedContig = p->GetAddedContigInAR();
	if( addedContig != NULL && addedContig->IfInAR() ){
		addedContig->SetIfInAR( false );
		m_activeRegion->pop_back();
	}

	list<Contig*> *removedContigs = p->GetRemovedContigsInAR();
	for( list<Contig*>::reverse_iterator iter = removedContigs->rbegin(); 
		iter != removedContigs->rend(); iter++ ){
			m_activeRegion->push_front( *iter );
			(*iter)->SetIfInAR( true );
	}

	// recover added happy and unhappy dangling edges
	list<PET*> *addedPets = p->GetAddedHappyDE();
	list<PET*>::iterator iter;
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		PET *edge = *iter;
		m_happyDanglingEdges->erase( edge->GetHappyDEIter() );
		edge->SetDE( false );
	}

	addedPets = p->GetAddedUnhappyDE();
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		PET *edge = *iter;
		m_unhappyDanglingEdges->erase( edge->GetUnhappyDEIter() );
		edge->SetDE( false );
		edge->SetHappy();
	}

	// recover removed happy and unhappy dangling edges
	list<PET*> *removedPets = p->GetRemovedHappyDE();
	for( iter = removedPets->begin(); iter != removedPets->end(); iter++ ){
		PET *edge = *iter;
		m_happyDanglingEdges->push_front( edge );
		edge->SetHappyDEIter( m_happyDanglingEdges->begin() );
		edge->SetDE( true );
	}

	removedPets = p->GetRemovedUnhappyDE();
	for( iter = removedPets->begin(); iter != removedPets->end(); iter++ ){
		PET *edge = *iter;
		m_unhappyDanglingEdges->push_front( edge );
		edge->SetUnhappyDEIter( m_unhappyDanglingEdges->begin() );
		edge->SetDE( true );
		edge->SetUnhappy();
	}

	// recover unhappy edges
	addedPets = p->GetAddedUnhappyEdges();
	for( iter = addedPets->begin(); iter != addedPets->end(); iter++ ){
		(*iter)->SetHappy();
	}

	// recover unassigned node sets
	list< pair<Contig*, int> >::iterator nodeIter;
	list< pair<Contig*, int> > *removedNodes = p->GetRemovedUnassignedNodes();
	for( nodeIter = removedNodes->begin(); nodeIter != removedNodes->end(); nodeIter++ ){
		pair<Contig*, int> newPair = *nodeIter;
		m_unassignedNodes->push_front( newPair );
		newPair.first->SetUnassignedNodeIter( newPair.second, m_unassignedNodes->begin() );
	}

	list< pair<Contig*, int> > *addedNodes = p->GetAddedUnassignedNodes();
	for( nodeIter = addedNodes->begin(); nodeIter != addedNodes->end(); nodeIter++ ){
		pair<Contig*, int> newPair = *nodeIter;
		m_unassignedNodes->erase( newPair.first->GetUnassignedNodeIter( newPair.second ) );
		newPair.first->ClearUnassignedNodeIter( newPair.second );
	}

	PartialScaffold *tempSca = p;
	p = p->GetParent();
	delete tempSca;
}

// generate results of scaffold
// save the scaffolds string
// and calculate the distance of each contigs in scaffolds
void opera::GenerateResults( PartialScaffold *p, ScaffoldResult **&results, int num ){
	// calculate gap
	PartialScaffold *tempPS = p;
	int scaffoldID = 0;
	int contigOrderID = 0;
	int numberOfGaps = 0;
	vector<Contig*> *contigInScaffold = new vector<Contig*>;		// save the contigs in one scaffold
	while( tempPS != NULL ){
		//if( tempPS->GetNumOfContigsInAR() == 0 && tempPS->GetNumOfHappyDE() == 0 ){
		if( tempPS->IfEmptyScaffold() ){
			// it is a new scaffold
			//if( tempPS->GetNumOfUnhappyDE() != 0 || tempPS->IfMannualEmptyScaffold() ){
				// calculate the gap sizes
				if( contigInScaffold->size() > 1 )
					CalculateGap( contigInScaffold, numberOfGaps - 1 );
				else{
#ifdef DEBUG
					cout<<"singleton scaffold"<<endl;
#endif
					contigInScaffold->at( 0 )->SetGap( 0 );
				}
				scaffoldID++;
			//}
			
			// start a new scaffold
			if( tempPS->IfMannualEmptyScaffold() )
				tempPS = tempPS->GetParent();

			contigOrderID = 0;
			numberOfGaps = 0;
			contigInScaffold->clear();
		}

		Contig *addedContig = tempPS->GetAddedContigInAR();
		addedContig->SetScaffoldPos( contigOrderID++ );
		addedContig->SetScaffoldID( scaffoldID );
		contigInScaffold->push_back( addedContig );
		numberOfGaps++;
		tempPS = tempPS->GetParent();
	}

	// calculate gap for last scaffold
	if( contigInScaffold->size() > 1 )
		CalculateGap( contigInScaffold, numberOfGaps - 1 );
	else{
#ifdef DEBUG
		cout<<"singleton scaffold"<<endl;
#endif
		contigInScaffold->at( 0 )->SetGap( 0 );
	}

	delete contigInScaffold;

	// generate scaffold
	PartialScaffold *tail = p;
	string scaffoldResult = "";
	double scaffoldLength = 0;
	double scaffoldCov = 0;		// the average coverage of scaffolds
	scaffoldID = 0;
	double endScaffoldPos = p->GetAddedContigInAR()->GetStartPosition()
					+ p->GetAddedContigInAR()->GetLength() - 1;
	// cout<<"start to generate scaffold results"<<endl;
	bool firstContig = true;
	while( p != NULL ){
		// save unhappy edges information
		list<PET*> *unhappyEdges = p->GetAddedUnhappyEdges();
		unhappyEdges->insert( unhappyEdges->begin(), p->GetAddedUnhappyDE()->begin(),
			p->GetAddedUnhappyDE()->end() );
		for( list<PET*>::iterator iter = unhappyEdges->begin(); iter != unhappyEdges->end(); iter++ ){
			m_unhappyEdgeString->push_back( (*iter)->GetOriString() );
		}
		//if( p->GetNumOfContigsInAR() == 0 && p->GetNumOfHappyDE() == 0 ){
		if( p->IfEmptyScaffold() ){
			// it is a new scaffold
			//if( p->GetNumOfUnhappyDE() != 0  || p->IfMannualEmptyScaffold() ){
				// calculate the gap sizes
				// not the end of solution, save previous results
#ifdef SPLIT
				fprintf( logFile, "new scaffold result %d\n", scaffoldID );
				if( scaffoldID >= num )
					fprintf( logFile, "Error: scaffoldID %d > total number of scaffolds (%d)\n", scaffoldID, num );
#endif

				results[ scaffoldID ] = new ScaffoldResult();
				results[ scaffoldID ]->SetScaffoldString( scaffoldResult );
				results[ scaffoldID ]->SetLength( scaffoldLength );
				scaffoldCov = scaffoldCov / scaffoldLength;
				results[ scaffoldID ]->SetCov( scaffoldCov );
				scaffoldID++;
			//}
			
			// start a new scaffold
			scaffoldResult = "";
			scaffoldLength = 0;
			scaffoldCov = 0;
			if( p->IfMannualEmptyScaffold() )
				p = p->GetParent();
			firstContig = true;

			endScaffoldPos = p->GetAddedContigInAR()->GetStartPosition()
					+ p->GetAddedContigInAR()->GetLength() - 1;
		}

		// generate the scaffold string
		Contig *addedContig = p->GetAddedContigInAR();
		if( !firstContig )
			scaffoldResult = "\n" + scaffoldResult;
		else
			firstContig = false;

		scaffoldResult = addedContig->GenScaffoldString( addedContig->GetOri() ) 
			+ itos( addedContig->GetGap() ) + scaffoldResult;
		scaffoldLength += addedContig->GetLength();
		scaffoldCov += addedContig->GetCov() * addedContig->GetLength();

		// calculate the left and right distance of contig
		addedContig->SetScaffoldID( scaffoldID );
		double leftDis = addedContig->GetStartPosition() - 1;
		double rightDis = endScaffoldPos - addedContig->GetStartPosition() - addedContig->GetLength() + 1;
#ifdef CHECK
		if( leftDis < 0 ){
			cout<<addedContig->GetScaffoldString()<<endl;
			cout<<"leftDis < 0: "<<leftDis<<endl<<endl;
		}
		
		if( rightDis < 0 ){
			cout<<addedContig->GetScaffoldString()<<endl;
			cout<<"end scaffold pos: "<<endScaffoldPos<<endl;
			cout<<"rightDis < 0: "<<rightDis<<endl<<endl;
		}
#endif
		addedContig->SetLeftDis( leftDis );
		addedContig->SetRightDis( rightDis );

		p = p->GetParent();
	}

	// save the last results
#ifdef SPLIT
	fprintf( logFile, "new scaffold result %d\n", scaffoldID );
#endif
	if( scaffoldID != num - 1 ){
		cout<<"Error: scaffoldID "<<scaffoldID<<" is not equal to total number of scaffolds "<<num<<endl;
		exit(0);
	}
	results[ scaffoldID ] = new ScaffoldResult();
	results[ scaffoldID ]->SetScaffoldString( scaffoldResult );
	results[ scaffoldID ]->SetLength( scaffoldLength );
	scaffoldCov = scaffoldCov / scaffoldLength;
	results[ scaffoldID ]->SetCov( scaffoldCov );
}

// calculate the gap sizes
void opera::CalculateGap( vector<Contig*> *contigs, int numberOfGaps ){
	// initialize gap sizes
	Matrix<double> G, CE, CI;
	Vector<double> g0, ce0, ci0, x;
	//int n, m, p;
	double sum = 0.0;

	G.resize( numberOfGaps, numberOfGaps );
	g0.resize( numberOfGaps );
	for (int i = 0; i < numberOfGaps; i++){
		for (int j = 0; j < numberOfGaps; j++)
			G[i][j] = 0;
		g0[ i ] = 0;
	}

	CE.resize( numberOfGaps, 0 );
	CI.resize( numberOfGaps, 0 );
	ce0.resize( 0 );
	ci0.resize( 0 );
	x.resize( numberOfGaps );

	// traverse each contig
#ifdef LOG
	//fprintf( logFile, "\nnew scaffolds\n" );
#endif
	set<PET*> *visitedEdges = new set<PET*>;
	for( vector<Contig*>::iterator contigIter = contigs->begin(); contigIter!= contigs->end(); contigIter++ ){
		Contig *contig = *contigIter;
#ifdef LOG
		//fprintf( logFile, "%s\t%d\t%.0f\n", contig->GetName().c_str(), contig->GetOri(), contig->GetLength() );
#endif
		// update parameters
		list<PET*> *edges = new list<PET*>;
		edges->insert( edges->begin(), contig->GetLeftEdges()->begin(), contig->GetLeftEdges()->end() );
		edges->insert( edges->begin(), contig->GetRightEdges()->begin(), contig->GetRightEdges()->end() );

		for( list<PET*>::iterator edgeIter = edges->begin(); edgeIter != edges->end(); edgeIter++ ){
			PET *edge = *edgeIter;
			if( !edge->isInSubgraph() || edge->IfUnhappy() )
				continue;

			// check if visited this edge before
			if( visitedEdges->find( edge ) == visitedEdges->end() )
				visitedEdges->insert( edge );
			else
				continue;

#ifdef SPLIT
			cout<<edge->GetOriString()<<endl;
#endif
			int firstPos = edge->GetStartContig()->GetScaffoldPos();
			int secondPos = edge->GetEndContig()->GetScaffoldPos();
			if( firstPos > secondPos ){
				int temp = secondPos;
				secondPos = firstPos;
				firstPos = temp;
			}

			double sumLength = 0;		// the total contig length in the middle
			for( int id = firstPos + 1; id < secondPos; id++ )
				sumLength += contigs->at( id )->GetLength();

			// update related gap coefficients
			for( int gapID = firstPos; gapID < secondPos; gapID++ ){
				g0[ gapID ] += 2.0 * ( sumLength - (double) edge->GetDis() ) 
					/ ( (double)edge->GetStd() * (double)edge->GetStd() );
				G[ gapID ][ gapID ] += 2.0 / (double)( edge->GetStd() * edge->GetStd() );
				for( int gapID2 = firstPos; gapID2 < secondPos; gapID2++ )
				{
					if( gapID != gapID2 )
						G[ gapID ][ gapID2 ] += 2.0 / (double)( edge->GetStd() * edge->GetStd() );
				}
			}
		}

		delete edges;
	}
	visitedEdges->clear();
	delete visitedEdges;
		
	solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

	// set the gap sizes to contigs
	contigs->at( 0 )->SetGap( 0 );
	for( int id = 1; id < contigs->size(); id++ ){
		//if( x[ id - 1 ] > 0 )
			contigs->at( id )->SetGap( (int)x[ id - 1 ] );
		//else
		//	contigs->at( id )->SetGap( 1 );
	}
}

// output scaffolds
// return -1 if fails
int opera::OutputScaffold( string fileName ){
	return m_graph->OutputScaffolds( fileName );
}

// output unhappy edges
int opera::OutputUnhappyEdges( string fileName ){
	ofstream unhappyWriter( fileName.c_str() );
	if( unhappyWriter == NULL ){
		cout<<"error writing unhappyEdges file!"<<endl;
		return -1;
	}

	string head = "First Contig\tOrientation\tSecond Contig\tOrientation\tDistance\tStandard Deviation\tSize\n";
	unhappyWriter.write( head.c_str(), head.length() );
	for( list<string>::iterator iter = m_unhappyEdgeString->begin(); iter != m_unhappyEdgeString->end();
		iter++ ){
			string line = (*iter) + "\n";
			unhappyWriter.write( line.c_str(), line.length() );
	}

	unhappyWriter.close();
	return 1;
}

// print scaffold result
void opera::PrintScaffoldResult( ScaffoldResult **results, int number ){
	for( int i = 0; i < number; i++ ){
		cout<<"Scaffold "<<i<<endl;
		cout<<results[ i ]->GetScaffoldString()<<endl;
	}
}

// Get statistics of assembly
string opera::GetStats(){
	return m_stats;
}

// sort all scaffold according to length
// return -1 if failed
int opera::SortScaffold(){
	multiset<FinalScaffold*, great_length> *scaffolds = new multiset<FinalScaffold*, great_length>;
	double totalLength = 0;

	// read small contigs
	ReadContigFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "smallContigs", scaffolds, totalLength );
	//cout<<"After read small contig, total Length is "<<totalLength<<endl;
	
	// read repeats
	ReadContigFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "repeatContigs", scaffolds, totalLength );
	//cout<<"After read repeat: "<<totalLength<<endl;

	// read scaffold
	int nonSingleton = 0;
	ReadScaffoldFile( Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "scaffolds.scaf", scaffolds, totalLength, nonSingleton );
	//cout<<"After read scaffold: "<<totalLength<<endl; 

	// calcualte N50
	double N50 = 0;
	multiset<FinalScaffold*, great_length>::iterator iter = scaffolds->begin();
	while( iter != scaffolds->end() ){
		N50 += (*iter)->GetLength();
		if( N50 > 0.5 * totalLength ){
			N50 = (*iter)->GetLength();
			break;
		}
		iter++;
	}
	iter = scaffolds->begin();
	string fileName = Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "statistics";

	// write in another way
	m_stats = "\nScaffold statistics:\n";
	m_stats += "\tN50: ";
	//cout<<endl<<"Scaffold statistics:"<<endl;
	//cout<<"\tN50: ";
	m_stats += PrintNumber( N50 );
	m_stats += "\n";
	//cout<<"\tTotal length: ";
	m_stats += "\tTotal length: ";
	m_stats += PrintNumber( totalLength );
	m_stats += "\n";
	//cout<<endl;
	//cout<<"\tLongest scaffold: ";
	m_stats += "\tLongest scaffold: ";
	m_stats += PrintNumber( (*iter)->GetLength() );
	m_stats += "\n";
	//cout<<endl;

	FILE *contigWriter = fopen( fileName.c_str(), "w" );
	fprintf( contigWriter, "N50: " );
	PrintNumberToFile( contigWriter, N50 );
	fprintf( contigWriter, "\n" );
	fprintf( contigWriter, "Total length of assembly: " );
	PrintNumberToFile( contigWriter, totalLength );
	fprintf( contigWriter, "\n" );
	fprintf( contigWriter, "Number of non-singleton scaffolds: %'i \n", nonSingleton );
	fprintf( contigWriter, "Length of the longest scaffold: " );
	PrintNumberToFile( contigWriter, (*iter)->GetLength() );
	fprintf( contigWriter, "\n" );
	fclose( contigWriter );

	// read original fasta file
	map<string, string> *contigs = new map<string, string>;
	ReadContigFasta( Configure::CONTIG_FILE, contigs );

	// output results
	fileName = Configure::OUTPUT_FOLDER + Configure::DIRECTORY_SIGN + "scaffoldSeq.fasta";
	ofstream resultWriter( fileName.c_str() );

	if( resultWriter == NULL ){
		cout<<"error writing "<<fileName<<endl;
		return -1;
	}
	iter = scaffolds->begin();
	vector<string> *scafVec = new vector<string>;
	vector<string> *scafLine = new vector<string>;
	while( iter != scaffolds->end() ){
		string resultString = ">" + (*iter)->GetName();
		resultWriter.write( resultString.c_str(), resultString.length() );

		if( !(*iter)->IsScaffold() ){
			// small contigs and repeats
			char tempName[500];
			sprintf( tempName, "\tlength: %.0f\tcov: %.1f\n", (*iter)->GetLength(), (*iter)->GetCov() );
			string seq( tempName );
			resultWriter.write( seq.c_str(), seq.length() );

			seq = (*contigs)[ (*iter)->GetName() ];
			resultWriter.write( seq.c_str(), seq.length() );
		}
		else{
			// scaffolds
			resultString = "\n";
			resultWriter.write( resultString.c_str(), resultString.length() );

			string scafString = (*iter)->GetScaffold();
			Split( scafString, "\n", scafVec );
			string content = "";
			for( int i = 0; i < scafVec->size(); i++ ){
				scafLine->clear();
				Split( (*scafVec)[ i ], "\t", scafLine );
				content = (*contigs)[ (*scafLine)[ 0 ] ];
				// output each contig
				if( (*scafLine)[ 1 ] == "BE" ){
					// output "+"
					resultWriter.write( content.c_str(), (*contigs)[ (*scafLine)[ 0 ] ].length() );
				}
				else{
					// output "-"
					string reverseSeq = "";
					for ( string::reverse_iterator reverseIter = content.rbegin(); reverseIter != content.rend(); reverseIter++ )
					{
						if( *reverseIter == 'A' )
							reverseSeq.append( "T" );
						else if( *reverseIter == 'C' )
							reverseSeq.append( "G" );
						else if( *reverseIter == 'G' )
							reverseSeq.append( "C" );
						else if( *reverseIter == 'T' )
							reverseSeq.append( "A" );
					}
					resultWriter.write( reverseSeq.c_str(), reverseSeq.size() );
				}
				int gapSize = atoi( (*scafLine)[ 3 ].c_str() );
				// add the gaps except for the last contig
				if( gapSize < 10 && i != scafVec->size() -1 )
					gapSize = 10;
				string gap = "";
				for( int num = 0; num < gapSize; num++ )
					gap.append( "N" );
				resultWriter.write( gap.c_str(), gap.size() );
			}
			scafVec->clear();
		}
		string newLine = "\n";
		resultWriter.write( newLine.c_str(), newLine.length() );
		iter++;
	}

	resultWriter.close();

	// clear
	iter = scaffolds->begin();
	while( iter != scaffolds->end() ){
		//cout<<(*iter)->GetName()<<"\t"<<(*iter)->GetLength()<<endl;
		delete *iter;
		scaffolds->erase( iter++ );
	}
	delete scaffolds;

	contigs->clear();
	delete contigs;
	scafLine->clear();
	delete scafLine;
	scafVec->clear();
	delete scafVec;
	return 0;
}

// read result contig files
int opera::ReadContigFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"error reading contig file: "<<fileName<<endl;
		return -1;
	}

	string line;
	vector<string> *contents = new vector<string>;
	getline( contigReader, line );		// remove the first title line
	while( getline( contigReader, line ) ){
		Split( line, "\t", contents );
		FinalScaffold *newFS = new FinalScaffold( (*contents)[ 0 ], atoi( (*contents)[ 1 ].c_str() ), "", atoi( (*contents)[ 2 ].c_str() ) );
		totalLength += newFS->GetLength();
		scaffolds->insert( newFS );
	}

	delete contents;
	contigReader.close();
	return 0;
}

// read scaffolds file
int opera::ReadScaffoldFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength, int &nonSingleton )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"error reading scaffold file: "<<fileName<<endl;
		return -1;
	}

	string line;
	vector<string> *contents = new vector<string>;
	bool isFirst = true;
	int length = 0;
	string name = "";
	string scaffold = "";
	int numOfContig = 0;
	while( getline( contigReader, line ) ){
		if( line.at( 0 ) == '>' ){
			// save previous scaffold
			if( isFirst )
				isFirst = false;
			else{
				FinalScaffold *newFS = new FinalScaffold( name, length, scaffold );
				totalLength += length;
				scaffolds->insert( newFS );
				if( numOfContig > 1 )
					nonSingleton++;
			}

			// start a new scaffold
			length = 0;
			name = line.substr( 1 );
			scaffold = "";
			numOfContig = 0;
		}
		else{
			Split( line, "\t", contents );
			scaffold += line + "\n";
			length += atoi( (*contents)[ 2 ].c_str() );
			numOfContig++;
		}
	}
	// save the last scaffold
	FinalScaffold *newFS = new FinalScaffold( name, length, scaffold );
	totalLength += length;
	scaffolds->insert( newFS );
	if( numOfContig > 1 )
		nonSingleton++;

	contigReader.close();
	delete contents;
	return 0;
}

// read contig fasta file
int opera::ReadContigFasta( string fileName, map<string, string> *contigs )
{
	ifstream contigReader( fileName.c_str() );

	if( contigReader == NULL )
	{
		cout<<"error reading scaffold file: "<<fileName<<endl;
		return -1;
	}

	string line;
	string seq = "";
	string name = "";
	vector<string> *column = new vector<string>;
	while( getline( contigReader, line ) ){
		if( line.length() == 0 )
			continue;
		if( line.at( 0 ) == '>' ){
			// save previous contig
			if( seq.length() > 0 ){
				contigs->insert( pair<string, string>( name, seq ) );
			}

			// start new contig
			Split( line, " ", column );
			name = (*column)[ 0 ].substr( 1, (*column)[ 0 ].length() - 1 );
			seq = "";
		}
		else{
			seq.append( line );
		}
	}
	column->clear();
	delete column;

	// save last contig
	contigs->insert( pair<string, string>( name, seq ) );
	return 0;
}

// check contig file format: velvet, soap or normal fasta file
int opera::CheckContigFormat()
{
	ifstream contigReader( Configure::CONTIG_FILE.c_str() );

	if( contigReader == NULL )
	{
		cout<<"error reading contig file"<<endl;
		return -1;
	}

	string line;
	getline( contigReader, line );

	if( line.find( "NODE_" ) != string::npos && line.find( "_length_" ) != string::npos && line.find( "_cov_" ) != string::npos ){
	//if( line.find( "contig_" ) != string::npos && line.find( "_length_" ) != string::npos && line.find( "_cov_" ) != string::npos ){
		// velvet format
		Configure::FILE_TYPE = VELVET;
		Configure::FILTER_REPEAT = true;
	}
	else{
		vector<string> *column = new vector<string>;
		Split( line, " ", column );
		if( ( column->size() == 2 ) && IsNumber( (*column)[ 1 ] ) ){  //">scaffold" ){
			// soap scaffold format
			//cout<<"soap format\n";
			Configure::FILE_TYPE = SOAP;
			Configure::FILTER_REPEAT = true;
		}
		else{
			if( line.find( "_tip_" ) != string::npos && line.find( " length " ) != string::npos && line.find( " cvg_" ) != string::npos ){
				// soap contig format
				//cout<<"soap contig format\n";
				Configure::FILE_TYPE = SOAP_CONTIG;
				Configure::FILTER_REPEAT = true;
			}
			else{
				// normal fasta format
				Configure::FILE_TYPE = NORMAL;
				Configure::FILTER_REPEAT = false;
				cout<<"Warning: The contig file is a normal fasta file. Please make sure all repeat contigs have been removed.\n";
			}
		}
		delete column;
	}

	return 1;
}

void opera::OpenLogFile()
{
	logFile = fopen( ( Configure::OUTPUT_FOLDER + "log" ).c_str(), "w" );
}

// move the edges in first library into graph
void opera::MoveEdgesToGraph( PetLibrary *lib ){
	m_graph->ClearEdge();

	// add edge to graph
	list<PET*> *edges = lib->GetEdges();
	list<PET*>::iterator edgeIter = edges->begin();
	//cout<<edges->size()<<" edges in this library\n";
	int tempNum = 0;
	while( edgeIter != edges->end() ){
		PET *edge = *edgeIter;

		if( edge->GetDis() + Configure::STD_TIMES * edge->GetStd() >= -Configure::KMER ){
			// a valid edge, save to graph
			m_graph->AddEdge( edge );
			tempNum++;
		}

		// remove edge from unused edges list in contigs
		edge->GetStartContig()->RemoveEdgeMultiLib( edge, START );
		edge->GetEndContig()->RemoveEdgeMultiLib( edge, END );

		edgeIter = edges->erase( edgeIter );

		if( edge->GetDis() + Configure::STD_TIMES * edge->GetStd() < -Configure::KMER )
			delete edge;
	}
}

// the pipeline of optimal scaffolding in opera
int main(int argc, char *argv[] )
{
	time_t start, end;
	struct timeval t_start,t_end;
	struct timeval t_startTemp,t_endTemp;
	time( &start );
#ifdef TIME
	gettimeofday( &t_start, NULL );
#endif
	
	// check parameters
	string configFileName = "";
	if( argc == 1 || argc == 3 || argc > 4 || argv[ 1 ] == "-h" ){
		cout<<"Usage:\n";
		cout<<"\topera <config file>\n";
		cout<<"\t\tOR\n";
		cout<<"\topera <contig file> <mapping file,(mapping file 2)> <output folder>\n";
		cout<<"\nNOTE: Please refer to test_dataset/multiLib.config for detailed settings.\n";
		return -1;
	}
	else if( argc == 2 ){
		configFileName = string( argv[ 1 ] );
		
		// step 0: read config file
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 1: Reading configuration file ..."<<endl;
		flush(cout);
		configureReader myConfig;
		if( myConfig.ReadConfigFile( configFileName ) == -1 ){
			cout<<"ERROR: The format of configuration file is not correct"<<endl; 
			return -1;
		}
#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		flush(cout);
#endif
	}
	else if( argc == 4 ){
		// get contig file name, mapping file name and output folder
		// step 0: read config file
#ifdef TIME
		gettimeofday( &t_startTemp, NULL );
#endif
		cout<<"Step 1: Setting parameters ..."<<endl;
		flush(cout);
		Configure::CONTIG_FILE = argv[ 1 ];
		Configure::MAP_FILE = argv[ 2 ];
			
		Configure::OUTPUT_FOLDER = argv[ 3 ];
		mkdir (Configure::OUTPUT_FOLDER.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
		if( Configure::OUTPUT_FOLDER.substr( Configure::OUTPUT_FOLDER.length() - 1, 1 ) != Configure::DIRECTORY_SIGN )
			Configure::OUTPUT_FOLDER += Configure::DIRECTORY_SIGN;
#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		flush(cout);
#endif
	}

	opera *m_opera = new opera();

	// check if all files exist
	if( !m_opera->CheckFileExist() ){
		return -1;
	}

	// step 0.5: check if the scaffold name will have conflicts
	m_opera->OpenLogFile();
	if( m_opera->CheckNameConfliction() ){
		// there is confliction
		cout<<"ERROR: There might be conflictions of scaffold name. Please specify another name using parameter: scaffold_name in configuraton file"<<endl;
		return -1;
	}

	// check contig file format
	if( m_opera->CheckContigFormat() == -1 ){
		return -1;
	}

	// step 1: convert contig file
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 2: Reading contig file ..."<<endl;
	flush(cout);
	ContigConverter myContigConverter;
	if( myContigConverter.ConvertContigFile( Configure::CONTIG_FILE, m_opera->m_graph ) == -1 ){
		cout<<"ERROR: Converting contig file error!"<<endl;
		return -1;
	}
#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	// step 2: convert mapping file
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 3: Reading mapping file ..."<<endl;
	flush(cout);
	MapConverter myMapConverter( m_opera->m_graph );
	/*if( myMapConverter.Analyze( Configure::MAP_FILE ) == -1 ){
		cout<<"ERROR: Converting mapping file error!"<<endl;
		return -1;
	}*/
	if( myMapConverter.AnalyzeMultiLib( Configure::MAP_FILE, m_opera->m_libraries ) == -1 ){
		cout<<"ERROR: Converting mapping file error!"<<endl;
		return -1;
	}

#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	// step 3: bundle
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 4: Bundling paired-end reads ..."<<endl;
	flush(cout);
	/*if( Configure::MAP_TYPE != OPERA ){
		if( myMapConverter.Bundle() == -1 ){
			cout<<"ERROR: Bundling error!"<<endl;
			return -1;
		}
	}*/

	if( Configure::MAP_TYPE != OPERA ){
		if( myMapConverter.BundleMultiLib( m_opera->m_libraries ) == -1 ){
			cout<<"ERROR: Bundling error!"<<endl;
			return -1;
		}
	}
#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	flush(cout);
#endif

	// traverse each library
	list<PetLibrary*>::iterator libIter = m_opera->m_libraries->begin();
	int libNum = 0;
	
	// step 4: scaffolder
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 5: Scaffolding ..."<<endl;
	flush(cout);

	while( libIter != m_opera->m_libraries->end() ){
		libNum++;
		cout<<"Dealing with library: "<<(*libIter)->GetFileName()<<endl;

		PetLibrary *currentLib = *libIter;
		Configure::LIB_MEAN = currentLib->GetMean();
		Configure::LIB_STD = currentLib->GetStd();

		// initialize opera for a new run
		m_opera->Initialize();

		// move the edges to graph
		m_opera->MoveEdgesToGraph( currentLib );

		// delete current library
		delete currentLib;

		// split the contigs according to their types
		m_opera->m_graph->InitializeContig( libNum );
		m_opera->m_graph->GenerateBorderAndScaffold();

		if( m_opera->StartScaffold() == -1 ){
			cout<<"ERROR: Scaffolding error!"<<endl;
			return -1;
		}
#ifdef TIME
		gettimeofday( &t_endTemp, NULL );
		cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
		flush(cout);
#endif

		// remove first library
		libIter = m_opera->m_libraries->erase( libIter );
	}

	// step 6: finalize results
#ifdef TIME
	gettimeofday( &t_startTemp, NULL );
#endif
	cout<<"Step 6: Outputing results ..."<<endl;
	flush(cout);
	if( m_opera->OutputScaffold( Configure::OUTPUT_FOLDER + "scaffolds.scaf" ) == -1 ){
		cout<<"ERROR: Output scaffolds result error!"<<endl;
		return -1;
	}
#ifdef TIME
	//cout<<"finish output .scaf file. Now generating sequence file..."<<endl;
	//flush(cout);
#endif

	m_opera->SortScaffold();

	if( m_opera->OutputUnhappyEdges( Configure::OUTPUT_FOLDER + "unhappyEdges" ) == -1 ){
		cout<<"ERROR: Output unhappy edges error!"<<endl;
		return -1;
	}
#ifdef TIME
	gettimeofday( &t_endTemp, NULL );
	cout<<"Time Taken: "<<(t_endTemp.tv_sec - t_startTemp.tv_sec) + (t_endTemp.tv_usec - t_startTemp.tv_usec)/1000000.0<<" seconds"<<endl;
	gettimeofday( &t_end, NULL);
#endif

	// print the statistics
	bool printTime = false;
	cout<<m_opera->GetStats()<<endl;

	time( &end );
	cout<<"Total running time is ";
	int days = 0, hours = 0, minutes = 0, seconds = 0;
	double dif = difftime( end, start );
	if( dif > 60*60*24 ){
		days = (int) (dif / 60*60*24);
		dif -= days * 60*60*24;
	}
	if( dif > 60*60 ){
		hours = (int) (dif / 60*60);
		dif -= hours * 60 * 60;
	}
	if( dif > 60 ){
		minutes = (int)(dif / 60);
		dif -= minutes * 60;
	}

	seconds = (int)dif;
	if( days != 0 ){
		printTime = true;
		cout<<days<<" days ";
	}
	if( hours != 0 ){
		printTime = true;
		cout<<hours<<" hours ";
	}
	if( minutes != 0 ){
		printTime = true;
		cout<<minutes<<" minutes ";
	}
	if( seconds != 0 ){
		printTime = true;
		cout<<seconds<<" seconds ";
	}
	if( !printTime ){
		cout<<"0 seconds";
	}
	cout<<endl;

	delete m_opera;
	
cout<<"The results are in "<<Configure::OUTPUT_FOLDER<<endl;

#ifdef DEBUG
	system("pause");
#endif
}

