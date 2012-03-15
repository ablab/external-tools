#include "PetLibrary.h"


PetLibrary::PetLibrary(string name)
{
	m_fileName = name;
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, lessDistance>*>;
	m_edges = new list<PET*>;
}


PetLibrary::~PetLibrary(void)
{
	delete m_singlePetsMap;
	delete m_edges;
}

void PetLibrary::SetMean( int mean ){
	m_mean = mean;
}
	
void PetLibrary::SetStd( int std ){
	m_std = std;
}

void PetLibrary::SetOri( int ori ){
	m_ori = ori;
}

int PetLibrary::GetOri(){
	return m_ori;
}

int PetLibrary::GetMean(){
	return m_mean;
}

int PetLibrary::GetStd(){
	return m_std;
}

// check if a contig pair exist in map
bool PetLibrary::InsertPair( int c1, int c2, SinglePet *pet ){
	pair<int, int> newPair( c1, c2 );
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*>::iterator mapIter = m_singlePetsMap->find( newPair );
	if( mapIter == m_singlePetsMap->end() ){
		// new pair, create a new set
		multiset<SinglePet*, lessDistance> *newSet = new multiset<SinglePet*, lessDistance>;
		newSet->insert( pet );
		m_singlePetsMap->insert( pair<pair<int, int>, multiset<SinglePet*, lessDistance>*>( newPair, newSet ) );
	}
	else{
		// insert the single pet to previous set
		mapIter->second->insert( pet );
	}
}

void PetLibrary::SetNameWithoutPath( string name ){
	m_fileNameWithoutPath = name;
}

string PetLibrary::GetNameWithoutPath(){
	return m_fileNameWithoutPath;
}

// get the pair map
map<pair<int, int>, multiset<SinglePet*, lessDistance>*>* PetLibrary::GetSinglePetsMap(){
	return m_singlePetsMap;
}

// add a cluster
void PetLibrary::AddCluster( PET *pet ){
	m_edges->push_front( pet );
	pet->AddIterInLib( this, m_edges->begin() );
}

// get the pet clusers
list<PET*>* PetLibrary::GetEdges(){
	return m_edges;
}

// get file name
string PetLibrary::GetFileName()
{
	return this->m_fileName;
}