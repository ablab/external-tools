#pragma once

#include <string>
#include <list>
#include <map>
#include <set>
class PET;
#include "PET.h"
#include "SinglePet.h"

using namespace std;

struct lessDistance {
  bool operator() (const SinglePet* p1, const SinglePet* p2) const
  {return (p1->GetDistance() < p2->GetDistance());}
};

class PetLibrary
{
public:
	PetLibrary( string name );
	virtual ~PetLibrary(void);

	// Attributes
private:
	int m_mean;		// the mean of this library
	int m_std;		// the standard deviation of this library
	int m_ori;		// orientation of paired reads
	string m_fileName;			// the name of mapping file
	string m_fileNameWithoutPath;		// the pure file name without path and extension
	// the map of single PET, first element is the contig pairs, second element is the set of all edges
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*> *m_singlePetsMap;
	list<PET*> *m_edges;			// the edges of this library

	// Methods
public:	
	void SetMean( int mean );
	void SetStd( int std );
	void SetOri( int ori );
	int GetOri();	
	int GetMean();
	int GetStd();
	void SetNameWithoutPath( string name );
	string GetNameWithoutPath();
	// insert a contig pair into map
	bool InsertPair( int c1, int c2, SinglePet *pet );
	// get the pair map
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*>* GetSinglePetsMap();
	// add a cluster
	void AddCluster( PET *pet );
	// get the pet clusers
	list<PET*>* GetEdges();
	// get file name
	string GetFileName();
};

