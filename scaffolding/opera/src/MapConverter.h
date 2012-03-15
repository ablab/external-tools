#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <list>
#include <set>
#include <algorithm>
#include <map>
#include <time.h>
#include <sys/time.h>
#include "Contig.h"
#include "Configure.h"
#include "Graph.h"
#include "SinglePet.h"
#include "PET.h"
#include "CommonFunction.h"
#include "PetLibrary.h"

using namespace std;

struct less_distance {
  bool operator() (const SinglePet* p1, const SinglePet* p2) const
  {return (p1->GetDistance() < p2->GetDistance());}
};

class MapConverter
{
public:
	MapConverter(void);
	MapConverter( Graph *graph );
	~MapConverter(void);

	// Methods
public:
	// analyze mapping file or opera edge file
	int Analyze( string fileName );
	// start bundle
	int Bundle();

	// for multiple libraries
	// analyze mapping file or opera edge file
	int AnalyzeMultiLib( string fileName, list<PetLibrary*> *libs  );
	// start bundle
	int BundleMultiLib( list<PetLibrary*> *libs );

private:
			// mapping file converter related functions
	// analyze bowtie file
	int AnalyzeBowtie( string fileName );
	// analyze opera edge file
	int AnalyzeOpera( string fileName );
	// calculate library mean and std
	int CalculateLibDis( string fileName );
	// check if two lines represent a pair of reads
	bool IsPair( string firstRead, string secondRead );
	// calculate the distance of two reads in the same contig
	int CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn );
	// convert bowtie format into opera format
	int ConvertBowtieFile( string fileName );
	// convert bowtie format into opera format
	int ConvertBowtieFile();
	// convert bowtie format
	string ConvertBowtieFormat( vector<string> *firstColumn, vector<string> *secondColumn, int id );
	
			// bundle related functions
	// bundle a certain cluster
	string BundleCluster( multiset<SinglePet*, less_distance> *group );

	// for multiple libraries
	// analyze bowtie file
	int AnalyzeBowtieMultiLib( string fileName, list<PetLibrary*> *libs );
	// calculate library mean and std
	int CalculateLibDisMultiLib( string fileName, list<PetLibrary*> *libs );
	// convert bowtie format into opera format
	int ConvertBowtieFileMultiLib( string fileName, list<PetLibrary*> *libs );
	// convert bowtie format
	string ConvertBowtieFormatMultiLib( vector<string> *firstColumn, vector<string> *secondColumn, int id, list<PetLibrary*> *libs );
	// bundle a certain cluster
	string BundleClusterMultiLib( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib );

	// Attributes
private:
	// graph
	Graph *m_graph;
	// number of single PET
	int m_numOfPet;
	// the array of single PET
	list<SinglePet*> *m_singlePets;
	// the map of single PET, first element is the contig pairs, second element is the set of all edges
	map<pair<int, int>, multiset<SinglePet*, less_distance>*> *m_singlePetsMap;
	// list of all paired reads mapped on different scaffold
	list<string> *m_pairMap;
	double m_totalTime;
	double m_pairTime;

	// for multiple libraries
	string m_libString;
};

