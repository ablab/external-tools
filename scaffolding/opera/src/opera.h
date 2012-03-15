#pragma once

#include <list>
#include <set>
#include <map>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>

#include "configureReader.h"
#include "ContigConverter.h"
#include "MapConverter.h"
#include "Configure.h"
#include "Graph.h"
#include "StartPoint.h"
#include "CommonFunction.h"
#include "PartialScaffold.h"
#include "QuadProg.h"
#include "Tree.h"
#include "FinalScaffold.h"
#include "PetLibrary.h"

using namespace std;
using namespace QuadProgPP;

struct great_length {
  bool operator() (const FinalScaffold *s1, const FinalScaffold *s2) const
  {return (s1->GetLength() > s2->GetLength());}
};

class opera
{
public:
	opera(void);
	~opera(void);

	// Attributes
public:
	Graph *m_graph;
	list<PetLibrary*> *m_libraries;		// all the libraries in the data

private:
	list<Contig*> *m_activeRegion;			// the global active region
	list<PET*> *m_happyDanglingEdges;		// the global happy dangling edges
	list<PET*> *m_unhappyDanglingEdges;		// the global unhappy dangling edges
	list< pair<Contig*, int> > *m_unassignedNodes;		// the global unassigned contigs
	list<string> *m_unhappyEdgeString;			// the string of all unhappy edges
	// initialize the tree saving all visited scaffolds
	Tree *m_visitedTree;
	FILE *logFile;						// the log file to record the output of program
	string m_stats;						// the statistics of assembly
	

	// Methods
public:
	// general pipeline of scaffolding
	int StartScaffold();
	// scaffolding
	int Scaffolding( PartialScaffold *&currentScaffold, int maxOfUnhappyEdges );
	// output scaffolds
	int OutputScaffold( string fileName );
	// output unhappy edges
	int OutputUnhappyEdges( string fileName );
	// sort all scaffold according to length
	int SortScaffold();
	// check if the names have confliction
	bool CheckNameConfliction();
	// check contig file format: velvet, soap or normal fasta file
	int CheckContigFormat();
	void OpenLogFile();
	// Get statistics of assembly
	string GetStats();
	// check if all files exist
	bool CheckFileExist();

	// multiple libraries related methods
	// move the edges in first library into graph
	void MoveEdgesToGraph( PetLibrary *lib );
	// initialize for another run
	void Initialize();

private:
	// select the starting contig
	void FindStartingContig( list<Contig*> *subgraph, int numOfContig, int numOfBorderContig,
		Contig *&startContig, int &startContigOri );
	// create a scaffold containing a start contig
	PartialScaffold* CreatScaffold( Contig *c, int ori );
	// add happy dangling edges to first scaffolds
	void AddHappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p );
	// add unhappy dangling edges to first scaffolds
	void AddUnhappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p, Contig *c );
	// generate active region string
	void GenARString( PartialScaffold *p );
	// generate unhappy dangling edges string
	void GenUDEString( PartialScaffold *p );
	// generate unassigned contig set
	void GenUnassigndNodes( PartialScaffold *p );
	// traverse edge and find unassigned node
	void TraverseEdges( Contig *c, list<PET*> *edges, list< pair<Contig*, int> > *&possibleContig,
		 set<Contig*> *&unassignedNodes, set<PET*> *&visitedEdges );
	// add one contig to partial scaffold
	bool AddContigToScaffold( PartialScaffold *&p, int maxUnhappyEdge, Contig *newContig, int newContigOri );
	// check edges of newly added contig
	void CheckEdgesOfNewContig( PartialScaffold *p, Contig *c, list<PET*> *edges, int ori );
	// check incident edge happiness
	bool CheckDanglingEdgeHappiness( PET *edge, Contig *newContig, int ori );
	// check if it is a left edge of current contig with certain orientation
	inline bool IsLeftEdge( PET *edge, Contig *newContig, int ori );
	// check connectivity and update scaffold active region & dangling edges
	bool CheckScaffold( PartialScaffold *p );
	// check if contig c could connect to active region
	bool CheckActiveRegionConnectivity( Contig *c );
	// trace back to parent partial scaffold
	void TraceBack( PartialScaffold *&p );
	// find start contig of new scaffold
	StartPoint* FindStartOfNewScaffold();
	// add connected contigs to corresponding list
	inline void AddToList( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, Contig *c, int ori );
	// check edges to find start point
	inline void CheckEdgesForStartPoint( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, 
		Contig *c, list<PET*> *edges, set<Contig*> *visitedContigs, list<Contig*> *possibleContigs );
	// Break scaffold
	void BreakScaffold( PartialScaffold *&p, int maxUnhappyEdges );
	// generate results of scaffold
	void GenerateResults( PartialScaffold *p, ScaffoldResult **&results, int num );
	// print scaffold result
	void PrintScaffoldResult( ScaffoldResult **results, int number );
	// clear scaffold variables
	inline void Clear();
	// delete all scaffold
	inline void DeleteAllScaffolds( PartialScaffold *p, PartialScaffold *head );
	// calculate gap sizes
	void CalculateGap( vector<Contig*> *contigs, int numberOfGaps );

		// read result contig files
	int ReadContigFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength );
	// read scaffolds file
	int ReadScaffoldFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength, int &nonSingleton );
	// read contig fasta file
	int ReadContigFasta( string fileName, map<string, string> *contigs );
};

bool SortStartPoints( const StartPoint *s1, const StartPoint *s2 );
