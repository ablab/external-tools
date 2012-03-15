#pragma once

#include <string>
#include <list>
#include <vector>

class PET;
#include "PET.h"
#include "CommonFunction.h"

using namespace std;

class Contig
{
public:
	Contig(void);
	Contig( string name, double length, double cov );
	~Contig(void);

	// Attributes
private:
	string m_name;		// contig name
	double m_length;		// contig length
	double m_cov;			// contig coverage
	int m_id;			// contig ID
	list<PET*> *m_leftEdges;		// the left edges, suppose ori is +
	list<PET*> *m_rightEdges;		// the right edges, suppose ori is +
	list<Contig*>::iterator m_listPos;		// the list iterator of current contig
	list<Contig*>::iterator m_borderListPos;		// the border list iterator of current contig
	bool m_ifInSubgraph;		// the flag to record if this contig is in current subgraph
	double m_leftDistance;			// the left distance to the scaffold, suppose ori is +
	double m_rightDistance;		// the right distance to the scaffold, suppose ori is +
	string m_scaffold;			// the scaffold content of this contig
	int m_scaffoldID;
	bool m_isBorderContig;		// flag to record if it is a border contig
	int m_extentionDirection;		// the direction during extension

	// scaffold related
	int m_ori;			// orientation in scaffold
	double m_startPosition;			// the startPosition in scaffold without gap
	double m_startPositionWithGap;		// the start position in scaffold with gap
	list< pair<Contig*, int> >::iterator m_unassignedNodeIter[2];		// iterator in unassigned node set, both plus and minus
	bool m_isUnassignedNode[2];			// record if it is unassigned node
	bool m_isInAR;				// recod if it is in active region
	list< pair<Contig*, int> >::iterator m_endOfList;	
	//*****unassigned nodes*****
	int m_unassignedOri;			// the untried orientation in unassigned nodes
	list< pair<Contig*, int> > *m_emptyList;			// the empty list, used for get end() iterator
	//*****end of unassigned nodes******
	//***********gap sizes related
	int m_scaffoldPos;			// the order in scaffolds
	int m_gapSize;				// the gap size after this scaffold
	//***********

	// heuristic parameters
	int m_step;			// the steps to find this contig

	// multiple libraries related
	list<PET*> *m_leftEdgesMultiLib;		// the left edges of all libraries, suppose ori is +
	list<PET*> *m_rightEdgesMultiLib;		// the right edges of all libraries, suppose ori is +

	// Methods
public:
	void SetName( string name );
	string GetName();
	void SetLength( double length );
	double GetLength();
	void SetCov( double cov );
	double GetCov();
	void SetID( int id );
	int GetID();
	void AddEdge( PET *p );			// add a pet cluster of this contig
	void SetListPos( list<Contig*>::iterator iter );
	void SetBorderListPos( list<Contig*>::iterator iter );
	bool IsSingleton();			// check if current contig is a singleton
	int GetScaffoldID();
	void SetScaffoldID( int id );
	list<PET*> * GetLeftEdges();
	list<PET*> * GetRightEdges();
	list<Contig*>::iterator GetListPos();
	list<Contig*>::iterator GetBorderListPos();
	bool HasEdge();			// check if current contig has any edge
	bool HasLeftEdge();
	bool HasRightEdge();
	bool IsBorderContig();			// check if current contig is a border contig
	double GetRightDis();
	double GetLeftDis();
	void SetLeftDis( double d );
	void SetRightDis( double d );
	void SetExtensionOri( int ori );
	int GetExtensionOri();
	void SetInSubgraph();
	void SetNotInSubgraph();
	bool isInSubgraph();
	void SetOri( int ori );
	int GetOri();
	void SetStartPosition( double pos, int gap );
	double GetStartPosition();
	double GetStartPositionWithGap();
	list< pair<Contig*, int> >::iterator GetUnassignedNodeIter( int type );
	void SetUnassignedNodeIter( int type, list< pair<Contig*, int> >::iterator iter );
	bool IsUnassignedNode( int type );
	void ClearUnassignedNodeIter( int type );			// set iter to end()
	void SetIfInAR( bool ar );		// set if it is in active region
	bool IfInAR();				// check if it is in active region
	bool IfHasHappyDE( int ori );		// check if it has happy dangling edges in right direction
	void SetEndOfList( list< pair<Contig*, int> >::iterator iter );
	// check the number of valid left Edge for certain orientation
	int GetNumOfValidLeftEdges( int ori );
	// generate scaffold String with certain orientation
	string GenScaffoldString( int ori );
	// set scaffold string
	void SetScaffoldString( string sca );
	// get original contig orientation in scaffold
	int GetContigOriInScaffold( string contigName );
	// reverse left and right edges
	void ReverseEdges();
	// get orientation of complex contig
	string GetOriOfFirstContig( int ori );
	string GetNameOfFirstContig();
	bool CheckOldContigOri( Contig *c );
	string GetScaffoldString();

	// gap related
	void SetScaffoldPos( int pos );
	int GetScaffoldPos();
	void SetGap( int gap );
	int GetGap();

	// heuristic related
	void SetStep( int s );
	int GetStep();

	// multiple library related
	void AddEdgeMultiLib( PET *p );			// add a pet cluster of this contig, from all libraries
	void RemoveEdgeMultiLib( PET *p, int contigPos );		// remove a pet cluster of this contig, from all libraries
	void Initialize();						// initialize the attributes
	bool HasMultiLibEdge();					// check if this contig has multiple libraries edge
	list<PET*>* GetLeftEdgeMultiLib();		// get the left edge of multiple libraries
	list<PET*>* GetRightEdgeMultiLib();		// get the right edge of multiple libraries

private:
	string toScaffoldString();
	inline void DeleteEdges( list<PET*> *edges );
};
