#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include "Configure.h"
#include "Contig.h"
#include "Graph.h"
#include "CommonFunction.h"

using namespace std;

class ContigConverter
{
public:
	ContigConverter(void);
	~ContigConverter(void);

	// Methods
public:
	// Convert contig file
	int ConvertContigFile( string fileName, Graph *graph );

private:
	// analyze velvet file
	int AnalyzeVelvet( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze soapdenovo file
	int AnalyzeSoapDenovo( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze soapdenovo file
	int AnalyzeSoapDenovoContig( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze normal fasta file
	int AnalyzeNormalFasta( ifstream *contigReader, vector<Contig*> *contigs );
	// filter repeat
	void FilterRepeat( vector<Contig*> *contigs, Graph *graph );
	// analyze opera's contig file
	int AnalyzeStatistic( ifstream *contigReader, vector<Contig*> *contigs );
	// remove small contigs only
	void RemoveSmallContigs( vector<Contig*> *contigs, Graph *graph );
	// print list of contigs
	int PrintContigs( list<Contig*> *contigs, string fileName );
	inline void DeleteContigs( list<Contig*> *contigs );

	// Attributes
private:
	// read contig file into a vector
	vector<Contig*> *myContigs;
	list<Contig*> *m_repeatContigs;
	list<Contig*> *m_smallContigs;
};
