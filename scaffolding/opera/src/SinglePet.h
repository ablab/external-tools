#pragma once

// save the paired-end reads

#include <string>
#include <iostream>
#include "Configure.h"

using namespace std;

class SinglePet
{
public:
	SinglePet(void);
	SinglePet( int id, int startContigID, string startContigOri, int endContigID, string endContigOri,
		int dis, int std );
	~SinglePet(void);

	// Methods
public:
	friend ostream& operator<<( ostream& out, const SinglePet& p );

	int GetStartID();
	int GetStartOri();
	int GetEndID();
	int GetEndOri();
	int GetDistance() const;
	int GetStd();

	// Attributes
private:
	int m_id;
	int m_startContigID;
	int m_startContigOri;
	int m_endContigID;
	int m_endContigOri;
	int m_distance;
	int m_std;
};
