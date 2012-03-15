#include "SinglePet.h"

SinglePet::SinglePet(void)
{
}

SinglePet::~SinglePet(void)
{
}

int SinglePet::GetStartID(){
	return m_startContigID;
}

int SinglePet::GetStartOri(){
	return m_startContigOri;
}

int SinglePet::GetEndID(){
	return m_endContigID;
}

int SinglePet::GetEndOri(){
	return m_endContigOri;
}

int SinglePet::GetDistance() const{
	return m_distance;
}

int SinglePet::GetStd(){
	return m_std;
}

// make sure startID < endID
SinglePet::SinglePet( int id, int startContigID, string startContigOri, int endContigID, string endContigOri,
		  int dis, int std ){
	if( startContigID < endContigID ){
		m_startContigID = startContigID;
		m_startContigOri = PLUS;
		if( startContigOri == "-" )
			m_startContigOri = MINUS;
		m_endContigID = endContigID;
		m_endContigOri = PLUS;
		if( endContigOri == "-" )
			m_endContigOri = MINUS;
	}
	else{
		// reverse two contigs
		m_startContigID = endContigID;;
		m_startContigOri = PLUS;
		if( endContigOri == "+" )
			m_startContigOri = MINUS;
		m_endContigID = startContigID;
		m_endContigOri = PLUS;
		if( startContigOri == "+" )
			m_endContigOri = MINUS;
	}
	m_id = id;
	m_distance = dis;
	m_std = std;
}

ostream& operator<<( ostream& out, const SinglePet& p ){
	out<<p.m_id<<"\t"<<p.m_startContigID<<"\t"<<p.m_startContigOri<<"\t"<<p.m_endContigID<<"\t"
		<<p.m_endContigOri<<"\t"<<p.m_distance<<"\t"<<p.m_std;
	return out;
}
