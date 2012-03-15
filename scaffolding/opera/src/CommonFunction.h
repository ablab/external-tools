#pragma once

#include <string.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "Configure.h"

using namespace std;

// convert int to string
string itos( int i );
// split the string
void Split( string line, string label, vector<string> *column );
void Split( string line, string label, list<string> *column );
// get the opposite orientation 
int GetOppositeOri( int ori );
// print thousand comma format of a number
string PrintNumber( double num );
// print thousand comma format of a number to file
void PrintNumberToFile( FILE *file, double num );
// check if it is a number (coverage)
bool IsNumber( string s );
