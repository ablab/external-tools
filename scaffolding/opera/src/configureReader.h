#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "Configure.h"

using namespace std;

class configureReader
{
public:
	configureReader(void);
	~configureReader(void);

	// Methods
public:
	// Read the configuration file
	int ReadConfigFile( string fileName );

private:
	// anylize the parameters
	int AnalyzeParameters( string line );
};
