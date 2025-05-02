#ifndef _BGGEN_ANA_INCLUDE_
#define _BGGEN_ANA_INCLUDE_

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSystem.h"
#include "TDirectory.h"

//---------------------------------------------//
// Structure for the configuration settings:

struct anaSettings_t {
	string versionString;
	int primexPhase;
	int analysisOption;
	int vetoOption;
	int minRunNumber;
	int maxRunNumber;
	string configFileName;
	
	// only used for analyzing phase 3 simulations:
	int batchNumber  = 0;
	int targetNumber = 0;
	
	// only used for testing code over single file:
	int singleRunNumber;
	string  inputFileName;
	string outputFileName;
};

int GetMinRunNumber(int phase) {
	if(phase==1)      return  61378;
	else if(phase==2) return  81396;
	else if(phase==3) return 110657;
	else return -1;
};
int GetMaxRunNumber(int phase) {
	if(phase==1)      return  61956;
	else if(phase==2) return  81716;
	else if(phase==3) return 111957;
	else return -1;
};

void printUsage(anaSettings_t, int goYes);

int main(int argc, char **argv);

//----------   Data Objects   ----------//

char rootTreePathName[256], rootFilePathName[256], runListPathName[256];

#endif
