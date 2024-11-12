#ifndef ETA_INCLUDE
#define ETA_INCLUDE

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <chrono>

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
	int primexPhase;
	int minRunNumber;
	int maxRunNumber;
	int runNumber;
	string  inputFileName;
	string outputFileName;
	string configFileName;
};

int GetMinRunNumber(int phase) {
	if(phase==1)      return  61462;
	else if(phase==2) return  81396;
	else if(phase==3) return 110641;
	else return -1;
};
int GetMaxRunNumber(int phase) {
	if(phase==1)      return  61956;
	else if(phase==2) return  81716;
	else if(phase==3) return 112001;
	else return -1;
};

void printUsage(anaSettings_t, int goYes);

int main(int argc, char **argv);

//----------   Data Objects   ----------//

char rootTreePathName[256];
char rootFilePathName[256];

#endif
