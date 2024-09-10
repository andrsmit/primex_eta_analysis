#ifndef _ETA_INCLUDE_
#define _ETA_INCLUDE_

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

struct genSettings_t {
	int min_run;
	int max_run;
	int run_number;       //
	string  input_fname;
	string output_fname;
};

void printUsage(genSettings_t, int goYes);

int main(int argc, char **argv);

//----------   Data Objects   ----------//

char rootTree_pathName[256], rootFile_pathName[256];

#endif
