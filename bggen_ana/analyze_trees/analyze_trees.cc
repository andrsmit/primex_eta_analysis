#include "eta.h"
#include "EtaAna.h"

int main(int argc, char **argv) {
	
	anaSettings_t locSettings;
	locSettings.versionString   = "";
	locSettings.primexPhase     = 1;
	locSettings.analysisOption  = 6;
	locSettings.vetoOption      = 6;
	locSettings.minRunNumber    = GetMinRunNumber(1);
	locSettings.maxRunNumber    = GetMaxRunNumber(1);
	locSettings.configFileName = "eta_ana.config";
	
	locSettings.singleRunNumber = 61434;
	locSettings.inputFileName  = "none";
	locSettings.outputFileName = "eta_ana.root";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'x':
				locSettings.versionString = ++argptr;
				break;
			case 'p':
				locSettings.primexPhase = atoi(++argptr);
				locSettings.minRunNumber = GetMinRunNumber(locSettings.primexPhase);
				locSettings.maxRunNumber = GetMaxRunNumber(locSettings.primexPhase);
				break;
			case 'a':
				locSettings.analysisOption = atoi(++argptr);
				break;
			case 'v':
				locSettings.vetoOption = atoi(++argptr);
				break;
			case 'c':
				locSettings.configFileName = ++argptr;
				break;
			case 'l':
				locSettings.minRunNumber = atoi(++argptr);
				break;
			case 'u':
				locSettings.maxRunNumber = atoi(++argptr);
				break;
			case 'r':
				locSettings.singleRunNumber = atoi(++argptr);
				break;
			case 'i':
				locSettings.inputFileName = ++argptr;
				break;
			case 'o':
				locSettings.outputFileName = ++argptr;
				break;
			case 't':
				locSettings.targetNumber = atoi(++argptr);
				break;
			case 'b':
				locSettings.batchNumber = atoi(++argptr);
				break;
			case 'h':
				printUsage(locSettings,0);
				break;
			default:
				fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
				printUsage(locSettings,0);
				break;
			}
		}
	}
	printUsage(locSettings, 1);
	
	char locPathName[256] = "/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana";
	
	// Construct analysis object:
	
	EtaAna locAna;
	
	TString analysisStr = "";
	switch(locSettings.analysisOption) {
		case 0:
			analysisStr = "";
			break;
		case 1:
			analysisStr = "_matrix";
			break;
		case 2:
			analysisStr = "_FCAL";
			break;
		case 3:
			analysisStr = "_BCAL";
			break;
		case 4:
			analysisStr = "_BEAM";
			break;
		case 5:
			analysisStr = "_TOF";
			break;
		case 6:
			analysisStr = "_bggen";
			break;
		default:
			std::cout << "Invalid analysis option provided." << std::endl;
			break;
	}
	
	TString vetoStr = "";
	if(locSettings.analysisOption>0) vetoStr = Form("_vetoOption%d", locSettings.vetoOption);
	
	// Read cut values from config file:
	TString locConfigFileName(locSettings.configFileName);
	locAna.SetCuts(locConfigFileName);
	locAna.DumpCuts();
	
	// Set veto option:
	locAna.SetVetoOption(locSettings.vetoOption);
	
	// Initialize histograms to be filled:
	locAna.InitHistograms(locSettings.analysisOption);
	
	int locPhase = locSettings.primexPhase;
	
	TString fieldString = locPhase==1 ? "nobfield" : "bfield";
	
	//----------------------------------------------------//
	// Get list of runs to process:
	
	// Directory where run lists are stored:
	sprintf(runListPathName, "/work/halld/home/andrsmit/primex_eta_analysis/run_lists");
	
	TString runListFileName = Form("%s/phase%d/he_full_%s.txt", runListPathName, locPhase, fieldString.Data());
	if(gSystem->AccessPathName(runListFileName.Data())) {
		printf("\n\nRun list for the specified run period does not exist.\n\n");
		return 1;
	}
	// store runs in a vector of integers:
	vector<int> runList; runList.clear();
	
	// objects used for reading txt file:
	ifstream runListStream(runListFileName.Data());
	string runNumberStr;
	
	// fill run list vector:
	while(getline(runListStream, runNumberStr)) {
		int locRun = stoi(runNumberStr);
		if(locRun>=locSettings.minRunNumber && locRun<=locSettings.maxRunNumber) {
			runList.push_back(stoi(runNumberStr));
		}
	}
	runListStream.close();
	
	// skip this setting if run list is empty:
	int locNRuns = (int)runList.size();
	if(locNRuns<1) {
		printf("\n\nNo runs found in run list.\n\n");
		return 1;
	}
	//----------------------------------------------------//
	
	TString locTargetStr = locSettings.targetNumber==0 ? "proton" : "neutron";
	
	// Directory where input ROOT trees are read from:
	if(locPhase==3) {
		sprintf(rootTreePathName, "%s/pass2/phase%d/batch%02d/Helium-%s/rootTrees", locPathName, locPhase, 
			locSettings.batchNumber, locTargetStr.Data());
	} else {
		sprintf(rootTreePathName, "%s/pass1/phase%d/%s/rootTrees", locPathName, locPhase, locSettings.versionString.c_str());
	}
	if(gSystem->AccessPathName(rootTreePathName)) {
		cout << "Input ROOT tree directory does not exist." << endl;
		cout << " directory name: " << rootTreePathName << endl;
		exit(1);
	}
	
	// Directory where output ROOT files will be stored:
	sprintf(rootFilePathName, "%s/analyze_trees/rootFiles/phase%d", locPathName, locPhase);
	if(gSystem->AccessPathName(rootFilePathName)) {
		cout << "Output ROOT file directory does not exist." << endl;
		cout << " directory name: " << rootFilePathName << endl;
		exit(1);
	}
	
	// Output file name:
	TString locOutputFileName;
	if(locPhase==3) {
		locOutputFileName = Form("%s/Helium-%s-batch%02d%s.root", rootFilePathName, locTargetStr.Data(), 
			locSettings.batchNumber, vetoStr.Data());
	} else {
		locOutputFileName = Form("%s/%s%s.root", rootFilePathName, locSettings.versionString.c_str(), vetoStr.Data());
	}
	locAna.SetOutputFileName(locOutputFileName.Data());
	
	int nRunsProcessed = 0;
	for(int iRun=0; iRun<locNRuns; iRun++) {
		
		int locRun = runList[iRun];
		TString locInputFileName = Form("%s/%06d.root", rootTreePathName, locRun);
		
		// check if file exists:
		if(gSystem->AccessPathName(locInputFileName.Data())) {
			printf("\n  ROOT Tree file does not exist for run %d. Skipping.\n", locRun);
			continue;
		}
		
		// check if flux has been processed for this run. If not, skip it:
		TString locFluxFileName = Form(
			"/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/rootFiles/phase%d/full/%d_flux.root",
			locPhase, locRun);
		if(gSystem->AccessPathName(locFluxFileName)) {
			printf("\n  Flux file does not exist for run %d. Skipping.\n", locRun);
			continue;
		}
		
		std::cout << "  processing run number " << locRun << std::endl;
		if(locAna.SetRunNumber(locRun)) continue;
		locAna.RunAnalysis(locInputFileName.Data(), locSettings.analysisOption);
		nRunsProcessed++;
	}
	if(nRunsProcessed>0) {
		locAna.WriteHistograms(locSettings.analysisOption);
	}
	
	return 0;
}

void printUsage(anaSettings_t anaSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-x<arg>\tVersion string (e.g. Helium-proton-v3895)\n");
		fprintf(stderr,"-p<arg>\tPrimEx-eta phase number (default=1, choose between 1,2,or3)\n");
		fprintf(stderr,"-a<arg>\tAnalysis Option (0: default, 1: matrix, 2: FCAL cuts, 5: TOF cuts, 6: bggen)\n");
		fprintf(stderr,"-v<arg>\tVeto Option (default: 6, shoud range between 0-8)\n");
		fprintf(stderr,"-c<arg>\tConfiguration file name\n");
		fprintf(stderr,"-l<arg>\tMinimum run number to process\n");
		fprintf(stderr,"-u<arg>\tMaximum run number to process\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Options for analyzing phase 3 simulations:\n");
		fprintf(stderr,"-t<arg>\tTarget (proton: 0 (default), neutron: 1)\n");
		fprintf(stderr,"-b<arg>\tBatch number\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Options for testing on a single file:\n");
		fprintf(stderr,"-r<arg>\tRun number for specified input file (for testing, not currently implemented)\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is eta_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(anaSettings.inputFileName!="none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", anaSettings.singleRunNumber, 
				anaSettings.inputFileName.c_str());
			cout << "" << endl;
		} else {
			printf("\nAnalyzing %s simulations for runs %d-%d:\n", anaSettings.versionString.c_str(), 
				anaSettings.minRunNumber, anaSettings.maxRunNumber);
			cout << "" << endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
