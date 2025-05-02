#include "eta.h"
#include "EtaAna.h"

int main(int argc, char **argv) {
	
	anaSettings_t locSettings;
	locSettings.primexPhase    =     1;
	locSettings.mcBatch        =     1;
	locSettings.analysisOption =     0;
	locSettings.vetoOption     =     6;
	locSettings.runNumber      =     0;
	locSettings.minRunNumber   = GetMinRunNumber(locSettings.primexPhase);
	locSettings.maxRunNumber   = GetMaxRunNumber(locSettings.primexPhase);
	locSettings.inputFileName  = "none";
	locSettings.outputFileName = "eta_ana.root";
	locSettings.configFileName = "eta_ana.config";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'p':
				locSettings.primexPhase  = atoi(++argptr);
				locSettings.minRunNumber = GetMinRunNumber(locSettings.primexPhase);
				locSettings.maxRunNumber = GetMaxRunNumber(locSettings.primexPhase);
				break;
			case 'b':
				locSettings.mcBatch = atoi(++argptr);
				break;
			case 'c':
				locSettings.configFileName = ++argptr;
				break;
			case 'a':
				locSettings.analysisOption = atoi(++argptr);
				break;
			case 'v':
				locSettings.vetoOption = atoi(++argptr);
				break;
			case 'i':
				locSettings.inputFileName = ++argptr;
				break;
			case 'o':
				locSettings.outputFileName = ++argptr;
				break;
			case 'r':
				locSettings.runNumber = atoi(++argptr);
				break;
			case 'l':
				locSettings.minRunNumber = atoi(++argptr);
				break;
			case 'u':
				locSettings.maxRunNumber = atoi(++argptr);
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
	
	// Directory where output ROOT files will be stored:
	sprintf(rootFilePathName, "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles");
	
	// Directory where run lists are stored:
	sprintf(runListPathName, "/work/halld/home/andrsmit/primex_eta_analysis/run_lists");
	
	// Initialize analysis object:
	EtaAna locAna;
	locAna.SetFillThrown(true);
	locAna.SetIsCohMC(true);
	
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
		case 8:
			analysisStr = "_angular";
			break;
		case 9:
			analysisStr = "_special";
			break;
		default:
			std::cout << "Invalid analysis option provided." << std::endl;
			break;
	}
	
	TString vetoStr = "";
	if(locSettings.analysisOption>0) vetoStr = Form("_VetoOption%d", locSettings.vetoOption);
	
	// Read cut values from config file:
	TString locConfigFileName(locSettings.configFileName);
	locAna.SetCuts(locConfigFileName);
	locAna.DumpCuts();
	
	// Set veto option:
	locAna.SetVetoOption(locSettings.vetoOption);
	
	// Initialize histograms to be filled:
	locAna.InitHistograms(locSettings.analysisOption);
	
	int locPhase = locSettings.primexPhase;
	
	TString fieldString = "bfield";
	if(locPhase==1) {
		fieldString = "nobfield";
	}
	
	// Directory where input ROOT trees will be read from:
	sprintf(rootTreePathName, "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/pass1/rootTrees/phase%d/full_target_%s/batch%02d", 
		locPhase, fieldString.Data(), locSettings.mcBatch);
	
	//----------------------------------------------------//
	// Get list of runs to process:
	
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
		runList.push_back(stoi(runNumberStr));
	}
	runListStream.close();
	
	// skip this setting if run list is empty:
	int locNRuns = (int)runList.size();
	if(locNRuns<1) {
		printf("\n\nNo runs found in run list.\n\n");
		return 1;
	}
	//----------------------------------------------------//
	
	TString outputFileName = Form("%s/phase%d/phase%d%s%s_batch%02d.root", rootFilePathName, locPhase, locPhase, 
		analysisStr.Data(), vetoStr.Data(), locSettings.mcBatch);
	
	//if(!gSystem->AccessPathName(outputFileName.Data())) return 0;
	locAna.SetOutputFileName(outputFileName.Data());
	
	int nRunsProcessed = 0;
	for(int iRun=0; iRun<locNRuns; iRun++) {
		int locRun = runList[iRun];
		if((locRun<locSettings.minRunNumber) || (locRun>locSettings.maxRunNumber)) continue;
		
		TString inputFileName = Form("%s/%d.root", rootTreePathName, locRun);
		
		// check if file exists:
		if(gSystem->AccessPathName(inputFileName.Data())) {
			printf("\n  ROOT Tree file does not exist for run %d. Skipping.\n", locRun);
			continue;
		}
		
		// check if flux has been processed for this run. If not, skip it:
		TString fluxFileName = Form(
			"/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/rootFiles/phase%d/full/%d_flux.root",
			locPhase, locRun);
		if(gSystem->AccessPathName(fluxFileName)) {
			printf("\n  Flux file does not exist for run %d. Skipping.\n", locRun);
			continue;
		}
		
		std::cout << "  processing run number " << locRun << std::endl;
		if(locAna.SetRunNumber(locRun)) continue;
		locAna.RunAnalysis(inputFileName.Data(), locSettings.analysisOption);
		nRunsProcessed++;
	}
	
	if(nRunsProcessed>0) {
		locAna.WriteHistograms(locSettings.analysisOption);
		locAna.ResetHistograms(locSettings.analysisOption);
	}
	
	return 0;
}

void printUsage(anaSettings_t anaSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-p<arg>\tPrimEx-eta phase number (default=1, choose between 1,2,or3)\n");
		fprintf(stderr,"-b<arg>\tMC Batch to analyze (0-9)\n");
		fprintf(stderr,"-c<arg>\tConfiguration file name\n");
		fprintf(stderr,"-a<arg>\tAnalysis Option (0: default, 1: matrix, 2: FCAL cuts, 5: TOF cuts, 8: angular systematics)\n");
		fprintf(stderr,"-v<arg>\tVeto Option (default: 6, shoud range between 0-8)\n");
		fprintf(stderr,"\nFor debug/testing:\n");
		fprintf(stderr,"-r<arg>\tRun number for specified input file\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is compton_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(anaSettings.inputFileName != "none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", anaSettings.runNumber, 
				anaSettings.inputFileName.c_str());
			std::cout << "" << std::endl;
		} else {
			printf("\nAnalyzing simulations for PrimEx-eta Phase %d (batch %d):\n", 
				anaSettings.primexPhase, anaSettings.mcBatch);
			std::cout << "" << std::endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
