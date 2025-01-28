#include "eta.h"
#include "EtaAna.h"

int main(int argc, char **argv) {
	
	anaSettings_t locSettings;
	locSettings.primexPhase    =     1;
	locSettings.minRunNumber   = 61355;
	locSettings.maxRunNumber   = 61956;
	locSettings.runNumber      =     0;
	locSettings.analysisOption =     0;
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
				locSettings.primexPhase = atoi(++argptr);
				locSettings.minRunNumber = GetMinRunNumber(locSettings.primexPhase);
				locSettings.maxRunNumber = GetMaxRunNumber(locSettings.primexPhase);
				break;
			case 'b':
				locSettings.minRunNumber = atoi(++argptr);
				break;
			case 'e':
				locSettings.maxRunNumber = atoi(++argptr);
				break;
			case 'c':
				locSettings.configFileName = ++argptr;
				break;
			case 'a':
				locSettings.analysisOption = atoi(++argptr);
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
	
	// Directory where input ROOT trees will be read from:
	int locPass = 2;
	if(locSettings.primexPhase==3) locPass = 1;
	sprintf(rootTreePathName, "/work/halld/home/andrsmit/primex_eta_analysis/data/pass%d/rootTrees", locPass);
	
	// Directory where output ROOT files will be stored:
	sprintf(rootFilePathName, "/work/halld/home/andrsmit/primex_eta_analysis/analyze_trees/rootFiles");
	
	// Directory where run lists are stored:
	sprintf(runListPathName, "/work/halld/home/andrsmit/primex_eta_analysis/run_lists");
	
	// Initialize analysis object:
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
			analysisStr = "_omega";
			break;
		default:
			std::cout << "Invalid analysis option provided." << std::endl;
			exit(0);
	}
	
	// Read cut values from config file:
	TString locConfigFileName(locSettings.configFileName);
	locAna.SetCuts(locConfigFileName);
	locAna.DumpCuts();
	
	// Initialize histograms to be filled:
	locAna.InitHistograms(locSettings.analysisOption);
	
	//
	// Check if an input filename was specificed at runtime. 
	// If not, we'll do a loop over the specified run numbers:
	//
	if(locSettings.inputFileName != "none") {
		
		if(locSettings.runNumber==0) {
			fprintf(stderr,"Need to specify the runnumber for the file you are processing.\n");
			exit(0);
		}
		TString inputFileName = Form("%s", locSettings.inputFileName.c_str());
		if(gSystem->AccessPathName(inputFileName.Data())) {
			fprintf(stderr,"Specified input filename is inaccessible.\n");
			exit(0);
		}
		if(locAna.SetRunNumber(locSettings.runNumber)) return 0;
		locAna.SetOutputFileName(Form("%s",locSettings.outputFileName.c_str()));
		locAna.RunAnalysis(inputFileName.Data(), locSettings.analysisOption);
		locAna.WriteHistograms(locSettings.analysisOption);
		
	} else {
		
		int locPhase = locSettings.primexPhase;
		
		vector<TString> targetStrings = {"empty","full"};
		
		vector<TString>  fieldStrings; fieldStrings.clear();
		if(locPhase==1) {
			fieldStrings.push_back("nobfield");
		} else if(locPhase==3) {
			fieldStrings.push_back("bfield");
		} else {
			fieldStrings.push_back("nobfield");
			fieldStrings.push_back("bfield");
		}
		
		for(int itarget=0; itarget<targetStrings.size(); itarget++) {
			for(int ifield=0; ifield<fieldStrings.size(); ifield++) {
				
				//----------------------------------------------------//
				// Get list of runs to process:
				
				TString runListFileName = Form("%s/phase%d/he_%s_%s.txt", 
					runListPathName, locPhase, targetStrings[itarget].Data(), fieldStrings[ifield].Data());
				
				// skip this target+field setting if no run list exists:
				if(gSystem->AccessPathName(runListFileName.Data())) continue;
				
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
				if(locNRuns<1) continue;
				
				//----------------------------------------------------//
				
				printf("Processing %s target %s runs:\n", targetStrings[itarget].Data(), fieldStrings[ifield].Data());
				
				TString outputFileName = Form("%s/phase%d/%s_target_%s%s.root", rootFilePathName, locPhase, 
					targetStrings[itarget].Data(), fieldStrings[ifield].Data(), analysisStr.Data());
				
				//if(!gSystem->AccessPathName(outputFileName.Data())) continue;
				locAna.SetOutputFileName(outputFileName.Data());
				
				char locRootTreePathName[256];
				sprintf(locRootTreePathName, "%s/phase%d/%s_target_%s", rootTreePathName, locPhase, 
					targetStrings[itarget].Data(), fieldStrings[ifield].Data());
				
				int nRunsProcessed = 0;
				for(int iRun=0; iRun<=locNRuns; iRun++) {
					
					int locRun = runList[iRun];
					TString inputFileName = Form("%s/%06d.root", locRootTreePathName, locRun);
					
					// check if flux has been processed for this run. If not, skip it:
					TString fluxFileName = Form(
						"/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/rootFiles/phase%d/%s/%d_flux.root",
						locPhase, targetStrings[itarget].Data(), locRun);
					
					// check if flux file exists, but not rootTree:
					if(!gSystem->AccessPathName(fluxFileName.Data()) && gSystem->AccessPathName(inputFileName.Data())) {
						std:cout << "Flux file exists, but root tree does not for run " << locRun << std::endl;
						continue;
					}
					
					// check if file exists:
					if(gSystem->AccessPathName(inputFileName.Data())) continue;
					if(gSystem->AccessPathName(fluxFileName.Data())) continue;
					
					std::cout << "  processing run number " << locRun << std::endl;
					if(locAna.SetRunNumber(locRun)) continue;
					locAna.RunAnalysis(inputFileName.Data(), locSettings.analysisOption);
					nRunsProcessed++;
				}
				
				if(nRunsProcessed>0) {
					locAna.WriteHistograms(locSettings.analysisOption);
					locAna.ResetHistograms(locSettings.analysisOption);
				}
			}
		}
	}
	
	return 0;
}

void printUsage(anaSettings_t anaSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-p<arg>\tPrimEx-eta phase number (default=1, choose between 1,2,or3)\n");
		fprintf(stderr,"-b<arg>\tMinimum run number to process\n");
		fprintf(stderr,"-e<arg>\tMaximum run number to process\n");
		fprintf(stderr,"-c<arg>\tConfiguration file name\n");
		fprintf(stderr,"-a<arg>\tAnalysis Option (0: default, 1: invmass matrix, 2: FCAL cuts, 3: BCAL cuts, 4: BEAM cuts, 5: TOF cuts)\n");
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
			printf("\nAnalyzing simulations for PrimEx-eta Phase %d (runs %d-%d):\n", 
				anaSettings.primexPhase, anaSettings.minRunNumber, anaSettings.maxRunNumber);
			std::cout << "" << std::endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
