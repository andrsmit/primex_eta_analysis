
#include "eta.h"
#include "EtaAna.h"

int main(int argc, char **argv) {
	
	genSettings_t genSettings;
	genSettings.min_run      = 61434;
	genSettings.max_run      = 61434;
	genSettings.run_number   = 61434;
	genSettings.input_fname  = "none";
	genSettings.output_fname = "eta_ana.root";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'b':
				genSettings.min_run  = atoi(++argptr);
				break;
			case 'e':
				genSettings.max_run = atoi(++argptr);
				break;
			case 'r':
				genSettings.run_number = atoi(++argptr);
				break;
			case 'i':
				genSettings.input_fname = ++argptr;
				break;
			case 'o':
				genSettings.output_fname = ++argptr;
				break;
			case 'h':
				printUsage(genSettings,0);
				break;
			default:
				fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
				printUsage(genSettings,0);
				break;
			}
		}
	}
	printUsage(genSettings, 1);
	
	char loc_path[256] = "/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana";
	
	// Construct analysis object:
	
	EtaAna locAna;
	
	// Initialize histograms to be filled:
	
	locAna.initHistograms();
	
	TString ver_str = "Helium-proton-v3895";
	
	//
	// Check if an input filename was specificed at runtime. 
	// If not, we'll do a loop over TAGGER counters from the rootTree directory above:
	//
	if(genSettings.input_fname!="none") {
		
		if(genSettings.run_number==0) {
			fprintf(stderr,"Need to specify the runnumber for the file you are processing.\n");
			exit(0);
		}
		
		TString input_fname = Form("%s", genSettings.input_fname.c_str());
		if(gSystem->AccessPathName(input_fname.Data())) {
			fprintf(stderr,"Specified input filename is inaccessible.\n");
			exit(0);
		}
		locAna.setOutputFileName(Form("%s",genSettings.output_fname.c_str()));
		locAna.resetHistograms();
		locAna.runAnalysis(input_fname.Data());
		locAna.writeHistograms();
		
	} else {
		
		int loc_phase = locAna.getPrimexPhase(genSettings.min_run);
		
		// Directory where input ROOT trees are read from:
		sprintf(rootTree_pathName, "%s/rootTrees/phase%d/%s", loc_path, loc_phase, ver_str.Data());
		if(gSystem->AccessPathName(rootTree_pathName)) {
			cout << "Input ROOT tree directory does not exist." << endl;
			cout << " directory name: " << rootTree_pathName << endl;
			exit(1);
		}
		
		// Directory where output ROOT files will be stored:
		sprintf(rootFile_pathName, "%s/analyze_trees/rootFiles/phase%d", loc_path, loc_phase);
		if(gSystem->AccessPathName(rootFile_pathName)) {
			cout << "Output ROOT file directory does not exist." << endl;
			cout << " directory name: " << rootFile_pathName << endl;
			exit(1);
		}
		
		// Output file name:
		TString output_fname = Form("%s/%s_%d_%d.root", rootFile_pathName, ver_str.Data(), 
			genSettings.min_run, genSettings.max_run);
		locAna.setOutputFileName(output_fname.Data());
		
		for(int loc_run = genSettings.min_run; loc_run <= genSettings.max_run; loc_run++) {
			
			locAna.setRunNumber(loc_run);
			
			auto start = chrono::high_resolution_clock::now();
			
			TString input_fname = Form("%s/%06d.root", rootTree_pathName, loc_run);
			
			// check if file exists:
			if(gSystem->AccessPathName(input_fname.Data())) continue;
			
			cout << "  processing run " << loc_run << endl;
			locAna.runAnalysis(input_fname.Data());
			
			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
			//cout << duration.count() << endl;
		}
		cout << "Done." << endl;
		
		locAna.writeHistograms();
	}
	
	return 0;
}

void printUsage(genSettings_t genSettings, int goYes) {
	
	if(goYes==0) {
		fprintf(stderr,"\nSWITCHES:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-b<arg>\tMinimum run number to process\n");
		fprintf(stderr,"-e<arg>\tMaximum run number to process\n");
		fprintf(stderr,"-r<arg>\tRun number for specified input file\n");
		fprintf(stderr,"-i<arg>\tInput file name (default is none)\n");
		fprintf(stderr,"-o<arg>\tOutput file name (default is eta_ana.root)\n\n\n");
	}
	
	if(goYes==1) {
		if(genSettings.input_fname!="none") {
			printf("\nAnalyzing simulations for run %d, input file: %s\n", genSettings.run_number, 
				genSettings.input_fname.c_str());
			cout << "" << endl;
		} else {
			printf("\nAnalyzing simulations for runs %d-%d:\n", genSettings.min_run, genSettings.max_run);
			cout << "" << endl;
		}
	}
	
	if(goYes==0) exit(0);
	
	return;
}
