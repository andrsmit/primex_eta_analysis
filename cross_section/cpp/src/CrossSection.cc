#include "CrossSection.h"

int main(int argc, char **argv) 
{
	configSettings_t locConfig;
	
	TString configFileName = "";
	
	// parse command line:
	char *argptr;
	for(int iarg=1; iarg<argc; iarg++) {
		argptr = argv[iarg];
		if(*argptr == '-') {
			argptr++;
			switch(*argptr) {
			case 'p':
				locConfig.primexPhase = atoi(++argptr);
				break;
			case 'a':
				locConfig.analysisOption = atoi(++argptr);
				break;
			case 'c':
				configFileName = ++argptr;
				break;
			case 'h':
				printUsage(locConfig,0);
				break;
			default:
				fprintf(stderr,"Unrecognized argument: [-%s]\n\n",argptr);
				printUsage(locConfig,0);
				break;
			}
		}
	}
	printUsage(locConfig, 1);
	
	//------------------------------------------------//
	
	// Object to extract angular yield:
	
	EtaAnalyzer locEtaAna(locConfig.primexPhase);
	
	// Initialize analyzer with settings from config file:
	
	if(LoadConfigSettings(locEtaAna, configFileName)) exit(1);
	
	locEtaAna.LoadLineshapes(locConfig.analysisOption);
	
	// Print to console:
	
	locEtaAna.DumpSettings();
	
	//------------------------------------------------//
	
	locEtaAna.CalcLuminosity();
	locEtaAna.CalcEmptyTargetFluxRatio();
	
	if(locEtaAna.GetDataHistograms(locConfig.analysisOption)) {
		cout << "\n\nProblem accessing data.\n\n" << endl;
		exit(1);
	}
	
	//------------------------------------------------//
	
	// needed for ROOT graphics:
	TApplication theApp("App", &argc, argv);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	//locEtaAna.FitInvariantMass(0.30, 0.42, 1);
	
	if(locEtaAna.GetAcceptance()) {
		cout << "\n\nProblem getting acceptance.\n\n" << endl;
	}
	
	locEtaAna.ExtractAngularYield(1);
	locEtaAna.PlotAngularYield();
	locEtaAna.PlotCrossSection();
	
	gSystem->ProcessEvents();
	theApp.Run();
	
	return 0;
}

int ReadConfigFile(TString fileName)
{
	if(gSystem->AccessPathName(fileName.Data())) return 1;
	
	return 0;
}
void printUsage(configSettings_t settings, int goYes) 
{
	if(goYes==0) {
		fprintf(stderr,"\nCommand line arguments:\n");
		fprintf(stderr,"-h\tPrintthis message\n");
		fprintf(stderr,"-p<arg>\tPrimEx-eta phase number (default=1, choose between 1,2,or3)\n");
		fprintf(stderr,"-a<arg>\tAnalysis Option (0: default, 1: invmass matrix, 2: FCAL cuts, 3: BCAL cuts, 4: BEAM cuts, 5: TOF cuts)\n");
		fprintf(stderr,"-c<arg>\tConfiguration file name\n");
	}
	
	return;
}

int LoadConfigSettings(EtaAnalyzer &anaObject, TString fileName)
{
	if(fileName=="") {
		// check if a config file exists in current directory:
		fileName = "eta_cs.config";
	} else {
		printf("config filename = [%s]\n", fileName.Data());
	}
	
	MyReadConfig *ReadFile = new MyReadConfig();
	
	if(gSystem->AccessPathName(fileName.Data())) {
		cout << "\nUnable to access config file. Use option '-c' to specify it.\n" << endl;
		return 1;
	} else {
		if(ReadFile->ReadConfigFile(fileName)) return 1;
	}
	
	// Load Configuration Settings:
	
	// Beam Energy:
	if(ReadFile->GetConfigName("minBeamEnergy") != "") {
		anaObject.minBeamEnergy = ReadFile->GetConfig1Par("minBeamEnergy")[0];
	}
	if(ReadFile->GetConfigName("maxBeamEnergy") != "") {
		anaObject.maxBeamEnergy = ReadFile->GetConfig1Par("maxBeamEnergy")[0];
	}
	
	// Bin sizing:
	if(ReadFile->GetConfigName("rebinsMgg") != "") {
		anaObject.SetRebinsMgg(ReadFile->GetConfig1Par("rebinsMgg")[0]);
	}
	if(ReadFile->GetConfigName("rebinsTheta") != "") {
		anaObject.SetRebinsTheta(ReadFile->GetConfig1Par("rebinsTheta")[0]);
	}
	
	//
	// Function options:
	//
	if(ReadFile->GetConfigName("signalFitOption") != "") {
		anaObject.SetFitOption_signal((int)ReadFile->GetConfig1Par("signalFitOption")[0]);
	}
	if(ReadFile->GetConfigName("bkgdFitOption") != "") {
		if(ReadFile->GetConfigName("polynomialOrder") != "") {
			anaObject.SetFitOption_bkgd((int)ReadFile->GetConfig1Par("bkgdFitOption")[0], 
				(int)ReadFile->GetConfig1Par("polynomialOrder")[0]);
		}
		anaObject.SetFitOption_bkgd((int)ReadFile->GetConfig1Par("bkgdFitOption")[0]);
	}
	if(ReadFile->GetConfigName("omegaFitOption") != "") {
		anaObject.SetFitOption_omega((int)ReadFile->GetConfig1Par("omegaFitOption")[0]);
	}
	if(ReadFile->GetConfigName("etapFitOption") != "") {
		anaObject.SetFitOption_etap((int)ReadFile->GetConfig1Par("etapFitOption")[0]);
	}
	
	// Fitting range:
	if(ReadFile->GetConfigName("fittingRange") != "") {
		anaObject.SetFitRange((double)ReadFile->GetConfig2Par("fittingRange")[0], 
				(double)ReadFile->GetConfig2Par("fittingRange")[1]);
	}
	
	return 0;
}
