#include "CrossSection.h"

int main(int argc, char **argv) 
{
	configSettings_t locConfig;
	
	TString configFileName = "";
	
	int locDrawOption = 0;
	
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
			case 'c':
				configFileName = ++argptr;
				break;
			case 'h':
				printUsage(locConfig,0);
				break;
			case 'd':
				locDrawOption = 1;
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
	
	EtaAnalyzer locEtaAna;
	
	// Initialize analyzer with settings from config file:
	
	if(LoadConfigSettings(locEtaAna, configFileName, locConfig.primexPhase)) exit(1);
	
	if(locEtaAna.LoadLineshapes()) {
		cout << "\n\nProblem loading MC lineshapes.\n\n" << endl;
		exit(1);
	}
	
	// Print to console:
	
	locEtaAna.DumpSettings();
	
	//------------------------------------------------//
	
	// needed for ROOT graphics:
	TApplication theApp("App", &argc, argv);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	if(locEtaAna.LoadLuminosity()) {
		printf("\n\nProblem getting integrated luminosity.\n\n");
		exit(1);
	}
	if(locEtaAna.LoadEmptyTargetFluxRatio()) {
		printf("\n\nProblem loading empty target photon flux ratio.\n\n");
		exit(1);
	}
	if(locEtaAna.LoadDataHistograms()) {
		printf("\n\nProblem accessing data.\n\n");
		exit(1);
	}
	if(locEtaAna.LoadAngularMatrix()) {
		printf("\n\nProblem loading angular matrices.\n\n");
		exit(1);
	}
	//------------------------------------------------//
	
	//locEtaAna.FitInvariantMass(0.4, 0.5, 1);
	if(1) {
		locEtaAna.ExtractAngularYield(locDrawOption);
		if(locEtaAna.CalcAcceptance()) {
			cout << "\n\nProblem getting acceptance.\n\n" << endl;
		}
		locEtaAna.PlotAngularYield();
		locEtaAna.PlotCrossSection();
	}
	
	YieldFitter locFitter;
	InitializeYieldFitter(locFitter, locEtaAna);
	locFitter.FitAngularYield();
	
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

int LoadConfigSettings(EtaAnalyzer &anaObject, TString fileName, int phase)
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
	
	// Primex phase:
	if(ReadFile->GetConfigName("primexPhase") != ""){
		int configPhase = (int)ReadFile->GetConfig1Par("primexPhase")[0];
		if((phase>0) && (phase != configPhase)) {
			anaObject.SetPhase(phase);
		} else {
			anaObject.SetPhase(configPhase);
		}
	}
	
	if(ReadFile->GetConfigName("analysisOption") != ""){
		anaObject.SetAnalysisOption((int)ReadFile->GetConfig1Par("analysisOption")[0]);
	}
	if(ReadFile->GetConfigName("histName") != ""){
		anaObject.SetMggHistName(ReadFile->GetConfigName("histName"));
	}
	if(ReadFile->GetConfigName("matrixName") != ""){
		anaObject.SetMatrixHistName(ReadFile->GetConfigName("matrixName"));
	}
	
	// Beam Energy:
	if(ReadFile->GetConfigName("beamEnergy") != "") {
		anaObject.SetBeamEnergy((double)ReadFile->GetConfig2Par("beamEnergy")[0], 
				(double)ReadFile->GetConfig2Par("beamEnergy")[1]);
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
	
	
	//
	// Empty target options:
	//
	if(ReadFile->GetConfigName("subtractEmptyTarget") != "") {
		anaObject.SetSubtractEmptyTarget((int)ReadFile->GetConfig1Par("subtractEmptyTarget")[0]);
	}
	if(ReadFile->GetConfigName("fitEmptyTarget") != "") {
		anaObject.SetFitEmptyTarget((int)ReadFile->GetConfig1Par("fitEmptyTarget")[0]);
	}
	
	if(ReadFile->GetConfigName("emptyFitOption_eta") != "") {
		anaObject.SetEmptyFitOption_eta((int)ReadFile->GetConfig1Par("emptyFitOption_eta")[0]);
	}
	if(ReadFile->GetConfigName("emptyFitOption_omega") != "") {
		anaObject.SetEmptyFitOption_omega((int)ReadFile->GetConfig1Par("emptyFitOption_omega")[0]);
	}
	if(ReadFile->GetConfigName("emptyFitOption_fdc") != "") {
		anaObject.SetEmptyFitOption_fdc((int)ReadFile->GetConfig1Par("emptyFitOption_fdc")[0]);
	}
	if(ReadFile->GetConfigName("emptyFitOption_bkgd") != "") {
		if(ReadFile->GetConfigName("emptyFitOption_poly") != "") {
			anaObject.SetEmptyFitOption_bkgd((int)ReadFile->GetConfig1Par("emptyFitOption_bkgd")[0], 
				(int)ReadFile->GetConfig1Par("emptyFitOption_poly")[0]);
		}
		anaObject.SetEmptyFitOption_bkgd((int)ReadFile->GetConfig1Par("emptyFitOption_bkgd")[0]);
	}
	
	if(ReadFile->GetConfigName("emptyFittingRange") != "") {
		anaObject.SetEmptyFitRange((double)ReadFile->GetConfig2Par("emptyFittingRange")[0], 
				(double)ReadFile->GetConfig2Par("emptyFittingRange")[1]);
	}
	
	return 0;
}
