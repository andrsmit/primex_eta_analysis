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
	
	// Print to console:
	
	locEtaAna.DumpSettings();
	
	//------------------------------------------------//
	
	// needed for ROOT graphics:
	TApplication theApp("App", &argc, argv);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	if(locEtaAna.LoadLineshapes()) {
		cout << "\n\nProblem loading MC lineshapes.\n\n" << endl;
		exit(1);
	}
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
	
	YieldFitter locFitter;
	InitializeYieldFitter(locFitter, locEtaAna);
	
	if(1) {
		locEtaAna.ExtractAngularYield(locDrawOption);
		if(locEtaAna.CalcAcceptance()) {
			cout << "\n\nProblem getting acceptance.\n\n" << endl;
		}
		locEtaAna.PlotAngularYield();
		locEtaAna.PlotCrossSection();
		if(locEtaAna.GetEmptyFitOption(0) && locEtaAna.GetEmptyFitOption(1)) {
			locEtaAna.PlotEmptyEtaRatio();
		}
		if(locEtaAna.GetFitOption(1)==7) {
			locEtaAna.PlotEtaPionFraction();
		}
		locEtaAna.WriteROOTFile(Form("output/yield_phase%d.root",locEtaAna.GetPhase()));
		locFitter.SetYield((TH1F*)locEtaAna.GetAngularYield(1));
	}
	else {
		TFile *fIn = new TFile("output/yield_phase3_fitEmpty.root", "READ");
		TH1F *hYield = (TH1F*)fIn->Get("AngularYieldFit");
		hYield->SetDirectory(0);
		fIn->Close();
		locFitter.SetYield(hYield);
	}
	
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
	/*
	Note: The order in which configuration options are read can be important.
	For example, if 'analysisOption'=0, we then need to look for 'vetoOption'. 
	But if 'analysisOption' is anything else, then 'vetoOption' is unimportant.
	Similarly, the beam energy range should only be changed from the hard-coded defaults
	if 'analysisOption'=1
	*/
	
	
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
	
	//----------------------------------------//
	// Primex phase:
	
	if(ReadFile->GetConfigName("primexPhase") != ""){
		int configPhase = (int)ReadFile->GetConfig1Par("primexPhase")[0];
		if((phase>0) && (phase != configPhase)) {
			anaObject.SetPhase(phase);
		} else {
			anaObject.SetPhase(configPhase);
		}
	}
	
	//----------------------------------------//
	//
	// Analysis and Veto Options:
	//   (mgg histogram and angular matrix names are specified here as well)
	//
	
	if(ReadFile->GetConfigName("analysisOption") != "")
	{
		anaObject.SetAnalysisOption((int)ReadFile->GetConfig1Par("analysisOption")[0]);
		if(ReadFile->GetConfigName("vetoOption") != "")
		{
			anaObject.SetVetoOption((int)ReadFile->GetConfig1Par("vetoOption")[0]);
		}
	}
	
	switch(anaObject.GetAnalysisOption()) {
		case 0:
			anaObject.SetMggHistName(Form("VetoOption%d/mgg_const_veto_%d",
				anaObject.GetVetoOption(), anaObject.GetVetoOption()));
			anaObject.SetMatrixHistName(Form("AngularMatrix/AngularMatrix_veto_%d",
				anaObject.GetVetoOption()));
			break;
		case 1:
			anaObject.SetMatrixHistName("AngularMatrix");
			break;
		default:
		{
			if(ReadFile->GetConfigName("histName") != "")
			{
				anaObject.SetMggHistName(ReadFile->GetConfigName("histName"));
			}
			else {
				printf("\nNo histogram name specified for analysis option: %d. Exiting.\n", 
					anaObject.GetAnalysisOption());
				exit(1);
			}
			
			if(ReadFile->GetConfigName("matrixName") != "")
			{
				anaObject.SetMatrixHistName(ReadFile->GetConfigName("matrixName"));
			}
			else {
				printf("\nNo matrix name specified for analysis option: %d. Exiting.\n", 
					anaObject.GetAnalysisOption());
				exit(1);
			}
		}
	}
	
	//----------------------------------------//
	// Beam Energy should only be taken from config file if analysisOption=1:
	
	if(anaObject.GetAnalysisOption()==1)
	{
		if(ReadFile->GetConfigName("beamEnergy") != "")
		{
			anaObject.SetBeamEnergy((double)ReadFile->GetConfig2Par("beamEnergy")[0], 
				(double)ReadFile->GetConfig2Par("beamEnergy")[1]);
		}
	}
	
	//----------------------------------------//
	// Bin sizing:
	
	if(ReadFile->GetConfigName("rebinsMgg") != "") {
		anaObject.SetRebinsMgg(ReadFile->GetConfig1Par("rebinsMgg")[0]);
	}
	if(ReadFile->GetConfigName("rebinsTheta") != "") {
		anaObject.SetRebinsTheta(ReadFile->GetConfig1Par("rebinsTheta")[0]);
	}
	
	//----------------------------------------//
	//
	// Fitting function options:
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
	
	//----------------------------------------//
	//
	// Empty target options:
	//
	
	if(ReadFile->GetConfigName("subtractEmptyTarget") != "") {
		anaObject.SetSubtractEmptyTarget((int)ReadFile->GetConfig1Par("subtractEmptyTarget")[0]);
	}
	if(!anaObject.GetEmptySubtractOption())
	{
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
		if(ReadFile->GetConfigName("rebinsEmptyMgg") != "") {
			anaObject.SetRebinsEmptyMgg(ReadFile->GetConfig1Par("rebinsEmptyMgg")[0]);
		}
	}
	
	return 0;
}
