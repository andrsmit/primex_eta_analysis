#include "CrossSection.h"

int main(int argc, char **argv) 
{
	configSettings_t locConfig;
	
	TString configFileName = "";
	
	int locDrawOption = 0;
	int batchMode     = 0;
	
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
			case 'b':
				batchMode = 1;
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
	gStyle->SetOptFit(111);
	gStyle->SetOptFit(0);
	
	if(locEtaAna.LoadLuminosity()) {
		printf("\n\nProblem getting integrated luminosity.\n\n");
		exit(1);
	}
	if(locEtaAna.LoadEmptyTargetFluxRatio()) {
		printf("\n\nProblem loading empty target photon flux ratio.\n\n");
		exit(1);
	}
	
	//------------------------------------------------//
	
	MggFitter locMggFitter;
	
	YieldFitter locYieldFitter;
	
	InitializeYieldFitter(locYieldFitter, locEtaAna);
	LoadModelType(locYieldFitter, configFileName);
	
	TString locPathName = "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/nll_fit";
	TString fileName = "";
	if(locEtaAna.GetAnalysisOption()==1) {
		double de, e1, e2;
		locEtaAna.GetBeamEnergyBinning(de, e1, e2);
		fileName = Form("%s/output/yield_phase%d_VetoOption%d_%.1fGeV_%.1fGeV.root",
			locPathName.Data(), locEtaAna.GetPhase(),locEtaAna.GetVetoOption(), e1, e2);
	} else {
		fileName = Form("%s/output/yield_phase%d_VetoOption%d_new.root",
			locPathName.Data(), locEtaAna.GetPhase(),locEtaAna.GetVetoOption());
		
		//fileName = "output/yield_phase3_VetoOption6_zqf_theory_omegaOption2.root";
		//fileName = "output/yield_phase3_VetoOption6_zqf_fixed_at_1.root";
		
		//fileName = "/work/halld/home/andrsmit/public/ForIgal/updated_results/loose-bcal-with-sc.root";
	}
	
	if(0)
	{
		if(!locEtaAna.IsMatrixLoaded()) {
			if(locEtaAna.LoadAngularMatrix()) {
				printf("\n\nProblem loading angular matrices.\n\n");
				exit(1);
			}
		}
		if(locEtaAna.CalcAcceptance()) {
			cout << "\n\nProblem getting acceptance.\n\n" << endl;
		}
		
		int loadVal = locEtaAna.LoadLineshapes();
		if(loadVal) {
			printf("\n\nProblem loading MC lineshapes. (error code %d)\n\n", loadVal);
			exit(1);
		}
		if(locEtaAna.LoadDataHistograms()) {
			printf("\n\nProblem accessing data.\n\n");
			exit(1);
		}
		locEtaAna.ExtractAngularYield(locMggFitter, locDrawOption);
		
		locEtaAna.PlotAngularYield();
		locEtaAna.PlotCrossSection();
		
		locEtaAna.PlotEmptyEtaRatio();
		
		locEtaAna.PlotHadronicBkgdFraction();
		locEtaAna.PlotEtaPionFraction();
		
		locEtaAna.PlotBackgrounds();
		locEtaAna.PlotOmegaFitPars();
		
		locEtaAna.WriteROOTFile(fileName);
		//return 0;
		locYieldFitter.SetYield((TH1F*)locEtaAna.GetAngularYield(1));
		locEtaAna.PlotLineshapeShift();
		locEtaAna.PlotQFFraction();
	}
	else if(1) {
		/*
		fileName = Form("%s/output/yield_phase%d_VetoOption%d_new.root",
			locPathName.Data(), locEtaAna.GetPhase(),locEtaAna.GetVetoOption());
		*/
		fileName = Form("%s/output/yield_phase3_VetoOption6_signal3_omega3_rho0.root", locPathName.Data());
		printf("\nGetting yield from file: %s\n", fileName.Data());
		
		TFile *fIn = new TFile(fileName.Data(), "READ");
		
		TH1F *hYield;
		
		if(1) {
			hYield = (TH1F*)fIn->Get("AngularYieldFit");
		}
		else {
			hYield = (TH1F*)fIn->Get("AngularYieldInclusive");
		}
		
		for(int ibin=1; ibin<=hYield->GetXaxis()->GetNbins(); ibin++) {
			double binC = hYield->GetBinContent(ibin);
			double binE = hYield->GetBinError(ibin);
			
			hYield->SetBinError(ibin, 2.0*sqrt(binC));
			/*
			// check for nans:
			if(binE!=binE) {
				hYield->SetBinError(ibin, 1.5*sqrt(binC));
				continue;
			}
			
			if(binE < sqrt(binC)) hYield->SetBinError(ibin, sqrt(binC));
			
			// make sure the error bars are not anomalously large:
			if(binE > 2.5*sqrt(binC)) hYield->SetBinError(ibin, 2.5*sqrt(binC));
			*/
		}
		
		hYield->SetDirectory(0);
		fIn->Close();
		locYieldFitter.SetYield(hYield);
	}
	else {
		// Combined fit across multiple energy bins:
		
		TFile *f1 = new TFile(Form("%s/output/yield_phase3_VetoOption6_8.0GeV_9.0GeV.root", locPathName.Data()), "READ");
		TH1F *hYield1 = (TH1F*)f1->Get("AngularYieldFit")->Clone("Yield1");
		hYield1->SetDirectory(0);
		f1->Close();
		
		TFile *f2 = new TFile(Form("%s/output/yield_phase3_VetoOption6_9.0GeV_10.0GeV.root", locPathName.Data()), "READ");
		TH1F *hYield2 = (TH1F*)f2->Get("AngularYieldFit")->Clone("Yield2");
		hYield2->SetDirectory(0);
		f2->Close();
		
		TFile *f3 = new TFile(Form("%s/output/yield_phase3_VetoOption6_10.0GeV_11.3GeV.root", locPathName.Data()), "READ");
		TH1F *hYield3 = (TH1F*)f3->Get("AngularYieldFit")->Clone("Yield3");
		hYield3->SetDirectory(0);
		f3->Close();
		
		pair<double,double> bin1 = { 8.0,  9.0};
		pair<double,double> bin2 = { 9.0, 10.0};
		pair<double,double> bin3 = {10.0, 11.3};
		
		vector<TH1F*> yields = {hYield1, hYield2, hYield3};
		vector<pair<double,double>> bins = {bin1, bin2, bin3};
		
		for(int ih=0; ih<yields.size(); ih++) {
			for(int ibin=1; ibin<=yields[ih]->GetXaxis()->GetNbins(); ibin++) {
				double binC = yields[ih]->GetBinContent(ibin);
				double binE = yields[ih]->GetBinError(ibin);
				
				// check for nans:
				if(binE!=binE) {
					yields[ih]->SetBinError(ibin, 1.5*sqrt(binC));
					continue;
				}
				
				if(binE < sqrt(binC)) yields[ih]->SetBinError(ibin, sqrt(binC));
			}
		}
		
		//locYieldFitter.SetYield(yields[2]);
		
		locYieldFitter.SetYield(yields, bins);
	}
	
	gStyle->SetStatX(0.50);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(0);
	
	double minFitRange = 0.0, maxFitRange = 2.0;
	locYieldFitter.FitAngularYield(minFitRange, maxFitRange);
	
	if(!batchMode) {
		gSystem->ProcessEvents();
		theApp.Run();
	}
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
		fprintf(stderr,"-h\tPrint this message\n");
		fprintf(stderr,"-p<arg>\tPrimEx-eta phase number (default=1, choose between 1,2,or3)\n");
		fprintf(stderr,"-a<arg>\tAnalysis Option (0: default, 1: invmass matrix, 2: FCAL cuts, 3: BCAL cuts, 4: BEAM cuts, 5: TOF cuts)\n");
		fprintf(stderr,"-c<arg>\tConfiguration file name\n");
		fprintf(stderr,"-b\tRun in batch mode\n");
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
	
	if(ReadFile->GetConfigName("angularRange") != "")
	{
		anaObject.SetAngularRange((double)ReadFile->GetConfig2Par("angularRange")[0], 
				(double)ReadFile->GetConfig2Par("angularRange")[1]);
	}
	
	//----------------------------------------//
	//
	// Analysis and Veto Options:
	//   (mgg histogram and angular matrix names are specified here as well)
	//
	
	if(ReadFile->GetConfigName("analysisOption") != "")
	{
		anaObject.SetAnalysisOption((int)ReadFile->GetConfig1Par("analysisOption")[0]);
		anaObject.SetAnalysisOption_MC((int)ReadFile->GetConfig1Par("analysisOption")[0]);
	}
	if(ReadFile->GetConfigName("vetoOption") != "")
	{
		anaObject.SetVetoOption((int)ReadFile->GetConfig1Par("vetoOption")[0]);
		anaObject.SetVetoOption_MC((int)ReadFile->GetConfig1Par("vetoOption")[0]);
	}
	
	// Default hist name for data+MC:
	anaObject.SetMggHistName(Form("VetoOption%d/mgg_const_cut_veto_%d", anaObject.GetVetoOption(), anaObject.GetVetoOption()));
	anaObject.SetMggHistName_MC(Form("VetoOption%d/mgg_const_cut_veto_%d", anaObject.GetVetoOption(), anaObject.GetVetoOption()));
	
	if(anaObject.GetAnalysisOption() > 0) {
		if(ReadFile->GetConfigName("histName") != ""){
			anaObject.SetMggHistName(ReadFile->GetConfigName("histName"));
			anaObject.SetMggHistName_MC(ReadFile->GetConfigName("histName"));
		}
		else {
			printf("\nNo histogram name specified for analysis option: %d. Exiting.\n", 
				anaObject.GetAnalysisOption());
			exit(1);
		}
	}
	
	// Override analysisOption, vetoOption, and histNames for MC if they were provided:
	
	if(ReadFile->GetConfigName("analysisOption_MC") != "")
	{
		anaObject.SetAnalysisOption_MC((int)ReadFile->GetConfig1Par("analysisOption_MC")[0]);
	}
	if(ReadFile->GetConfigName("vetoOption_MC") != "")
	{
		anaObject.SetVetoOption_MC((int)ReadFile->GetConfig1Par("vetoOption_MC")[0]);
	}
	
	if(anaObject.GetAnalysisOption_MC() > 0) {
		if(ReadFile->GetConfigName("histName_MC") != ""){
			anaObject.SetMggHistName(ReadFile->GetConfigName("histName"));
		}
		else {
			printf("\nNo histogram name specified for MC analysis option: %d. Exiting.\n", 
				anaObject.GetAnalysisOption_MC());
			exit(1);
		}
	}
	
	// Default matrix name:
	anaObject.SetMatrixHistName(Form("AngularMatrix/AngularMatrix_veto_%d", anaObject.GetVetoOption()));
	
	if(anaObject.GetAnalysisOption_MC() > 0) {
		if(ReadFile->GetConfigName("matrixName") != "")
		{
			anaObject.SetMatrixHistName(ReadFile->GetConfigName("matrixName"));
		}
		else {
			printf("\nNo matrix name specified for MC analysis option: %d. Exiting.\n", 
				anaObject.GetAnalysisOption_MC());
			exit(1);
		}
	}
	
	//----------------------------------------//
	// Beam Energy should only be taken from config file if analysisOption=1:
	
	if(anaObject.GetAnalysisOption()==1)
	{
		// beam energy is set according to configuration file:
		if(ReadFile->GetConfigName("beamEnergy") != "")
		{
			anaObject.SetBeamEnergy((double)ReadFile->GetConfig2Par("beamEnergy")[0], 
				(double)ReadFile->GetConfig2Par("beamEnergy")[1]);
		}
	}
	else {
		// beam energy is set according to phase number (hard coded):
		anaObject.SetBeamEnergy();
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
		if(anaObject.GetFitOption(1)>=5) {
			if(ReadFile->GetConfigName("lineshapeOption") != "") {
				anaObject.SetLineshapeOption((int)ReadFile->GetConfig1Par("lineshapeOption")[0]);
			}
		}
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
	if(ReadFile->GetConfigName("rhoFitOption") != "") {
		anaObject.SetFitOption_rho((int)ReadFile->GetConfig1Par("rhoFitOption")[0]);
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
	
	return 0;
}

int LoadModelType(YieldFitter &fitterObject, TString fileName)
{
	// Parse Config File for "ModelType"
	
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
	
	if(ReadFile->GetConfigName("ModelType") != "") {
		fitterObject.SetModel(StringToModelType(ReadFile->GetConfigName("ModelType").Data()));
		
		if(fitterObject.GetModel()==SGEVORKYAN_SIGMA_VAR) {
			if(ReadFile->GetConfigName("ModelParameter_sigma") != "") {
				fitterObject.SetModel_sigmaVer((int)ReadFile->GetConfig1Par("ModelParameter_sigma")[0]);
			}
		}
		else if(fitterObject.GetModel()==SGEVORKYAN_AP_VAR) {
			if(ReadFile->GetConfigName("ModelParameter_ap") != "") {
				fitterObject.SetModel_apVer((int)ReadFile->GetConfig1Par("ModelParameter_ap")[0]);
			}
		}
		else if(fitterObject.GetModel()==SGEVORKYAN_STRONG_RADIUS_VAR) {
			if(ReadFile->GetConfigName("ModelParameter_strongRadiusVer") != "") {
				fitterObject.SetModel_strongRadiusStr(ReadFile->GetConfigName("ModelParameter_strongRadiusVer"));
			}
		}
	}
	return 0;
}

ModelType StringToModelType(std::string str) {
	static const std::unordered_map<std::string, ModelType> lookup = {
		{"AFIX",                 ModelType::AFIX}, 
		{"SGEVORKYAN",           ModelType::SGEVORKYAN},
		{"MIXED_V1",             ModelType::MIXED_V1},
		{"MIXED_V2",             ModelType::MIXED_V2},
		{"SGEVORKYAN_FERMI",     ModelType::SGEVORKYAN_FERMI},
		{"SGEVORKYAN_UPD_V0",    ModelType::SGEVORKYAN_UPD_V0},
		{"SGEVORKYAN_UPD_V1",    ModelType::SGEVORKYAN_UPD_V1},
		{"SGEVORKYAN_UPD_V2",    ModelType::SGEVORKYAN_UPD_V2},
		{"SGEVORKYAN_UPD_V3",    ModelType::SGEVORKYAN_UPD_V3},
		{"SGEVORKYAN_UPD_FERMI", ModelType::SGEVORKYAN_UPD_FERMI},
		{"SGEVORKYAN_SIGMA_VAR", ModelType::SGEVORKYAN_SIGMA_VAR},
		{"SGEVORKYAN_AP_VAR",    ModelType::SGEVORKYAN_AP_VAR},
		{"SGEVORKYAN_STRONG_RADIUS_VAR", ModelType::SGEVORKYAN_STRONG_RADIUS_VAR}
	};
	
	auto it = lookup.find(str);
	if(it != lookup.end())
		return it->second;
	return ModelType::UNKNOWN;
}
