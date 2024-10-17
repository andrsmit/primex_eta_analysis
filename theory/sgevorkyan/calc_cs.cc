#include "PrimExCS.h"
#include "FormFactor.h"
#include "TApplication.h"
#include "MyReadConfig.cc"

TH2F *h_absorptionMatrix;

modelParameters_t modelParameters;

int main(int argc, char* argv[])
{
	TApplication theApp("App", &argc, argv);
	
	TString configFileName = "eta_he4.config";
	TString outputFileName = "primex_eta_cs.root";
	TString  inputFileName = "";
	
	double beamEnergy = 10.0;
	double minAngle = 0.0, maxAngle = 4.5;
	double angleBinSize = 0.01;
	
	// parse command line:
	std::cout << argc << std::endl;
	for(int i=1; i<argc; i++) {
		string arg(argv[i]);
		if(arg=="-c") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else { configFileName = argv[i+1]; std::cout << configFileName << std::endl; }
		}
		if(arg=="-o") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else outputFileName = argv[++i]; 
		}
		if(arg=="-i") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else inputFileName = argv[++i];
		}
		if(arg=="-e") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else beamEnergy = atof(argv[++i]);
		}
		if(arg=="-a") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else minAngle = atof(argv[++i]);
		}
		if(arg=="-b") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else maxAngle = atof(argv[++i]);
		}
		if(arg=="-w") {
			if((i+1==argc) || (argv[i+1][0] == '-')) arg = "-h";
			else angleBinSize = atof(argv[++i]);
		}
		if (arg == "-h") {
			cout << endl << " Usage for: " << argv[0] << endl << endl;
			cout << "\t -c <file_name>\t Configuration file name (default is eta_he4.config)" << endl;
			cout << "\t -o <file_name>\t Output ROOT file name (default is primex_eta_cs.root)" << endl;
			cout << "\t -i <file_name>\t  Input ROOT file name to skip form factor calculations [optional]" << endl;
			cout << "\t -e <file_name>\t Beam Energy in GeV (default is 10GeV)" << endl;
			exit(1);
		}
	}
	
	//-------------------------------------------------------------//
	// Initialize PrimExCS object with meson and target specified in config file:
	
	PrimExCS locCS; // default constructor assumes meson=Eta and target=Helium
	
	if(gSystem->AccessPathName(configFileName.Data())) {
		std::cout << "Specified configuration file is not accessible. Quitting program." << std::endl;
		exit(1);
	}
	
	MyReadConfig *ReadFile = new MyReadConfig();
	ReadFile->ReadConfigFile(configFileName);
	
	if(ReadFile->GetConfigName("meson") != "") {
		TString mesonStr = ReadFile->GetConfigName("meson");
		if((mesonStr=="eta") || (mesonStr=="Eta")) locCS.setMeson(Eta);
		else if((mesonStr=="pi0") || (mesonStr=="Pi0")) locCS.setMeson(Pi0);
		else if((mesonStr="etap") || (mesonStr=="Etap")) locCS.setMeson(EtaPrime);
		else {
			std::cout << "Please specify either Eta, Pi0, or Etap." << std::endl;
			exit(1);
		}
	}
	if(ReadFile->GetConfigName("target") != "") {
		TString targetStr = ReadFile->GetConfigName("target");
		if((targetStr=="Helium") || (targetStr=="helium")) locCS.setTarget(Helium);
		else if((targetStr=="C12") || (targetStr=="Carbon") || (targetStr=="carbon")) locCS.setTarget(C12);
		else {
			std::cout << "Unknown target specified. Supported options: Helium, Carbon" << std::endl;
			exit(1);
		}
	}
	
	locCS.setBeamEnergy(beamEnergy);
	locCS.setAngularRange(minAngle, maxAngle);
	locCS.setAngularBinSize(angleBinSize);
	locCS.setOutputFileName(outputFileName);
	
	locCS.initialize();
	
	//-------------------------------------------------------------//
	// Initialize model parameters from configuration file:
	
	modelParameters.targetA = locCS.getTargetA();
	
	if(ReadFile->GetConfigName("targetRadius") != "") {
		modelParameters.targetRadius = ReadFile->GetConfig1Par("targetRadius")[0];
	}
	
	if(ReadFile->GetConfigName("mesonNucleonCS") != "") {
		modelParameters.sigmaby2 = 0.05*ReadFile->GetConfig1Par("mesonNucleonCS")[0];
		// note: the factor of 0.05 is to convert the provided cross section
		//  in units of mb to units of fm2, and divide by 2.
	}
	
	if(ReadFile->GetConfigName("elementaryProductionSlope") != "") {
		modelParameters.slopeAp = ReadFile->GetConfig1Par("elementaryProductionSlope")[0];
	}
	
	if(ReadFile->GetConfigName("elementaryScatteringSlope") != "") {
		modelParameters.slopeAs = ReadFile->GetConfig1Par("elementaryScatteringSlope")[0];
	}
	
	if(ReadFile->GetConfigName("realImaginaryScatteringRatio") != "") {
		modelParameters.alphaIm = ReadFile->GetConfig1Par("realImaginaryScatteringRatio")[0];
	}
	
	if(ReadFile->GetConfigName("shadowingParameter") != "") {
		modelParameters.shadowingParameter = ReadFile->GetConfig1Par("shadowingParameter")[0];
	}
	
	if(ReadFile->GetConfigName("densityModel") != "") {
		TString model = ReadFile->GetConfigName("densityModel");
		if((model!="HO") && (model!="SOG")) {
			std::cout << "Unknown density model specified. Supported options: HO, SOG" << std::endl;
			exit(1);
		}
		else {
			if(model=="HO") {
				modelParameters.densityModel = 0;
			}
			else if(model=="SOG") {
				modelParameters.densityModel = 1;
				InitializeSOGDensityParameters(locCS);
			}
		}
	} else {
		std::cout << "No density model specified. Using the default harmonic oscillator model." << std::endl;
	}
	
	//-------------------------------------------------------------//
	
	printf("\n\n");
	printf("===============================================\n");
	printf("Calculating %s photoproduction cross section on %s with %.1f GeV photons\n", 
		locCS.getMesonName(), locCS.getTargetName(), locCS.getBeamEnergy());
	string density_str = modelParameters.densityModel==0 ? "Harmonic Oscillator model" : "Sum-of-Gaussians";
	printf("  Nuclear density parameterized using %s\n\n", density_str.c_str());
	
	//=================================================================================================//
	// The code below calculates the electromagnetic and strong form factors based on the PrimExCS object
	
	CalculateCoulombFF(locCS);
	locCS.PlotCoulombFF();
	
	CalculateStrongFF(locCS);
	locCS.PlotStrongFF();
	
	locCS.CalculateCrossSection();
	locCS.PlotCrossSection();
	
	getchar();
	
	locCS.WriteToROOTFile();
	
	return 0;
}

int InitializeSOGDensityParameters(PrimExCS csObj) {
	
	if(csObj.getTargetA()==4.0) {
		modelParameters.gammaSOG = sqrt(2.0/3.0);
		modelParameters.densityParametersSOG.clear();
		modelParameters.densityParametersSOG.push_back({0.2, 0.034724});
		modelParameters.densityParametersSOG.push_back({0.6, 0.430761});
		modelParameters.densityParametersSOG.push_back({0.9, 0.203166});
		modelParameters.densityParametersSOG.push_back({1.4, 0.192986});
		modelParameters.densityParametersSOG.push_back({1.9, 0.083866});
		modelParameters.densityParametersSOG.push_back({2.3, 0.033007});
		modelParameters.densityParametersSOG.push_back({2.6, 0.014201});
		modelParameters.densityParametersSOG.push_back({3.1, 0.000000});
		modelParameters.densityParametersSOG.push_back({3.5, 0.006860});
		modelParameters.densityParametersSOG.push_back({4.2, 0.000000});
		modelParameters.densityParametersSOG.push_back({4.9, 0.000438});
		modelParameters.densityParametersSOG.push_back({5.2, 0.000000});
		modelParameters.nParametersSOG = (int)modelParameters.densityParametersSOG.size();
	} else {
		std::cout << "Sum-Of-Gaussians density specified for unsupported target." << std::endl;
		return 1;
	}
	return 0;
}

double Interpolate2D(TH2F *h, double x, double y) {
	
	int xbin = h->GetXaxis()->FindBin(x);
	int ybin = h->GetYaxis()->FindBin(y);
	
	int Nbinsx = h->GetXaxis()->GetNbins();
	int Nbinsy = h->GetYaxis()->GetNbins();
	
	double interpolation = 0.0;
	if(xbin<1) {
		if(ybin<1) {
			interpolation = h->GetBinContent(1,1);
		}
		else if(ybin>Nbinsy) {
			interpolation = h->GetBinContent(1, Nbinsy);
		}
		else {
			interpolation = h->Interpolate(Nbinsx, y);
		}
	}
	else if(xbin>Nbinsx) {
		if(ybin<1) {
			interpolation = h->GetBinContent(Nbinsx,1);
		}
		else if(ybin>h->GetYaxis()->GetNbins()) {
			interpolation = h->GetBinContent(Nbinsx, Nbinsy);
		}
		else {
			interpolation = h->Interpolate(h->GetXaxis()->GetBinCenter(Nbinsx), y);
		}
	}
	else {
		if(ybin<1) {
			interpolation = h->Interpolate(x, h->GetYaxis()->GetBinCenter(1));
		}
		else if(ybin>h->GetYaxis()->GetNbins()) {
			interpolation = h->Interpolate(x, h->GetYaxis()->GetBinCenter(Nbinsy));
		}
		else {
			interpolation = h->Interpolate(x, y);
		}
	}
	
	return interpolation;
}

std::complex<double> IntegrateAbsorption(double b, double z, double exponent) {
	
	if(modelParameters.sigmaby2<=0.0) return {1.0, 0.0};
	
	double integral = 0.;
	if(modelParameters.densityModel==0) {
		
		// Analytical integration (only for HO density model):
		integral = GetGbzIntegral_HO(b, z);
	}
	else {
		
		// Check that absorption matrix is initialized. If not, initialize it:
		if(h_absorptionMatrix == NULL) InitializeAbsorptionMatrix();
		
		integral = Interpolate2D(h_absorptionMatrix, b, z);
	}
	std::complex<double> G_abs(1.0-integral, -modelParameters.alphaIm*integral);
	
	complex<double> F_abs = pow(G_abs, modelParameters.targetA - exponent);
	return F_abs;
}

void InitializeAbsorptionMatrix() {
	
	h_absorptionMatrix = new TH2F ("absorptionMatrix", ";b [fm]; z [fm]", 
		500, 0.0, 5.0*sqrt(2.0*modelParameters.slopeAs), 
		500, -5.0*modelParameters.targetRadius, 5.0*modelParameters.targetRadius);
	
	TF2 *f_Gbz = new TF2("Gbz", Integrand_Gbz, -1.e2, 1.e2, 0.0, 1.e2, 1);
	
	std::cout << "Calculating absorption matrix...";
	for(int xbin=1; xbin<=h_absorptionMatrix->GetXaxis()->GetNbins(); xbin++) {
		double b = h_absorptionMatrix->GetXaxis()->GetBinCenter(xbin);
		std::cout << "  b = " << b << std::endl;
		f_Gbz->SetParameter(0, b);
		for(int ybin=1; ybin<=h_absorptionMatrix->GetYaxis()->GetNbins(); ybin++) {
			double z = h_absorptionMatrix->GetYaxis()->GetBinCenter(ybin);
			h_absorptionMatrix->SetBinContent(xbin, ybin, 
				f_Gbz->Integral(z, 5.0*modelParameters.targetRadius, 0.0, 10.0*sqrt(2.0*modelParameters.slopeAs)));
		}
	}
	std::cout << "done." << std::endl;
	
	f_Gbz->Delete();
	
	return;
}

double Integrand_Gbz(double *x, double *par) {
	
	double zp = x[0];
	double sp = x[1];
	double b  = par[0];
	
	double density  = GetNuclearDensity(sqrt(pow(sp,2.0)+pow(zp,2.0)));
	double besselI0 = ROOT::Math::cyl_bessel_i(0.0, b*sp/modelParameters.slopeAs);
	
	double S_b_sp  = (modelParameters.sigmaby2/modelParameters.slopeAs) *
		sp * exp(-(pow(b,2.0)+pow(sp,2.0))/(2.0*modelParameters.slopeAs)) * besselI0 * density;
	
	return S_b_sp;
}

double GetGbzIntegral_HO(double b, double z) {
	
	// This function returns the analytical expression for G(b,z) 
	// (for the harmonic oscillator density model)
	
	double a_eff = pow(modelParameters.targetRadius,2.0) 
		/ (pow(modelParameters.targetRadius,2.0) + 2.0*modelParameters.slopeAs);
	
	double erfx = std::erf(z/modelParameters.targetRadius);
	double Coef = 0.5*modelParameters.targetRadius*SQRT_PI;
	
	double P1z = Coef*(1.0 + (modelParameters.targetA-4.0)/12.0)*(1.0-erfx) 
		+ ((modelParameters.targetA-4.0)/12.0)*z*exp(-pow(z/modelParameters.targetRadius,2.0));
	double P2z = Coef*(modelParameters.targetA-4.0)/(6.0*pow(modelParameters.targetRadius,2.0))*(1.0-erfx);
	
	double g1 = modelParameters.sigmaby2 * (4.0/modelParameters.targetA) / pow(SQRT_PI*modelParameters.targetRadius,3.0) * a_eff 
		* exp(-pow(b,2.0)/(pow(modelParameters.targetRadius,2.0) + 2.0*modelParameters.slopeAs));
	double g2 = P1z + a_eff*(2.0*modelParameters.slopeAs + a_eff*pow(b,2.0))*P2z;
	
	return g1*g2;
}

double GetNuclearDensity(double r) {
	
	switch(modelParameters.densityModel) {
		case 0 :
			return density_HO(r);
		case 1:
			return density_SOG(r);
		default:
			return 0.0;
	}
}

double density_SOG(double r) {
	
	double density = 0.0;
	for(int ipar=0; ipar<modelParameters.nParametersSOG; ipar++) {
		double Ri = modelParameters.densityParametersSOG[ipar].first;
		double Ai = modelParameters.densityParametersSOG[ipar].second 
			/ (2.0*pow(SQRT_PI*modelParameters.gammaSOG,3.0)*(1.0 + 2.0*pow(Ri/modelParameters.gammaSOG,2.0)));
		double locRho = Ai * (exp(-pow((r-Ri)/modelParameters.gammaSOG,2.0)) 
			+ exp(-pow((r+Ri)/modelParameters.gammaSOG,2.0)));
		density += locRho;
	}
	return density;
}

double density_HO(double r) {
	
	double rho0    = 4.0/pow(SQRT_PI*modelParameters.targetRadius,3.0)/modelParameters.targetA;
	double density = rho0 * (1.0 + pow(r/modelParameters.targetRadius,2.0)*(modelParameters.targetA-4.0)/6.0) 
		* exp(-pow(r/modelParameters.targetRadius,2.0));
	return density;
}
