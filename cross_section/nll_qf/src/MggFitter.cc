#include "CrossSection.h"
#include "EtaAnalyzer.h"
#include "MggFitter.h"

//================================================================================================================//
// Constructor:

MggFitter::MggFitter() :
	f_full(nullptr), f_empty(nullptr), f_emptyWide(nullptr),
	h_cohLineshape(nullptr),          h_qfLineshape(nullptr),
	f_cohLineshape(nullptr),          f_qfLineshape(nullptr),
	f_etaLineshape(nullptr),
	h_hadronicBkgdLineshape(nullptr), f_hadronicBkgdLineshape(nullptr), 
	h_etaPionLineshape(nullptr),      f_etaPionLineshape(nullptr), 
	h_omegaLineshape(nullptr),        f_omegaLineshape(nullptr),
	h_rhoLineshape(nullptr)
{
	f_chebyshev = new TF1("chebyshev", "cheb5", 0.0, 1.0);
	
	for(int i=0; i<2; i++) h_full[i] = nullptr;
	h_empty     = nullptr;
	h_emptyWide = nullptr;
}

//================================================================================================================//
// Initialization of parameter arrays:

void MggFitter::InitializeParameterArrays()
{
	/*------------------------------------------------------------------------------------------------------//
	
	This function initializes the private member vectors:
		- 'm_parametersFull': Full list of parameter names
		- 'm_parIndexFull':   Index of each parameter for the full-target fit function within m_parametersFull
		- 'm_parIndexEmpty':  Index of each parameter for the empt-target fit function within m_parametersFull
	
	The latter two vectors are needed to do a simultaneous fit of the empty and full target data.
	
	The full list of parameters includes:
		
		- Signal + Hadronic Background:
			- nominally: 1 normalization and 1 shift parameter for (quasi)elastic signal
			- nominally: 2 normalization (eta+pi and eta+2pi) and 1 shift parameter for hadronic bkgd
		
		- Omega(->pi0+gamma)::
			- nominally: 5 for crystal ball function
		
		- Smooth (non-peaking) Background:
			- nominally: 4 parameters for exponential background
		
		- Empty target Background:
			- nominally: 5x3 for FDC peaking structures (Gaussians), 4 for smooth background
			- 1 for scaling shapes from eta, eta+pi, and omega signals
		
		- Scale factors:
			- empty target flux ratio
			- accidental scaling factor
	
	------------------------------------------------------------------------------------------------------*/
	
	vector<TString> locSignalParameters;
	int nParsSignal = GetSignalParameters(locSignalParameters);
	
	vector<TString> locOmegaParameters;
	int nParsOmega = GetOmegaParameters(locOmegaParameters);
	
	vector<TString> locBkgdParameters;
	int nParsBkgd = GetBkgdParameters(locBkgdParameters);
	
	vector<TString> locEtaPrimeParameters;
	int nParsEtaPrime = GetEtaPrimeParameters(locEtaPrimeParameters);
	
	vector<TString> locEmptyParameters;
	int nParsEmpty = GetEmptyParameters(locEmptyParameters);
	
	// combine all parameters into a single vector:
	
	m_parametersFull.clear();
	m_parametersFull.insert(m_parametersFull.end(), locSignalParameters.begin(),   locSignalParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locOmegaParameters.begin(),    locOmegaParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locBkgdParameters.begin(),     locBkgdParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locEtaPrimeParameters.begin(), locEtaPrimeParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locEmptyParameters.begin(),    locEmptyParameters.end());
	m_parametersFull.push_back("A_{empty}");
	m_parametersFull.push_back("#alpha_{empty}");
	m_parametersFull.push_back("#alpha_{acc}");
	m_parametersFull.push_back("#alpha_{acc,switch}");
	m_parametersFull.push_back("binSize");
	m_parametersFull.push_back("binSize_empty");
	
	int nTotalParameters = (int)m_parametersFull.size();
	
	// Now we need to initial vectors which store the index of each parameter that is used across both functions:
	
	m_parIndexFull.clear();
	m_parIndexEmpty.clear();
	
	for(int ipar=0; ipar<nTotalParameters; ipar++)
	{
		// all but two parameters are used when fitting the full target mgg spectrum:
		if((m_parametersFull[ipar]!="A_{empty}") && (m_parametersFull[ipar]!="binSize_empty")) m_parIndexFull.push_back(ipar);
		
		// all but three parameters are used when fitting the empty target mgg spectrum:
		if((m_parametersFull[ipar]!="#alpha_{empty}") && (m_parametersFull[ipar]!="binSize")
			&& (m_parametersFull[ipar]!="#alpha_{acc}") && (m_parametersFull[ipar]!="#alpha_{acc,switch}")) m_parIndexEmpty.push_back(ipar);
	}
}

struct MggFitter::CombinedNLL {
	
	MggFitter *fitter;
	
	double x1, x2;
	TH1F *hFull, *hEmpty;
	vector<int> indexFull, indexEmpty;
	
	CombinedNLL(MggFitter *f, const vector<int> &index1, const vector<int> &index2) 
		: fitter(f), indexFull(index1), indexEmpty(index2) {
		x1 = fitter->minFitRange;
		x2 = fitter->maxFitRange;
	}
	CombinedNLL(MggFitter *f, double a, double b, const vector<int> &index1, const vector<int> &index2) 
		: fitter(f), x1(a), x2(b), indexFull(index1), indexEmpty(index2) {}
	CombinedNLL() {};
	
	double operator()(const double* p) const {
		vector<double> p1(indexFull.size()), p2(indexEmpty.size());
		for(size_t ipar=0; ipar<indexFull.size(); ipar++)
			p1[ipar] = p[indexFull[ipar]];
		for(size_t jpar=0; jpar<indexEmpty.size(); jpar++)
			p2[jpar] = p[indexEmpty[jpar]];
		
		double nll = 0.0;
		
		//------------------------------------------------------//
		// Likelihood contribution from full target:
		
		int minBin = fitter->h_full[0]->FindBin(x1);
		int maxBin = fitter->h_full[0]->FindBin(x2);
		
		bool fitFull = false, fitEmpty = false;
		switch(fitter->combinedFit) {
			case 0:
				fitFull  = true;
				fitEmpty = false;
				break;
			case 1:
				fitFull  = true;
				fitEmpty = true;
				break;
			case 2:
				fitFull  = false;
				fitEmpty = true;
				break;
		}
		if(fitFull) {
			for(int ibin=minBin; ibin<=maxBin; ibin++) {
				double mgg = fitter->h_full[0]->GetBinCenter(ibin);
				
				int skipBin = 0;
				for(int iexc=0; iexc<(int)fitter->excludeRegions.size(); iexc++) {
					if(fitter->excludeRegions[iexc].first < mgg && mgg < fitter->excludeRegions[iexc].second) {
						skipBin = 1;
						break;
					}
				}
				if(skipBin) continue;
				
				double model_full = fitter->MggFitFunction(&mgg, p1.data());
				if(model_full<=0.0) continue;
				
				nll += (model_full - fitter->h_full[0]->GetBinContent(ibin)*log(model_full));
			}
		}
		
		//------------------------------------------------------//
		// Likelihood contribution from empty target:
		
		if(fitEmpty) 
		{
			minBin = fitter->h_empty->FindBin(x1);
			maxBin = fitter->h_empty->FindBin(x2);
			
			for(int ibin=minBin; ibin<=maxBin; ibin++) {
				double mgg = fitter->h_empty->GetBinCenter(ibin);
				
				int skipBin = 0;
				for(int iexc=0; iexc<(int)fitter->excludeRegions.size(); iexc++) {
					if(fitter->excludeRegions[iexc].first < mgg && mgg < fitter->excludeRegions[iexc].second) {
						skipBin = 1;
						break;
					}
				}
				if(skipBin) continue;
				
				double model_empty = fitter->EmptyMggFitFunction(&mgg, p2.data());
				if(model_empty<=0.0) continue;
				
				nll += (model_empty - fitter->h_emptyWide->GetBinContent(ibin)*log(model_empty));
			}
		}
		
		//------------------------------------------------------//
		// Likelihood contribution from Gaussian constraints:
		
		// add constraint on empty target flux ratio:
		double alpha_beam = 1.0;
		int parInd1 = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), "#alpha_{empty}") - 
			fitter->m_parametersFull.begin());
		if(parInd1 < fitter->m_parametersFull.size()) alpha_beam = p[parInd1];
		double constr1 = 0.5*pow((alpha_beam-1.0)/0.02,2.0);
		
		// add constraint on accidental scaling factor:
		double acc_scale_factor = 0.1;
		int parInd2 = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), "#alpha_{acc}") - 
			fitter->m_parametersFull.begin());
		if(parInd2 < fitter->m_parametersFull.size()) acc_scale_factor = p[parInd2];
		double constr2 = 0.5*pow((acc_scale_factor-1.0)/0.01,2.0);
		
		nll += (constr1 + constr2);
		return nll;
	}
};

//================================================================================================================//
// Main fit routine:

void MggFitter::FitData()
{
	int debug = 0;
	
	/*-----------------------------*/
	// First thing to do is set up vector's storing parameters for full target and empty target fits:
	
	if(debug) printf("\nInitializing parameter arrays...\n");
	InitializeParameterArrays();
	
	// define some useful parameters we will play with later:
	
	// parameter to account for uncertainty in full/empty flux ratio:
	int alpha_empty_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{empty}") - m_parametersFull.begin());
	
	// parameter to account for density of cold gas and target walls:
	int A_empty_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{empty}") - m_parametersFull.begin());
	
	// parameter to account for offset of simulated lineshape compared to dta.
	int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
	
	/*-----------------------------*/
	// Next, initialize the fit functions themselves:
	
	f_etaLineshape = new TF1("f_etaLineshape", this, &MggFitter::EtaLineshape, minFitRange, maxFitRange, 1);
	f_etaLineshape->SetParameter(0, 0.9);
	f_etaLineshape->SetLineColor(kRed);
	f_etaLineshape->SetNpx(1000);
	
	if(debug) printf("Initializing fit functions...\n");
	InitializeFitFunction(&f_full);
	InitializeEmptyFitFunction(&f_empty);
	
	//=======================================================================================================//
	
	/*-----------------------------*/
	// To set up initial guesses for empty bkgd parameters, let's first fit h_emptyWide by itself:
	
	if(debug) printf("Fitting empty target...\n");
	FitEmptyWide();
	if(debug) DumpEmptyFitParameters();
	
	excludeRegions.clear();
	
	/*-----------------------------*/
	// Define two objects that can be used across all subsequent fits:
	
	int locNPars = (int)m_parametersFull.size();
	
	CombinedNLL locNLL(this, m_parIndexFull, m_parIndexEmpty);
	
	double *dummyPars = new double[locNPars];
	for(int ipar=0; ipar<locNPars; ipar++) dummyPars[ipar] = 0.0;
	
	//=======================================================================================================//
	// For the first fit we only let the smooth, in-target background float:
	
	printf("  fit 1: In-target EM background floating...\n");
	
	ROOT::Math::Functor locFCN1(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter1;
	locFitter1.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter1.Config().SetMinimizer("Minuit", "Migrad");
	
	locFitter1.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter1.SetFCN(locFCN1);
	InitializeFitParameters(locFitter1);
	
	ReleaseEmptyParameters(locFitter1);
	
	combinedFit = 2;
	locFitter1.FitFCN();
	
	excludeRegions.push_back({0.50,0.64});
	excludeRegions.push_back({0.68,0.85});
	
	FixEmptyParameters(locFitter1);
	ReleaseBkgdParameters(locFitter1);
	
	combinedFit = 0;
	locFitter1.FitFCN();
	
	auto result1 = locFitter1.Result();
	if(debug) {
		printf("\n\nFit Results (1st fit):\n");
		result1.Print(std::cout);
	}
	//excludeRegions.clear();
	//UpdateFitFunctions(result1);
	//return;
	
	//=======================================================================================================//
	// Now, let the empty target background parameters float if they are within this fitting range:
	
	printf("  fit 2: Empty target floating as well...\n");
	
	ROOT::Math::Functor locFCN2(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter2;
	locFitter2.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter2.Config().SetMinimizer("Minuit", "Migrad");
	
	locFitter2.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter2.SetFCN(locFCN2);
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter2, locFitter1);
	
	// Only do the combined fit in the region where empty background is large:
	if(angle<2.0) {
		
		ReleaseEmptyParameters(locFitter2);
		
		locFitter2.Config().ParSettings(alpha_empty_par).Release();
		locFitter2.Config().ParSettings(alpha_empty_par).SetLimits(0.9, 1.1);
		
		combinedFit = 1;
		locFitter2.FitFCN();
		
		auto result2 = locFitter2.Result();
		if(debug) {
			printf("\n\nFit Results (2nd fit):\n");
			result2.Print(std::cout);
		}
		//excludeRegions.clear();
		//UpdateFitFunctions(result2);
		//return;
	}
	
	//=======================================================================================================//
	// Next, fit the region to the right of the peak, allowing the omega parameters to float:
	
	printf("  fit 3: omega parameters floating...\n");
	
	excludeRegions.clear();
	excludeRegions.push_back({minFitRange,0.65});
	
	locNLL.x1 = 0.65;
	locNLL.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN3(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter3;
	locFitter3.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter3.Config().SetMinimizer("Minuit", "Migrad");
	
	locFitter3.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter3.SetFCN(locFCN3);
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter3, locFitter2);
	
	// Fix parameters that were floating in previous step:
	FixBkgdParameters(locFitter3);
	FixEmptyParameters(locFitter3);
	ReleaseOmegaParameters(locFitter3);
	
	combinedFit = 0;
	locFitter3.FitFCN();
	
	auto result3 = locFitter3.Result();
	if(debug) {
		printf("\n\nFit Results (3rd fit):\n");
		result3.Print(std::cout);
	}
	//excludeRegions.clear();
	//UpdateFitFunctions(result3);
	//return;
	
	//=======================================================================================================//
	// Next, fit the eta mass region with everything else fixed:
	
	printf("  fit 4: eta parameters floating...\n");
	
	excludeRegions.clear();
	excludeRegions.push_back({minFitRange,0.5});
	excludeRegions.push_back({0.6,maxFitRange});
	
	// only fit around eta peak:
	locNLL.x1 = 0.5;
	locNLL.x2 = 0.6;
	
	ROOT::Math::Functor locFCN4(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter4;
	locFitter4.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter4.Config().SetMinimizer("Minuit", "Migrad");
	
	locFitter4.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter4.SetFCN(locFCN4);
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter4, locFitter3);
	
	FixOmegaParameters(locFitter4);
	ReleaseEmptyFDCParameters(locFitter4);
	ReleaseEtaParameters(locFitter4);
	/*
	locFitter4.Config().ParSettings(A_empty_par).Release();
	locFitter4.Config().ParSettings(A_empty_par).SetLimits(0.0, 0.1);
	locFitter4.Config().ParSettings(A_empty_par).SetStepSize(0.005);
	*/
	/*
	locFitter4.Config().ParSettings(offsetPar).Release();
	locFitter4.Config().ParSettings(offsetPar).SetLimits(-0.005, 0.005);
	locFitter4.Config().ParSettings(offsetPar).SetStepSize(0.0001);
	*/
	
	combinedFit = 1;
	locFitter4.FitFCN();
	
	auto result4 = locFitter4.Result();
	if(debug) {
		printf("\n\nFit Results (4th fit):\n");
		result4.Print(std::cout);
	}
	//result4.Print(std::cout);
	//excludeRegions.clear();
	//UpdateFitFunctions(result4);
	//return;
	
	//=======================================================================================================//
	// Fix Eta parameters, and re-fit the empty target background coming from the 2nd FDC package:
	
	printf("  fit 5: eta parameters fixed. Empty parameters floating again...\n");
	
	excludeRegions.clear();
	
	locNLL.x1 = minFitRange;
	locNLL.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN5(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter5;
	locFitter5.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter5.Config().SetMinimizer("Minuit", "Migrad");
	
	locFitter5.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter5.SetFCN(locFCN5);
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter5, locFitter4);
	
	FixEtaParameters(locFitter5);
	ReleaseEmptyParameters(locFitter5);
	
	combinedFit = 2;
	locFitter5.FitFCN();
	
	auto result5 = locFitter5.Result();
	//UpdateFitFunctions(result5);
	//return;
	
	FixEmptyParameters(locFitter5);
	
	//=======================================================================================================//
	// Finally, let everything float:
	
	printf("  fit 6: all parameters floating...\n");
	
	excludeRegions.clear();
	
	locNLL.x1 = minFitRange;
	locNLL.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN6(locNLL, locNPars);
	
	ROOT::Fit::Fitter locFitter6;
	locFitter6.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter6.Config().SetParamsSettings(locNPars, dummyPars);
	
	locFitter6.SetFCN(locFCN6);
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter6, locFitter5);
	
	ReleaseEtaParameters(locFitter6);
	ReleaseBkgdParameters(locFitter6);
	ReleaseOmegaParameters(locFitter6);
	
	combinedFit = 0;
	if(angle<2.0) {
		combinedFit = 1;
		ReleaseEmptyParameters(locFitter6);
	}
	locFitter6.Config().MinimizerOptions().SetMaxFunctionCalls(100000);
	locFitter6.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	
	//locFitter6.Config().SetMinimizer("Minuit", "Simplex");
	//locFitter6.FitFCN();
	locFitter6.Config().SetMinimizer("Minuit", "Migrad");
	locFitter6.FitFCN();
	
	FixEmptyParameters(locFitter6);
	FixOmegaParameters(locFitter6);
	FixBkgdParameters(locFitter6);
	combinedFit = 0;
	
	bool ok = locFitter6.FitFCN();
	
	// which parameters do we want to run Minos errors on:
	unsigned int yieldPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
	unsigned int etaPiPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
	
	//locFitter6.Config().SetMinosErrors({yieldPar, etaPiPar});
	
	/*
	if(!ok) {
		std::cerr << "Migrad failed, retrying with all empty target parameters fixed...\n";
		FixEmptyParameters(locFitter6);
		ok = locFitter6.FitFCN();
	}
	if(!ok) {
		std::cerr << "Migrad failed, retrying with all other background parameters fixed...\n";
		FixOmegaParameters(locFitter6);
		FixBkgdParameters(locFitter6);
		ok = locFitter6.FitFCN();
	}
	*/
	auto result = locFitter6.Result();
	result.Print(std::cout);
	
	UpdateFitFunctions(result);
	/*
	if(ok) {
		// Run HESSE (improves covariance estimates):
		locFitter6.CalculateHessErrors();
		
		locFitter6.CalculateMinosErrors();
		
		// update fit functions with minos errors:
		
		double yieldErrLow  = result.LowerError(yieldPar);
		double yieldErrHigh = result.UpperError(yieldPar);
		double yieldErr     = fabs(yieldErrLow) > fabs(yieldErrHigh) ? fabs(yieldErrLow) : fabs(yieldErrHigh);
		
		double etapiErrLow  = result.LowerError(etaPiPar);
		double etapiErrHigh = result.UpperError(etaPiPar);
		double etapiErr     = fabs(etapiErrLow) > fabs(etapiErrHigh) ? fabs(etapiErrLow) : fabs(etapiErrHigh);
		
		f_full->SetParError(f_full->GetParNumber(m_parametersFull[yieldPar].Data()), yieldErr);
		f_full->SetParError(f_full->GetParNumber(m_parametersFull[etaPiPar].Data()), etapiErr);
	}
	*/
	
	// parameter to account for offset of simulated lineshape compared to dta.
	int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
	double z_qf_fit = result.Parameter(zPar);
	f_etaLineshape->SetParameter(0, z_qf_fit);
	
	return;
}

void MggFitter::UpdateFitFunctions(ROOT::Fit::FitResult result) {
	
	//=======================================================================================================//
	// Copy the parameter errors over to f_full and f_empty:
	
	int locNPars = (int)m_parametersFull.size();
	
	for(int ipar=0; ipar<locNPars; ipar++) {
		TString parName = m_parametersFull[ipar];
		
		double locParVal = result.Parameter(ipar);
		double locParErr = result.ParError(ipar);
		
		// check to see if this parameter is used for f_full:
		if(find(m_parIndexFull.begin(), m_parIndexFull.end(), ipar) != m_parIndexFull.end()) {
			int locParNumber = f_full->GetParNumber(m_parametersFull[ipar].Data());
			f_full->SetParameter(locParNumber, locParVal);
			f_full->SetParError(locParNumber, locParErr);
		}
		
		// Do the same for f_empty (even though I don't think it's needed anywhere):
		if(find(m_parIndexEmpty.begin(), m_parIndexEmpty.end(), ipar) != m_parIndexEmpty.end()) {
			int locParNumber = f_empty->GetParNumber(m_parametersFull[ipar].Data());
			f_empty->SetParameter(locParNumber, locParVal);
			f_empty->SetParError(locParNumber, locParErr);
		}
	}
}

void MggFitter::UpdateFitFunctions(ROOT::Fit::Fitter &fitter) {
	
	//=======================================================================================================//
	// Copy the parameter errors over to f_full and f_empty:
	
	int locNPars = (int)m_parametersFull.size();
	
	for(int ipar=0; ipar<locNPars; ipar++) {
		TString parName = m_parametersFull[ipar];
		
		double locParVal = fitter.Config().ParSettings(ipar).Value();
		
		// check to see if this parameter is used for f_full:
		if(find(m_parIndexFull.begin(), m_parIndexFull.end(), ipar) != m_parIndexFull.end()) {
			int locParNumber = f_full->GetParNumber(m_parametersFull[ipar].Data());
			f_full->SetParameter(locParNumber, locParVal);
		}
		
		// Do the same for f_empty (even though I don't think it's needed anywhere):
		if(find(m_parIndexEmpty.begin(), m_parIndexEmpty.end(), ipar) != m_parIndexEmpty.end()) {
			int locParNumber = f_empty->GetParNumber(m_parametersFull[ipar].Data());
			f_empty->SetParameter(locParNumber, locParVal);
		}
	}
}


void MggFitter::ReleaseBkgdParameters(ROOT::Fit::Fitter &fitter) {
	
	// Let the smooth, in-target background float:
	
	switch(fitOption_bkgd) {
		case 0:
		{
			break;
		}
		case 1:
		{
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Release();
				fitter.Config().ParSettings(pPar).SetLimits(-1.e6, 1.e6);
			}
			break;
		}
		case 2:
		{
			int p0Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0") - m_parametersFull.begin());
			int p1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p1") - m_parametersFull.begin());
			int p2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p2") - m_parametersFull.begin());
			int p3Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p3") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(p0Par).Release();
			fitter.Config().ParSettings(p0Par).SetLimits(0.00, 1.e5);
			
			//fitter.Config().ParSettings(p1Par).Release();
			//fitter.Config().ParSettings(p1Par).SetLimits(0.00, 1.e6);
			
			fitter.Config().ParSettings(p2Par).Release();
			fitter.Config().ParSettings(p2Par).SetLimits(-1.e2, 1.e2);
			
			fitter.Config().ParSettings(p3Par).Release();
			//fitter.Config().ParSettings(p3Par).SetLimits(-1.e2, 1.e2);
			fitter.Config().ParSettings(p3Par).SetLimits(0.0, 1.e2);
			break;
		}
		case 3:
		{
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Release();
				fitter.Config().ParSettings(pPar).SetLimits(-1.e6, 1.e6);
			}
			break;
		}
		case 4:
		{
			// Parametric approximation to dsigma/dM for e+e- pair production
			int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0") - m_parametersFull.begin());
			fitter.Config().ParSettings(pPar).Release();
			fitter.Config().ParSettings(pPar).SetLimits(0., 1.e6);
			break;
		}
	}
}

void MggFitter::FixBkgdParameters(ROOT::Fit::Fitter &fitter)
{
	switch(fitOption_bkgd) {
		case 0:
		{
			break;
		}
		case 1:
		{
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Fix();
			}
			break;
		}
		case 2:
		{
			int p0Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0") - m_parametersFull.begin());
			int p1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p1") - m_parametersFull.begin());
			int p2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p2") - m_parametersFull.begin());
			int p3Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p3") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(p0Par).Fix();
			fitter.Config().ParSettings(p1Par).Fix();
			fitter.Config().ParSettings(p2Par).Fix();
			fitter.Config().ParSettings(p3Par).Fix();
			break;
		}
		case 3:
		{
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Fix();
			}
			break;
		}
		case 4:
		{
			int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0") - m_parametersFull.begin());
			fitter.Config().ParSettings(pPar).Fix();
			break;
		}
	}
}

void MggFitter::ReleaseEmptyParameters(ROOT::Fit::Fitter &fitter) {
	
	// Let the smooth (non-peaking), beamline background float:
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d,empty",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Release();
				fitter.Config().ParSettings(pPar).SetLimits(-1.e7, 1.e7);
			}
			break;
		}
		case 2:
		{
			int p0Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0,empty") - m_parametersFull.begin());
			int p1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p1,empty") - m_parametersFull.begin());
			int p2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p2,empty") - m_parametersFull.begin());
			int p3Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p3,empty") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(p0Par).Release();
			fitter.Config().ParSettings(p0Par).SetLimits(0.00, 1.e6);
			
			//fitter.Config().ParSettings(p1Par).Release();
			//fitter.Config().ParSettings(p1Par).SetLimits(0.00, 1.e6);
			
			fitter.Config().ParSettings(p2Par).Release();
			fitter.Config().ParSettings(p2Par).SetLimits(-1.e3, 1.e3);
			
			fitter.Config().ParSettings(p3Par).Release();
			fitter.Config().ParSettings(p3Par).SetLimits(-1.e3, 1.e3);
			
			break;
		}
		case 3:
		{
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d,empty",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Release();
				fitter.Config().ParSettings(pPar).SetLimits(-1.e7, 1.e7);
			}
			break;
		}
	}
	
	// Next, let the FDC peaking structures float:
	
	ReleaseEmptyFDCParameters(fitter);
}

void MggFitter::ReleaseEmptyFDCParameters(ROOT::Fit::Fitter &fitter) {
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				
				bool skip = false;
				for(int ireg=0; ireg<(int)excludeRegions.size(); ireg++) {
					double locMin = excludeRegions[ireg].first;
					double locMax = excludeRegions[ireg].second;
					if((locMin<m_muFDC[i]) && (m_muFDC[i]<locMax)) skip = true;
				}
				if(skip) continue;
				
				int locNPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("N_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locSigPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#sigma_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				
				fitter.Config().ParSettings(locNPar).Release();
				fitter.Config().ParSettings(locNPar).SetLimits(0.0, 1.e5);
				
				fitter.Config().ParSettings(locMuPar).Release();
				fitter.Config().ParSettings(locMuPar).SetLimits(m_muFDC[i]-0.005, m_muFDC[i]+0.005);
				
				if(0) {
					fitter.Config().ParSettings(locSigPar).Release();
					fitter.Config().ParSettings(locSigPar).SetLimits(0.0075, 0.035);
				}
			}
			break;
		}
		case 2:
		{
			int locNPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{fdc}") - m_parametersFull.begin());
			fitter.Config().ParSettings(locNPar).Release();
			fitter.Config().ParSettings(locNPar).SetLimits(0.0, 1.e5);
			
			for(int i=0; i<m_muFDC.size(); i++) {
				
				bool skip = false;
				for(int ireg=0; ireg<(int)excludeRegions.size(); ireg++) {
					double locMin = excludeRegions[ireg].first;
					double locMax = excludeRegions[ireg].second;
					if((locMin<m_muFDC[i]) && (m_muFDC[i]<locMax)) skip = true;
				}
				if(skip) continue;
				
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#Delta#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				
				fitter.Config().ParSettings(locMuPar).Release();
				fitter.Config().ParSettings(locMuPar).SetLimits(m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				
				bool skip = false;
				for(int ireg=0; ireg<(int)excludeRegions.size(); ireg++) {
					double locMin = excludeRegions[ireg].first;
					double locMax = excludeRegions[ireg].second;
					if((locMin<m_muFDC[i]) && (m_muFDC[i]<locMax)) skip = true;
				}
				if(skip) continue;
				
				int locNPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("N_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#Delta#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				
				fitter.Config().ParSettings(locNPar).Release();
				fitter.Config().ParSettings(locNPar).SetLimits(0.0, 1.e5);
				fitter.Config().ParSettings(locMuPar).Release();
				fitter.Config().ParSettings(locMuPar).SetLimits(m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
	}
}

void MggFitter::FixEmptyParameters(ROOT::Fit::Fitter &fitter) {
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d,empty",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Fix();
			}
			break;
		}
		case 2:
		{
			int p0Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p0,empty") - m_parametersFull.begin());
			int p1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p1,empty") - m_parametersFull.begin());
			int p2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p2,empty") - m_parametersFull.begin());
			int p3Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "p3,empty") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(p0Par).Fix();
			fitter.Config().ParSettings(p1Par).Fix();
			fitter.Config().ParSettings(p2Par).Fix();
			fitter.Config().ParSettings(p3Par).Fix();
			break;
		}
		case 3:
		{
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d,empty",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Fix();
			}
			break;
		}
	}
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("N_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locSigPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#sigma_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				
				fitter.Config().ParSettings(locNPar).Fix();
				fitter.Config().ParSettings(locMuPar).Fix();
				fitter.Config().ParSettings(locSigPar).Fix();
			}
			break;
		}
		case 2:
		{
			int locNPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{fdc}") - m_parametersFull.begin());
			fitter.Config().ParSettings(locNPar).Release();
			fitter.Config().ParSettings(locNPar).SetLimits(0.0, 1.e5);
			
			for(int i=0; i<m_muFDC.size(); i++) {
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#Delta#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				fitter.Config().ParSettings(locMuPar).Fix();
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("N_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				int locMuPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("#Delta#mu_{fdc,%d}",i+1)) 
					- m_parametersFull.begin());
				fitter.Config().ParSettings(locNPar).Fix();
				fitter.Config().ParSettings(locMuPar).Fix();
			}
			break;
		}
	}
	int A_empty_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{empty}") - m_parametersFull.begin());
	fitter.Config().ParSettings(A_empty_par).Fix();
	
	int alpha_empty_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{empty}") - m_parametersFull.begin());
	fitter.Config().ParSettings(alpha_empty_par).Fix();
}


void MggFitter::ReleaseOmegaParameters(ROOT::Fit::Fitter &fitter) {
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
		{
			int NPar     = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int MuPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#mu_{#omega}") - m_parametersFull.begin());
			int SigmaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega}") - m_parametersFull.begin());
			int AlphaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{#omega}") - m_parametersFull.begin());
			int nPar     = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "n_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			fitter.Config().ParSettings(MuPar).Release();
			fitter.Config().ParSettings(MuPar).SetLimits(0.750, 0.800);
			
			fitter.Config().ParSettings(SigmaPar).Release();
			fitter.Config().ParSettings(SigmaPar).SetLimits(0.030, 0.035);
			
			fitter.Config().ParSettings(AlphaPar).Release();
			fitter.Config().ParSettings(AlphaPar).SetLimits(0.200, 9.99);
			
			fitter.Config().ParSettings(nPar).Release();
			fitter.Config().ParSettings(nPar).SetLimits(1.1, 99.99);
			break;
		}
		case 2:
		{
			int NPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			fitter.Config().ParSettings(dMuPar).Release();
			fitter.Config().ParSettings(dMuPar).SetLimits(-0.02, 0.02);
			break;
		}
		case 3:
		{
			int NPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			fitter.Config().ParSettings(dMuPar).Release();
			fitter.Config().ParSettings(dMuPar).SetLimits(-0.02, 0.02);
			break;
		}
		case 4:
		{
			int NomegaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int NrhoPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#rho}") - m_parametersFull.begin());
			int dMuPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NomegaPar).Release();
			fitter.Config().ParSettings(NomegaPar).SetLimits(0.0, 1.e6);
			
			fitter.Config().ParSettings(NrhoPar).Release();
			fitter.Config().ParSettings(NrhoPar).SetLimits(0.0, 1.e5);
			
			fitter.Config().ParSettings(dMuPar).Release();
			fitter.Config().ParSettings(dMuPar).SetLimits(-0.02, 0.02);
			break;
		}
		case 5:
		{
			// allow mean and sigmas to float from lineshape fit:
			
			int NPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int MuPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#mu_{#omega}") - m_parametersFull.begin());
			int Sig1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega,1}") - m_parametersFull.begin());
			int Sig2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega,2}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			fitter.Config().ParSettings(MuPar).Release();
			fitter.Config().ParSettings(MuPar).SetLimits(0.76, 0.80);
			
			fitter.Config().ParSettings(Sig1Par).Release();
			fitter.Config().ParSettings(Sig1Par).SetLimits(0.02, 0.10);
			
			fitter.Config().ParSettings(Sig2Par).Release();
			fitter.Config().ParSettings(Sig2Par).SetLimits(0.02, 0.10);
			break;
		}
	}
}

void MggFitter::FixOmegaParameters(ROOT::Fit::Fitter &fitter) {
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
		{
			int NPar     = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int MuPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#mu_{#omega}") - m_parametersFull.begin());
			int SigmaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega}") - m_parametersFull.begin());
			int AlphaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{#omega}") - m_parametersFull.begin());
			int nPar     = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "n_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(MuPar).Fix();
			fitter.Config().ParSettings(SigmaPar).Fix();
			fitter.Config().ParSettings(AlphaPar).Fix();
			fitter.Config().ParSettings(nPar).Fix();
			break;
		}
		case 2:
		{
			int NPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(dMuPar).Fix();
			break;
		}
		case 3:
		{
			int NPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(dMuPar).Fix();
			break;
		}
		case 4:
		{
			int NomegaPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int NrhoPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#rho}") - m_parametersFull.begin());
			int dMuPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NomegaPar).Fix();
			fitter.Config().ParSettings(NrhoPar).Fix();
			fitter.Config().ParSettings(dMuPar).Fix();
			break;
		}
		case 5:
		{
			// allow mean and sigmas to float from lineshape fit:
			
			int NPar    = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#omega}") - m_parametersFull.begin());
			int MuPar   = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#mu_{#omega}") - m_parametersFull.begin());
			int Sig1Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega,1}") - m_parametersFull.begin());
			int Sig2Par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#sigma_{#omega,2}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(MuPar).Fix();
			fitter.Config().ParSettings(Sig1Par).Fix();
			fitter.Config().ParSettings(Sig2Par).Fix();
			break;
		}
	}
}

void MggFitter::ReleaseEtaParameters(ROOT::Fit::Fitter &fitter) {
	switch(fitOption_signal) {
		case 0:
			break;
		case 1:
		{
			//int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			//fitter.Config().ParSettings(offsetPar).Release();
			//fitter.Config().ParSettings(offsetPar).SetLimits(-0.002, 0.002);
			
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			if(angle<2.5) {
				fitter.Config().ParSettings(zPar).Release();
				fitter.Config().ParSettings(zPar).SetLimits(0.0, 1.0);
			}
			
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NEtaPiPar).Release();
			fitter.Config().ParSettings(NEtaPiPar).SetLimits(0.0, 1.e6);
			
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Release();
				fitter.Config().ParSettings(NHadBgPar).SetLimits(0.0, 10.0);
			}
			if(vetoOption==6) {
				fitter.Config().ParSettings(NHadBgPar).SetValue(0.0);
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
		case 2:
		{
			//int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			//fitter.Config().ParSettings(offsetPar).Release();
			//fitter.Config().ParSettings(offsetPar).SetLimits(-0.002, 0.002);
			
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			fitter.Config().ParSettings(zPar).Release();
			fitter.Config().ParSettings(zPar).SetLimits(0.0, 1.0);
			
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NEtaPiPar).Release();
			fitter.Config().ParSettings(NEtaPiPar).SetLimits(0.0, 1.e6);
			
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Release();
				fitter.Config().ParSettings(NHadBgPar).SetLimits(0.0, 10.0);
			}
			if(vetoOption==6) {
				fitter.Config().ParSettings(NHadBgPar).SetValue(0.0);
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
	}
}

void MggFitter::FixEtaParameters(ROOT::Fit::Fitter &fitter) {
	switch(fitOption_signal) {
		case 0:
			break;
		case 1:
		{
			//int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			//fitter.Config().ParSettings(offsetPar).Release();
			//fitter.Config().ParSettings(offsetPar).SetLimits(-0.002, 0.002);
			
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(zPar).Fix();
			fitter.Config().ParSettings(NEtaPiPar).Fix();
			if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			if(vetoOption==6) {
				//fitter.Config().ParSettings(NHadBgPar).SetValue(0.0);
				//fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
		case 2:
		{
			//int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			//fitter.Config().ParSettings(offsetPar).Release();
			//fitter.Config().ParSettings(offsetPar).SetLimits(-0.002, 0.002);
			
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(zPar).Fix();
			fitter.Config().ParSettings(NEtaPiPar).Fix();
			if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			if(vetoOption==6) {
				//fitter.Config().ParSettings(NHadBgPar).SetValue(0.0);
				//fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
	}
}


//==============================================================//
// Lineshape Fits:

void MggFitter::FitCohLineshape(int drawOption)
{
	if(h_cohLineshape==NULL) return;
	
	h_cohLineshape->SetMarkerStyle(8);
	h_cohLineshape->SetMarkerSize(0.7);
	h_cohLineshape->SetMarkerColor(kBlue);
	h_cohLineshape->SetLineColor(kBlue);
	h_cohLineshape->GetYaxis()->SetTitle(Form("Normalized Counts / %d MeV/c^{2}", (int)(1.e3 * h_cohLineshape->GetBinWidth(1))));
	h_cohLineshape->GetYaxis()->SetTitleSize(0.05);
	h_cohLineshape->GetYaxis()->SetTitleOffset(1.2);
	h_cohLineshape->GetYaxis()->CenterTitle(true);
	h_cohLineshape->GetXaxis()->SetTitleSize(0.05);
	h_cohLineshape->GetXaxis()->SetTitleOffset(1.0);
	h_cohLineshape->GetXaxis()->CenterTitle(true);
	h_cohLineshape->SetTitle("");
	
	// Two Crystal Ball Functions + One Gaussian:
	
	TF1 *locfCoh = new TF1("locfCoh", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
	
	// Initialize parameters:
	
	locfCoh->SetParName(0, "#mu_{1}");
	locfCoh->SetParameter(0, 0.540);
	locfCoh->SetParLimits(0, 0.530,  0.560);
	
	locfCoh->SetParName(1, "#sigma_{1}");
	locfCoh->SetParameter(1, 0.0075);
	locfCoh->SetParLimits(1, 0.004,  0.050);
	
	locfCoh->SetParName(2, "#alpha_{1}");
	locfCoh->SetParameter(2, 1.400);
	locfCoh->SetParLimits(2, 0.200,  9.999);
	
	locfCoh->SetParName(3, "n_{1}");
	locfCoh->SetParameter(3, 9.000);
	locfCoh->SetParLimits(3, 1.100, 49.999);
	
	locfCoh->SetParName(4, "#mu_{2}-#mu_{1}");
	locfCoh->SetParameter(4, 0.015);
	locfCoh->SetParLimits(4,-0.050,  0.050);
	
	locfCoh->SetParName(5, "#sigma_{2}");
	locfCoh->SetParameter(5, 0.010);
	locfCoh->SetParLimits(5, 0.004,  0.050);
	
	locfCoh->SetParName(6, "#alpha_{2}");
	locfCoh->SetParameter(6, 1.400);
	locfCoh->SetParLimits(6, 0.200,  9.999);
	
	locfCoh->SetParName(7, "n_{2}");
	locfCoh->SetParameter(7, 13.50);
	locfCoh->SetParLimits(7, 1.100, 49.999);
	
	locfCoh->SetParName(8, "#mu_{3}-#mu_{1}");
	locfCoh->FixParameter(8, 0.0);
	
	locfCoh->SetParName(9, "#sigma_{3}");
	locfCoh->FixParameter(9, 0.007);
	
	locfCoh->SetParName(10, "frac1");
	locfCoh->SetParameter(10, 0.63);
	locfCoh->SetParLimits(10, 0.0, 1.0);
	
	locfCoh->SetParName(11, "frac2");
	locfCoh->FixParameter(11, 0.0);
	
	locfCoh->SetParName(12, "BinWidth");
	locfCoh->FixParameter(12, h_cohLineshape->GetXaxis()->GetBinWidth(1));
	
	// Force both crystal ball shapes to have the same mean:
	//locfCoh->FixParameter(4, 0.0);
	//locfCoh->FixParameter(8, 0.0);
	
	// Fit:
	h_cohLineshape->Fit(locfCoh, "R0QI");
	
	// Initialize private member, 'f_cohLineshape' based on this fit result:
	f_cohLineshape = new TF1("f_cohLineshape", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
	f_cohLineshape->SetParameters(locfCoh->GetParameters());
	f_cohLineshape->SetLineColor(kCyan);
	
	if(drawOption) {
		TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
		cEtaLS->SetLeftMargin(0.13); cEtaLS->SetRightMargin(0.07);
		cEtaLS->SetBottomMargin(0.13); cEtaLS->SetTopMargin(0.07);
		//gStyle->SetOptFit(0);
		//cEtaLS->SetLogy();
		
		h_cohLineshape->GetXaxis()->SetRangeUser(0.5,0.6);
		h_cohLineshape->Draw();
		f_cohLineshape->SetRange(0.4,0.7);
		f_cohLineshape->SetNpx(1000);
		f_cohLineshape->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = locfCoh->Integral(0.5,0.6) / h_cohLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
		
		TLatex locLat;
		locLat.SetTextFont(42);
		locLat.SetTextSize(0.045);
		
		double locMinAngle = angle - angleWidth;
		double locMaxAngle = angle + angleWidth;
				
		locLat.SetTextColor(kBlue);
		locLat.DrawLatexNDC(0.155, 0.875, 
			Form("#scale[1.0]{%.2f#circ < #theta_{#gamma#gamma} < %.2f#circ}", locMinAngle, locMaxAngle));
		
		locLat.SetTextColor(kBlack);
		//locLat.DrawLatexNDC(0.165, 0.725, Form("#int_{%.2f}^{%.2f}#color[632]{f_{#eta}}dm_{#gamma#gamma} = %.4f", 
		//	minMggCut, maxMggCut, fracAccepted));
		
		/*
		TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9);
		leg1->AddEntry(h_cohLineshape, "Simulation", "PE");
		leg1->AddEntry(f_cohLineshape, "f_{Coh}",    "l");
		leg1->AddEntry(f_qfLineshape, "f_{QF}", "l");
		leg1->AddEntry(f_etaLineshape, "f_{#eta} (z=0.9)", "l");
		leg1->Draw();
		*/
		cEtaLS->Update();
		//cEtaLS->SaveAs(Form("mgg_lineshape_%.2fdeg_%.2fdeg_log.pdf", locMinAngle, locMaxAngle));
		getchar();
		/*
		gStyle->SetOptFit(0);
		cEtaLS->SetLogy();
		cEtaLS->Update();
		//cEtaLS->SaveAs(Form("mgg_lineshape_%.2fdeg_%.2fdeg_log.pdf", locMinAngle, locMaxAngle));
		getchar();
		*/
		delete cEtaLS;
	}
	f_cohLineshape->FixParameter(12, 1.0);
	return;
}


void MggFitter::FitQFLineshape(int drawOption)
{
	if(h_qfLineshape==NULL) return;
	
	h_qfLineshape->SetMarkerStyle(8);
	h_qfLineshape->SetMarkerSize(0.7);
	h_qfLineshape->SetMarkerColor(kBlue);
	h_qfLineshape->SetLineColor(kBlue);
	h_qfLineshape->GetYaxis()->SetTitle(Form("Normalized Counts / %d MeV/c^{2}", (int)(1.e3 * h_qfLineshape->GetBinWidth(1))));
	h_qfLineshape->GetYaxis()->SetTitleSize(0.05);
	h_qfLineshape->GetYaxis()->SetTitleOffset(1.2);
	h_qfLineshape->GetYaxis()->CenterTitle(true);
	h_qfLineshape->GetXaxis()->SetTitleSize(0.05);
	h_qfLineshape->GetXaxis()->SetTitleOffset(1.0);
	h_qfLineshape->GetXaxis()->CenterTitle(true);
	h_qfLineshape->SetTitle("");
	
	// Two Crystal Ball Functions + One Gaussian:
	
	TF1 *locfQF = new TF1("locfQF", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
	
	// Initialize parameters:
	
	locfQF->SetParName(0, "#mu_{1}");
	locfQF->SetParameter(0, 0.540);
	locfQF->SetParLimits(0, 0.530,  0.560);
	
	locfQF->SetParName(1, "#sigma_{1}");
	locfQF->SetParameter(1, 0.0075);
	locfQF->SetParLimits(1, 0.004,  0.050);
	
	locfQF->SetParName(2, "#alpha_{1}");
	locfQF->SetParameter(2, 1.400);
	locfQF->SetParLimits(2, 0.200,  9.999);
	
	locfQF->SetParName(3, "n_{1}");
	locfQF->SetParameter(3, 9.000);
	locfQF->SetParLimits(3, 1.100, 49.999);
	
	locfQF->SetParName(4, "#mu_{2}-#mu_{1}");
	locfQF->SetParameter(4, 0.015);
	locfQF->SetParLimits(4,-0.050,  0.050);
	
	locfQF->SetParName(5, "#sigma_{2}");
	locfQF->SetParameter(5, 0.010);
	locfQF->SetParLimits(5, 0.004,  0.050);
	
	locfQF->SetParName(6, "#alpha_{2}");
	locfQF->SetParameter(6, 1.400);
	locfQF->SetParLimits(6, 0.200,  9.999);
	
	locfQF->SetParName(7, "n_{2}");
	locfQF->SetParameter(7, 13.50);
	locfQF->SetParLimits(7, 1.100, 49.999);
	
	locfQF->SetParName(8, "#mu_{3}-#mu_{1}");
	locfQF->FixParameter(8, 0.0);
	
	locfQF->SetParName(9, "#sigma_{3}");
	locfQF->FixParameter(9, 0.007);
	
	locfQF->SetParName(10, "frac1");
	locfQF->SetParameter(10, 0.63);
	locfQF->SetParLimits(10, 0.0, 1.0);
	
	locfQF->SetParName(11, "frac2");
	locfQF->FixParameter(11, 0.0);
	
	locfQF->SetParName(12, "BinWidth");
	locfQF->FixParameter(12, h_qfLineshape->GetXaxis()->GetBinWidth(1));
	
	// Force both crystal ball shapes to have the same mean:
	//locfQF->FixParameter(4, 0.0);
	//locfQF->FixParameter(8, 0.0);
	
	// Fit:
	h_qfLineshape->Fit(locfQF, "R0QI");
	
	// Initialize private member, 'f_qfLineshape' based on this fit result:
	f_qfLineshape = new TF1("f_qfLineshape", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
	f_qfLineshape->SetParameters(locfQF->GetParameters());
	f_qfLineshape->SetLineColor(kGreen);
	
	if(drawOption) {
		TCanvas *cQFLS = new TCanvas("cQFLS", "cQFLS", 950, 700);
		cQFLS->SetLeftMargin(0.13); cQFLS->SetRightMargin(0.07);
		cQFLS->SetBottomMargin(0.13); cQFLS->SetTopMargin(0.07);
		//gStyle->SetOptFit(0);
		//cQFLS->SetLogy();
		
		h_qfLineshape->GetXaxis()->SetRangeUser(0.5,0.6);
		h_qfLineshape->Draw();
		f_qfLineshape->SetRange(0.4,0.7);
		f_qfLineshape->SetNpx(1000);
		f_qfLineshape->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = locfQF->Integral(0.5,0.6) / h_qfLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
		
		TLatex locLat;
		locLat.SetTextFont(42);
		locLat.SetTextSize(0.045);
		
		double locMinAngle = angle - angleWidth;
		double locMaxAngle = angle + angleWidth;
				
		locLat.SetTextColor(kBlue);
		locLat.DrawLatexNDC(0.155, 0.875, 
			Form("#scale[1.0]{%.2f#circ < #theta_{#gamma#gamma} < %.2f#circ}", locMinAngle, locMaxAngle));
		
		locLat.SetTextColor(kBlack);
		//locLat.DrawLatexNDC(0.165, 0.725, Form("#int_{%.2f}^{%.2f}#color[632]{f_{#eta}}dm_{#gamma#gamma} = %.4f", 
		//	minMggCut, maxMggCut, fracAccepted));
		
		cQFLS->Update();
		getchar();
		
		delete cQFLS;
	}
	f_qfLineshape->FixParameter(12, 1.0);
	return;
}

double MggFitter::EtaLineshape(double *x, double *par) {
	
	double f_coh = f_cohLineshape->Eval(x[0]);
	double f_qf  =  f_qfLineshape->Eval(x[0]);
	return ((1.0-par[0])*f_coh + par[0]*f_qf);
}

void MggFitter::FitHadronicBkgdLineshape(int drawOption)
{
	if(h_hadronicBkgdLineshape==NULL || m_hadronicBkgdYieldBGGEN<1.0) {
		f_hadronicBkgdLineshape = new TF1("f_hadronicBkgdLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
		f_hadronicBkgdLineshape->SetParameters(0.57, 0.015, 1.0, 5.0, 0.02, 0.020, 1.0, 5.0, 0.0, 1.0);
		return;
	}
	
	TF1 *fHadronicBkgd1 = new TF1("fHadronicBkgd1", CrystalBallPDF_flip, 0.5, 0.7, 5);
	
	double    muGuess = h_hadronicBkgdLineshape->GetBinCenter(h_hadronicBkgdLineshape->GetMaximumBin());
	double sigmaGuess = 0.015;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fHadronicBkgd1->SetParameter(0,    muGuess);
	fHadronicBkgd1->SetParameter(1, sigmaGuess);
	fHadronicBkgd1->SetParameter(2, alphaGuess);
	fHadronicBkgd1->SetParameter(3,     nGuess);
	
	fHadronicBkgd1->SetParLimits(0, 0.540,  0.610);
	fHadronicBkgd1->SetParLimits(1, 0.010,  0.050);
	fHadronicBkgd1->SetParLimits(2, 0.200,  9.999);
	fHadronicBkgd1->SetParLimits(3, 1.100, 49.999);
	
	fHadronicBkgd1->FixParameter(4, h_hadronicBkgdLineshape->GetBinWidth(1));
	
	h_hadronicBkgdLineshape->Fit(fHadronicBkgd1, "R0QL");
	
	TF1 *fHadronicBkgd2 = new TF1("fHadronicBkgd2", DoubleCrystalBallPDF_flip, 0.5, 0.75, 10);
	fHadronicBkgd2->SetParameter(0, fHadronicBkgd1->GetParameter(0));
	fHadronicBkgd2->SetParameter(1, fHadronicBkgd1->GetParameter(1));
	fHadronicBkgd2->SetParameter(2, fHadronicBkgd1->GetParameter(2));
	fHadronicBkgd2->SetParameter(3, fHadronicBkgd1->GetParameter(3));
	fHadronicBkgd2->SetParameter(4, 0.02);
	fHadronicBkgd2->SetParameter(5, fHadronicBkgd1->GetParameter(1)*2.0);
	fHadronicBkgd2->SetParameter(6, fHadronicBkgd1->GetParameter(2));
	fHadronicBkgd2->SetParameter(7, fHadronicBkgd1->GetParameter(3));
	fHadronicBkgd2->SetParameter(8, 0.0);
	
	fHadronicBkgd2->SetParLimits(0,  0.540,  0.610);
	fHadronicBkgd2->SetParLimits(1,  0.010,  0.050);
	fHadronicBkgd2->SetParLimits(2,  0.200,  9.999);
	fHadronicBkgd2->SetParLimits(3,  1.100, 49.999);
	fHadronicBkgd2->SetParLimits(4, -0.050,  0.050);
	fHadronicBkgd2->SetParLimits(5,  0.010,  0.050);
	fHadronicBkgd2->SetParLimits(6,  0.200,  9.999);
	fHadronicBkgd2->SetParLimits(7,  1.100, 49.999);
	fHadronicBkgd2->SetParLimits(8,  0.000,  1.000);
	
	fHadronicBkgd2->FixParameter(9, h_hadronicBkgdLineshape->GetBinWidth(1));
	
	h_hadronicBkgdLineshape->Fit(fHadronicBkgd2, "R0Q");
	
	if(drawOption) {
		TCanvas *cHadronicBkgdLS = new TCanvas("cHadronicBkgdLS", "cHadronicBkgdLS", 950, 700);
		//cHadronicBkgdLS->SetLogy();
		h_hadronicBkgdLineshape->Draw();
		fHadronicBkgd2->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = fHadronicBkgd2->Integral(0.5,0.6) / h_hadronicBkgdLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of hadronic bkgd lineshape within mgg cut: %f\n", fracAccepted);
		
		cHadronicBkgdLS->Update();
		getchar();
		delete cHadronicBkgdLS;
	}
	
	f_hadronicBkgdLineshape = new TF1("f_hadronicBkgdLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
	f_hadronicBkgdLineshape->SetParameters(fHadronicBkgd2->GetParameters());
	f_hadronicBkgdLineshape->FixParameter(9, 1.0);
	
	fHadronicBkgd1->Delete();
	fHadronicBkgd2->Delete();
	
	return;
}

void MggFitter::FitEtaPionLineshape(int drawOption)
{
	if(h_etaPionLineshape==NULL || m_etaPionYieldBGGEN<10.0) {
		f_etaPionLineshape = new TF1("f_etaPionLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
		f_etaPionLineshape->SetParameters(0.57, 0.015, 1.0, 5.0, 0.02, 0.020, 1.0, 5.0, 0.0, 1.0);
		return;
	}
	
	TF1 *fEtaPion1 = new TF1("fEtaPion1", CrystalBallPDF_flip, 0.5, 0.7, 5);
	
	double    muGuess = h_etaPionLineshape->GetBinCenter(h_etaPionLineshape->GetMaximumBin());
	double sigmaGuess = 0.015;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fEtaPion1->SetParameter(0,    muGuess);
	fEtaPion1->SetParameter(1, sigmaGuess);
	fEtaPion1->SetParameter(2, alphaGuess);
	fEtaPion1->SetParameter(3,     nGuess);
	
	fEtaPion1->SetParLimits(0, 0.540,  0.610);
	fEtaPion1->SetParLimits(1, 0.010,  0.050);
	fEtaPion1->SetParLimits(2, 0.200,  9.999);
	fEtaPion1->SetParLimits(3, 1.100, 49.999);
	
	fEtaPion1->FixParameter(4, h_etaPionLineshape->GetBinWidth(1));
	
	h_etaPionLineshape->Fit(fEtaPion1, "R0QL");
	
	TF1 *fEtaPion2 = new TF1("fEtaPion2", DoubleCrystalBallPDF_flip, 0.5, 0.75, 10);
	fEtaPion2->SetParameter(0, fEtaPion1->GetParameter(0));
	fEtaPion2->SetParameter(1, fEtaPion1->GetParameter(1));
	fEtaPion2->SetParameter(2, fEtaPion1->GetParameter(2));
	fEtaPion2->SetParameter(3, fEtaPion1->GetParameter(3));
	fEtaPion2->SetParameter(4, 0.02);
	fEtaPion2->SetParameter(5, fEtaPion1->GetParameter(1)*2.0);
	fEtaPion2->SetParameter(6, fEtaPion1->GetParameter(2));
	fEtaPion2->SetParameter(7, fEtaPion1->GetParameter(3));
	fEtaPion2->SetParameter(8, 0.0);
	
	fEtaPion2->SetParLimits(0,  0.540,  0.610);
	fEtaPion2->SetParLimits(1,  0.010,  0.050);
	fEtaPion2->SetParLimits(2,  0.200,  9.999);
	fEtaPion2->SetParLimits(3,  1.100, 49.999);
	fEtaPion2->SetParLimits(4, -0.050,  0.050);
	fEtaPion2->SetParLimits(5,  0.010,  0.050);
	fEtaPion2->SetParLimits(6,  0.200,  9.999);
	fEtaPion2->SetParLimits(7,  1.100, 49.999);
	fEtaPion2->SetParLimits(8,  0.000,  1.000);
	
	fEtaPion2->FixParameter(9, h_etaPionLineshape->GetBinWidth(1));
	
	h_etaPionLineshape->Fit(fEtaPion2, "R0Q");
	
	if(drawOption) {
		TCanvas *cEtaPionLS = new TCanvas("cEtaPionLS", "cEtaPionLS", 950, 700);
		//cEtaPionLS->SetLogy();
		h_etaPionLineshape->Draw();
		fEtaPion2->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = fEtaPion2->Integral(0.5,0.6) / h_etaPionLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of eta+pion lineshape within mgg cut: %f\n", fracAccepted);
		
		cEtaPionLS->Update();
		getchar();
		delete cEtaPionLS;
	}
	
	f_etaPionLineshape = new TF1("f_etaPionLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
	f_etaPionLineshape->SetParameters(fEtaPion2->GetParameters());
	f_etaPionLineshape->FixParameter(9, 1.0);
	
	fEtaPion1->Delete();
	fEtaPion2->Delete();
	
	return;
}

void MggFitter::FitOmegaLineshape(int drawOption)
{
	if(h_omegaLineshape==nullptr) return;
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBallPDF, 0.2, 0.9, 5);
	
	double    muGuess = 0.78;
	double sigmaGuess = 0.025;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	double locMax = 0.0;
	for(int ibin=h_omegaLineshape->FindBin(0.75); ibin<=h_omegaLineshape->FindBin(0.85); ibin++) {
		if(h_omegaLineshape->GetBinContent(ibin) > locMax) {
			locMax  = h_omegaLineshape->GetBinContent(ibin);
			muGuess = h_omegaLineshape->GetBinCenter(ibin);
		}
	}
	
	fOmega1->SetParameter(0,    muGuess);
	fOmega1->SetParameter(1, sigmaGuess);
	fOmega1->SetParameter(2, alphaGuess);
	fOmega1->SetParameter(3,     nGuess);
	
	fOmega1->SetParLimits(0, 0.750,  0.800);
	fOmega1->SetParLimits(1, 0.015,  0.050);
	fOmega1->SetParLimits(2, 0.200,  9.999);
	fOmega1->SetParLimits(3, 1.100, 49.999);
	
	fOmega1->FixParameter(4, h_omegaLineshape->GetBinWidth(1));
	
	h_omegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", DoubleCrystalBallPDF, 0.2, 0.9, 10);
	fOmega2->SetParameter(0, fOmega1->GetParameter(0));
	fOmega2->SetParameter(1, fOmega1->GetParameter(1));
	fOmega2->SetParameter(2, fOmega1->GetParameter(2));
	fOmega2->SetParameter(3, fOmega1->GetParameter(3));
	fOmega2->SetParameter(4, 0.0);
	fOmega2->SetParameter(5, fOmega1->GetParameter(1)*2.0);
	fOmega2->SetParameter(6, fOmega1->GetParameter(2));
	fOmega2->SetParameter(7, fOmega1->GetParameter(3));
	fOmega2->SetParameter(8, 0.0);
	
	fOmega2->SetParName(0, "#mu_{1}");
	fOmega2->SetParName(1, "#sigma_{1}");
	fOmega2->SetParName(2, "#alpha_{1}");
	fOmega2->SetParName(3, "n_{1}");
	fOmega2->SetParName(4, "#mu_{2}-#mu_{1}");
	fOmega2->SetParName(5, "#sigma_{2}");
	fOmega2->SetParName(6, "#alpha_{2}");
	fOmega2->SetParName(7, "n_{2}");
	fOmega2->SetParName(8, "frac");
	
	fOmega2->SetParLimits(0,  0.750,  0.800);
	fOmega2->SetParLimits(1,  0.015,  0.050);
	fOmega2->SetParLimits(2,  0.200,  9.999);
	fOmega2->SetParLimits(3,  1.100, 49.999);
	fOmega2->SetParLimits(4, -0.050,  0.050);
	fOmega2->SetParLimits(5,  0.010,  0.100);
	fOmega2->SetParLimits(6,  0.200,  9.999);
	fOmega2->SetParLimits(7,  1.100, 49.999);
	fOmega2->SetParLimits(8,  0.000,  1.000);
	
	fOmega2->FixParameter(9, h_omegaLineshape->GetBinWidth(1));
	
	h_omegaLineshape->Fit(fOmega2, "R0QL");
	h_omegaLineshape->SetTitle("#gamma+p(n)#rightarrow#omega+p(n)");
	
	if(drawOption) {
		TCanvas *cOmegaLS = new TCanvas("cOmegaLS", "cOmegaLS", 950, 700);
		//cOmegaLS->SetLogy();
		h_omegaLineshape->Draw();
		fOmega2->Draw("same");
		
		cOmegaLS->Update();
		/*
		gStyle->SetOptStat(0);
		cOmegaLS->Modified();
		cOmegaLS->SaveAs("omega_linshape_fit.pdf");
		*/
		getchar();
		delete cOmegaLS;
	}
	
	f_omegaLineshape = new TF1("f_omegaLineshape", DoubleCrystalBallPDF, minFitRange, maxFitRange, 10);
	f_omegaLineshape->SetParameters(fOmega2->GetParameters());
	f_omegaLineshape->FixParameter(9, 1.0);
	
	fOmega1->Delete();
	fOmega2->Delete();
	
	return;
}

void MggFitter::GetYield(double &yield, double &yieldErr, int useFitPars, int subtractHadBkgd, int verbose) {
	
	// if useFitPars==0 (default), extract yield by integrating counts and subtracting bkgd parameters
	// if useFitPars==1, extract yield by integrating signal lineshape
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut+0.0001*binSize);
	int maxMggBin = h_full[0]->FindBin(maxMggCut-0.0001*binSize);
	
	double locMinMggCut = h_full[0]->GetBinCenter(minMggBin) - 0.5*binSize;
	double locMaxMggCut = h_full[0]->GetBinCenter(maxMggBin) + 0.5*binSize;
	
	// check that locMinMggCut and locMaxMggCut align with desired cut range:
	if((fabs(locMinMggCut-minMggCut)>1.e-6) || (fabs(locMaxMggCut-maxMggCut)>1.e-6)) {
		printf("\n\nWarning: Set mgg cut range does not overlap with histogram bin edges.\n");
		printf("  Cut range: %f GeV - %f GeV\n", minMggCut, maxMggCut);
		printf("  Bin edges used in signal integration: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
	}
	//-----------------------------------------------//
	
	TF1 *locfBkgd;
	int nParameters = InitializeFitFunction(&locfBkgd, "locBkgdClone");
	locfBkgd->SetParameters(f_full->GetParameters());
	ZeroSignalPars(locfBkgd, subtractHadBkgd);
	
	if(useFitPars==0) {
		double dataCounts = 0.0;
		double  fitCounts = 0.0;
		for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
			double locData = h_full[0]->GetBinContent(ibin);
			double locBkgd = locfBkgd->Eval(h_full[0]->GetBinCenter(ibin));
			double locFit  = f_full->Eval(h_full[0]->GetBinCenter(ibin));
			fitCounts  += locFit;
			dataCounts += locData;
			yield    += locData - locBkgd;
			yieldErr += pow(h_full[0]->GetBinError(ibin),2.0);
			//printf("  mgg, data, fit: %f %f %f\n", h_full[0]->GetBinCenter(ibin), locData, locFit);
		}
		yieldErr = sqrt(yieldErr);
		
		if(verbose) printf("\nDifference between data and fit: %f sigma\n", (dataCounts-fitCounts)/yieldErr);
	}
	else {
		switch(fitOption_signal) {
			case 1:
			{
				// Eta (histogram) + Eta+Pion Background (histogram) + Other Hadronic Bkgd (histogram):
				
				double lsShift       = f_full->GetParameter("#Delta#mu_{#eta}");
				double locQFFraction = f_full->GetParameter("z_{qf}");
				
				int yieldPar      = f_full->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_full->GetParNumber("N_{#eta#pi}");
				int yieldParBkgd  = f_full->GetParNumber("A_{#eta#pi#pi}");
				
				double N_eta      = f_full->GetParameter(yieldPar);
				double N_etaErr   = f_full->GetParError(yieldPar);
				
				double N_etapi    = f_full->GetParameter(yieldParEtaPi);
				double N_etapiErr = f_full->GetParError(yieldParEtaPi);
				
				double N_bkgd     = f_full->GetParameter(yieldParBkgd) * N_etapi;
				double N_bkgdErr  = sqrt(pow(f_full->GetParError(yieldParBkgd),2.0) + pow(N_etapiErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrectionCoh = h_cohLineshape->Integral(h_cohLineshape->FindBin(minMggCut-lsShift), 
					h_cohLineshape->FindBin(maxMggCut-lsShift));
				double locCorrectionQF  =  h_qfLineshape->Integral(h_qfLineshape->FindBin(minMggCut-lsShift), 
					h_qfLineshape->FindBin(maxMggCut-lsShift));
				
				double locCorrection  = (1.0-locQFFraction)*locCorrectionCoh + locQFFraction*locCorrectionQF;
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift-0.002),
					h_etaPionLineshape->FindBin(maxMggCut-lsShift-0.002));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift-0.002), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift-0.002));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 2:
			{
				// Eta (parameterized lineshape) + Eta+Pion Background (histogram) + Other Hadronic Bkgd (histogram):
				
				double lsShift = f_full->GetParameter("#Delta#mu_{#eta}");
				
				int yieldPar      = f_full->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_full->GetParNumber("N_{#eta#pi}");
				int yieldParBkgd  = f_full->GetParNumber("A_{#eta#pi#pi}");
				
				double N_eta      = f_full->GetParameter(yieldPar);
				double N_etaErr   = f_full->GetParError(yieldPar);
				
				double N_etapi    = f_full->GetParameter(yieldParEtaPi);
				double N_etapiErr = f_full->GetParError(yieldParEtaPi);
				
				double N_bkgd     = f_full->GetParameter(yieldParBkgd) * N_etapi;
				double N_bkgdErr  = sqrt(pow(f_full->GetParError(yieldParBkgd),2.0) + pow(N_etapiErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and 'other' 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift-0.002),
					h_etaPionLineshape->FindBin(maxMggCut-lsShift-0.002));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift-0.002), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift-0.002));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
		}
	}
	
	delete locfBkgd;
	return;
}

void MggFitter::GetEmptyYield(double &yield, double &yieldErr, int excludeNonPeaking) {
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeEmptyFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_empty->GetParameters());
	
	// zero the parameters not associated with peaking structure in eta mass region:
	
	if(excludeNonPeaking){
		
		ZeroOmegaPars(locfEta);
		ZeroEmptyBkgdPars(locfEta);
		ZeroEmptyFDCPars(locfEta);
		ZeroAccPars(locfEta);
		
		// the above function zeroed the peaking structure from the 2nd FDC package. Let's turn it back on:
		locfEta->SetParameter(locfEta->GetParNumber("N_{fdc,2}"), f_empty->GetParameter("N_{fdc,2"));
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_empty->FindBin(minMggCut);
	int maxMggBin = h_empty->FindBin(maxMggCut)-1;
	
	// integrate empty function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locEtas = locfEta->Eval(h_empty->GetBinCenter(ibin));
		yield += locEtas;
	}
	
	// correct this yield by normalization parameter found by fit to full data:
	//double A_empty = f_full->GetParameter("N_{empty}");
	//yield *= (A_empty * m_emptyRatio);
	
	yieldErr = sqrt(yield);
	
	delete locfEta;
	return;
}

void MggFitter::GetOmegaYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	if(fitOption_omega==0) return;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the omega background:
	
	ZeroSignalPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroEmptyPars(locfEta);
	ZeroAccPars(locfEta);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int omegaYieldPar = f_full->GetParNumber("N_{#omega}");
	
	double locRelErr = f_full->GetParError(omegaYieldPar) / f_full->GetParameter(omegaYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}



void MggFitter::GetOmegaFitPars(double &mu, double &muErr, double &sigma, double &sigmaErr,
	double &alpha, double &alphaErr, double &n, double &nErr)
{
	if(fitOption_omega!=1) {
		mu       = 0.0;
		muErr    = 0.0;
		sigma    = 0.0;
		sigmaErr = 0.0;
		alpha    = 0.0;
		alphaErr = 0.0;
		n        = 0.0;
		nErr     = 0.0;
		return;
	}
	
	int muPar    = f_full->GetParNumber("#mu_{#omega}");
	mu           = f_full->GetParameter(muPar);
	muErr        = f_full->GetParError(muPar);
	
	int sigmaPar = f_full->GetParNumber("#sigma_{#omega}");
	sigma        = f_full->GetParameter(sigmaPar);
	sigmaErr     = f_full->GetParError(sigmaPar);
	
	int alphaPar = f_full->GetParNumber("#alpha_{#omega}");
	alpha        = f_full->GetParameter(alphaPar);
	alphaErr     = f_full->GetParError(alphaPar);
	
	int nPar     = f_full->GetParNumber("n_{#omega}");
	n            = f_full->GetParameter(nPar);
	nErr         = f_full->GetParError(nPar);
}

void MggFitter::GetBkgdYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	if(fitOption_bkgd==0) return;
	
	//-----------------------------------------------//
	
	TF1 *locfBkgd;
	int nParameters = InitializeFitFunction(&locfBkgd, "locfBkgd");
	locfBkgd->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the background:
	
	ZeroSignalPars(locfBkgd);
	ZeroOmegaPars(locfBkgd);
	ZeroEtaPrimePars(locfBkgd);
	ZeroEmptyPars(locfBkgd);
	ZeroAccPars(locfBkgd);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfBkgd->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// (placeholder) uncertainty:
	
	yieldErr = sqrt(yield);
	
	delete locfBkgd;
	return;
}

void MggFitter::GetHadronicBkgdYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the hadronic background:
	
	ZeroSignalPars(locfEta, 1);
	ZeroOmegaPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroEmptyPars(locfEta);
	ZeroAccPars(locfEta);
	
	locfEta->SetParameter("A_{#eta#pi}", 0.0);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int bkgdYieldPar = f_full->GetParNumber("A_{#eta#pi#pi}");
	
	double locRelErr = f_full->GetParError(bkgdYieldPar) / f_full->GetParameter(bkgdYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}

void MggFitter::GetEtaPionYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the eta+pion background:
	
	ZeroSignalPars(locfEta, 1);
	ZeroOmegaPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroEmptyPars(locfEta);
	ZeroAccPars(locfEta);
	
	locfEta->SetParameter("A_{#eta#pi#pi}", 0.0);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int etaPionYieldPar = f_full->GetParNumber("N_{#eta#pi}");
	
	double locRelErr = f_full->GetParError(etaPionYieldPar) / f_full->GetParameter(etaPionYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}

void MggFitter::GetLineshapeShift(double &shift, double &shiftErr)
{
	shift    = 0.0;
	shiftErr = 0.0;
	
	int ipar = f_full->GetParNumber("#Delta#mu_{#eta}");
	shift    = f_full->GetParameter(ipar);
	shiftErr = f_full->GetParError(ipar);
	return;
}

void MggFitter::GetQFFraction(double &frac, double &fracErr)
{
	frac    = 0.0;
	fracErr = 0.0;
	
	int ipar = f_full->GetParNumber("z_{qf}");
	frac     = f_full->GetParameter(ipar);
	fracErr  = f_full->GetParError(ipar);
	return;
}

void MggFitter::ZeroSignalPars(TF1 *f1, int excludeHadronicBkgd)
{
	switch(fitOption_signal) {
		case 1:
			if(excludeHadronicBkgd) {
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("N_{#eta#pi}", 0.0);
			}
			break;
		case 2:
			if(excludeHadronicBkgd) {
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("N_{#eta#pi}", 0.0);
			}
			break;
	}
	return;
}

void MggFitter::ZeroHadronicBkgdPars(TF1 *f1)
{
	f1->SetParameter("N_{#eta#pi}", 0.0);
	return;
}

void MggFitter::ZeroOmegaPars(TF1 *f1)
{
	f1->SetParameter("N_{#omega}", 0.0);
	if(fitOption_omega==4) f1->SetParameter("N_{#rho}", 0.0);
	return;
}

void MggFitter::ZeroBkgdPars(TF1 *f1)
{
	switch(fitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
		case 2:
			for(int ipar=0; ipar<4; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
		case 3:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
		case 4:
			f1->SetParameter("p0", 0.0);
			break;
	}
	return;
}

void MggFitter::ZeroEtaPrimePars(TF1 *f1)
{
	if(fitOption_etap==0) return;
	f1->SetParameter("N_{#eta'}", 0.0);
	return;
}

void MggFitter::ZeroEmptyPars(TF1 *f1)
{
	if(fitOption_empty==0) return;
	f1->SetParameter("#alpha_{empty}",0.0);
	return;
}

void MggFitter::ZeroAccPars(TF1 *f1)
{
	f1->SetParameter("#alpha_{acc}",0.0);
	return;
}

void MggFitter::ZeroEmptyBkgdPars(TF1 *f1)
{
	if(fitOption_empty==0) return;
	switch(emptyFitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d,empty",ipar), 0.0);
			}
			break;
		case 2:
			for(int ipar=0; ipar<4; ipar++) {
				f1->SetParameter(Form("p%d,empty",ipar), 0.0);
			}
			break;
		case 3:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d,empty",ipar), 0.0);
			}
			break;
	}
	return;
}

void MggFitter::ZeroEmptyFDCPars(TF1 *f1)
{
	if(fitOption_empty==0) return;
	switch(emptyFitOption_fdc) {
		case 0:
			return;
		case 1:
			for(int ipar=0; ipar<(int)m_muFDC.size(); ipar++) {
				f1->SetParameter(Form("N_{fdc,%d}",ipar+1), 0.0);
			}
			break;
		case 2:
			f1->SetParameter("N_{fdc}",0.0);
			break;
		case 3:
			for(int ipar=0; ipar<(int)m_muFDC.size(); ipar++) {
				f1->SetParameter(Form("N_{fdc,%d}",ipar+1), 0.0);
			}
			break;
	}
	return;
}

void MggFitter::FillPull(TH1F *h_pull)
{
	if(h_pull==NULL) {
		return;
	}
	
	// Just check that the binning of supplied histogram and h_full are consistent:
	if(h_pull->GetXaxis()->GetNbins() != h_full[0]->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of pull histogram.\n" << endl;
	}
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_full[0]->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_full[0]->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_full[0]->GetBinContent(ibin) - f_full->Eval(loc_mgg))/loc_unc);
		h_pull->SetBinError(ibin, 1.0);
	}
	h_pull->GetYaxis()->SetRangeUser(-6.5, 6.5);
	
	h_pull->GetXaxis()->SetTitleSize(0.15);
	h_pull->GetXaxis()->SetLabelSize(0.10);
	h_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma}");
	h_pull->GetYaxis()->SetTitleSize(0.125);
	h_pull->GetYaxis()->SetTitleOffset(0.3);
	h_pull->GetYaxis()->SetLabelSize(0.10);
	
	return;
}

void MggFitter::FillEmptyPull(TH1F *h_pull)
{
	if(h_pull==NULL) {
		return;
	}
	
	// Just check that the binning of supplied histogram and h_full are consistent:
	if(h_empty->GetXaxis()->GetNbins() != h_empty->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of empty pull histogram.\n" << endl;
	}
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_empty->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_empty->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_empty->GetBinContent(ibin) - f_empty->Eval(loc_mgg))/loc_unc);
		h_pull->SetBinError(ibin, 1.0);
	}
	h_pull->GetYaxis()->SetRangeUser(-6.5, 6.5);
	
	h_pull->GetXaxis()->SetTitleSize(0.15);
	h_pull->GetXaxis()->SetLabelSize(0.10);
	h_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma}");
	h_pull->GetYaxis()->SetTitleSize(0.125);
	h_pull->GetYaxis()->SetTitleOffset(0.3);
	h_pull->GetYaxis()->SetLabelSize(0.10);
	
	return;
}

double MggFitter::IntegrateGaussian(double mu, double sigma, double x1, double x2)
{
	// Integrates normalized Gaussian PDF of mean, mu, and width, sigma, between x1 and x2
	// Normalized means that if x1=-inf and x2=+inf, Integral = 1.
	
	double erf1 = TMath::Erf((x1-mu)/(sqrt(2.0)*sigma));
	double erf2 = TMath::Erf((x2-mu)/(sqrt(2.0)*sigma));
	double integral = 0.5 * (erf2 - erf1);
	return integral;
}

void MggFitter::DumpFitParameters()
{
	printf("\n\nFit parameters:\n");
	for(int ipar=0; ipar<f_full->GetNpar()-1; ipar++) {
		printf("  p%d (%s): %f +/- %f\n", ipar, f_full->GetParName(ipar), f_full->GetParameter(ipar), f_full->GetParError(ipar));
	}
	printf("\n\n");
	return;
}

void MggFitter::GetEmptyEtaFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	// get yield of etas from fit to empty target:
	
	double locYieldEmpty, locYieldEmptyErr;
	GetEmptyYield(locYieldEmpty, locYieldEmptyErr);
	
	// get yield of etas from fit to full target:
	
	double locYield, locYieldErr;
	GetYield(locYield, locYieldErr);
	
	fraction    = locYieldEmpty / locYield;
	fractionErr = sqrt(pow(locYieldEmptyErr/locYield,2.0) 
		+ pow(locYieldErr*locYieldEmpty/(locYield*locYield),2.0));
	return;
}

void MggFitter::GetHadronicBkgdFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	int fractionPar;
	
	// For fit options 1 and 2, we don't have the fraction as a direct parameter of the fit function.
	// Instead, we need to calculate it by the ratio of the yields:
	
	double excYield, excYieldErr;
	double bkgYield, bkgYieldErr;
	
	// get yield (data - bkgd fit function) of exclusive part:
	GetYield(excYield, excYieldErr, 0, 1);
	
	// get yield of background:
	GetHadronicBkgdYield(bkgYield, bkgYieldErr);
	
	if((excYield==0) || (excYieldErr==0.0)) return;
	fraction    = bkgYield / excYield;
	fractionErr = sqrt(pow(bkgYieldErr/excYield,2.0) + pow(excYieldErr*bkgYield/(excYield*excYield),2.0));
	return;
}

void MggFitter::GetEtaPionBkgdFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	int fractionPar;
	
	// For fit options 1 and 2, we don't have the fraction as a direct parameter of the fit function.
	// Instead, we need to calculate it by the ratio of the yields:
	
	double excYield, excYieldErr;
	double bkgYield, bkgYieldErr;
	
	// get yield (data - bkgd fit function) of exclusive part:
	GetYield(excYield, excYieldErr, 0, 1);
	
	// get yield of background:
	GetEtaPionYield(bkgYield, bkgYieldErr);
	
	if((excYield==0) || (excYieldErr==0.0)) return;
	fraction    = bkgYield / excYield;
	fractionErr = sqrt(pow(bkgYieldErr/excYield,2.0) + pow(excYieldErr*bkgYield/(excYield*excYield),2.0));
	return;
}

void MggFitter::CheckBinSize(TH1F *h1, TString histTitle)
{
	double locBinWidth = h1->GetXaxis()->GetBinWidth(1);
	if(fabs(locBinWidth-binSize)>1.e-6) {
		printf("\n\nWARNING IN FIT OF %s:\n", histTitle.Data());
		printf("  Bin width does not match MggFitter object (%.3f)\n\n", locBinWidth);
	}
	
	return;
}

int MggFitter::GetSignalParameters(vector<TString> &parameters)
{
	int nParameters = 0;
	switch(fitOption_signal) {
		case 1:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("z_{qf}");
			parameters.push_back("N_{#eta#pi}");
			parameters.push_back("A_{#eta#pi}");
			parameters.push_back("A_{#eta#pi#pi}");
			nParameters += 6;
			break;
		}
		case 2:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("z_{qf}");
			parameters.push_back("N_{#eta#pi}");
			parameters.push_back("A_{#eta#pi}");
			parameters.push_back("A_{#eta#pi#pi}");
			nParameters += 6;
			break;
		}
	}
	return nParameters;
}

int MggFitter::GetOmegaParameters(vector<TString> &parameters)
{
	int nParameters = 0;
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
			parameters.push_back("N_{#omega}");
			parameters.push_back("#mu_{#omega}");
			parameters.push_back("#sigma_{#omega}");
			parameters.push_back("#alpha_{#omega}");
			parameters.push_back("n_{#omega}");
			nParameters += 5;
			break;
		case 2:
			parameters.push_back("N_{#omega}");
			parameters.push_back("#Delta#mu_{#omega}");
			nParameters += 2;
			break;
		case 3:
			parameters.push_back("N_{#omega}");
			parameters.push_back("#Delta#mu_{#omega}");
			nParameters += 2;
			break;
		case 4:
			parameters.push_back("N_{#omega}");
			parameters.push_back("N_{#rho}");
			parameters.push_back("#Delta#mu_{#omega}");
			nParameters += 3;
			break;
		case 5:
			parameters.push_back("N_{#omega}");
			parameters.push_back("#mu_{#omega}");
			parameters.push_back("#sigma_{#omega,1}");
			parameters.push_back("#alpha_{#omega,1}");
			parameters.push_back("#n_{#omega,1}");
			parameters.push_back("#Delta#mu_{#omega}");
			parameters.push_back("#sigma_{#omega,2}");
			parameters.push_back("#alpha_{#omega,2}");
			parameters.push_back("#n_{#omega,2}");
			parameters.push_back("frac_{#omega}");
			nParameters += 10;
			break;
	}
	return nParameters;
}

int MggFitter::GetBkgdParameters(vector<TString> &parameters)
{
	int nParameters = 0;
	switch(fitOption_bkgd) {
		case 0:
			nParameters = 0;
			break;
		case 1:
			nParameters = fitOption_poly + 1;
			break;
		case 2:
			nParameters = 4;
			break;
		case 3:
			nParameters = fitOption_poly + 1;
			break;
		case 4:
			nParameters = 1;
			break;
	}
	for(int ipar=0; ipar<nParameters; ipar++) parameters.push_back(Form("p%d",ipar));
	return nParameters;
}

int MggFitter::GetEtaPrimeParameters(vector<TString> &parameters)
{
	int nParameters = 0;
	if(fitOption_etap==1) {
		parameters.push_back("N_{#eta'}");
		parameters.push_back("#mu_{#eta'}");
		parameters.push_back("#sigma_{#eta'}");
		nParameters = 3;
	}
	return nParameters;
}

int MggFitter::GetEmptyParameters(vector<TString> &parameters)
{
	int nParameters = 0;
	
	// Smooth (non-peaking) background:
	int locNBkgdPars = 0;
	switch(emptyFitOption_bkgd) {
		case 1:
			locNBkgdPars = emptyFitOption_poly + 1;
			break;
		case 2:
			locNBkgdPars = 4;
			break;
		case 3:
			locNBkgdPars = emptyFitOption_poly + 1;
			break;
		case 4:
			locNBkgdPars = 0;
			break;
	}
	nParameters += locNBkgdPars;
	for(int ipar=0; ipar<locNBkgdPars; ipar++) parameters.push_back(Form("p%d,empty",ipar));
	
	// Peaking structures from FDC packages:
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
			for(int i=0; i<m_muFDC.size(); i++) {
				nParameters += 3;
				parameters.push_back(Form("N_{fdc,%d}",i+1));
				parameters.push_back(Form("#mu_{fdc,%d}",i+1));
				parameters.push_back(Form("#sigma_{fdc,%d}",i+1));
			}
			break;
		case 2:
			nParameters += 1;
			parameters.push_back("N_{fdc}");
			for(int i=0; i<m_muFDC.size(); i++) {
				nParameters++;
				parameters.push_back(Form("#Delta#mu_{fdc,%d}",i+1));
			}
			break;
		case 3:
			for(int i=0; i<m_muFDC.size(); i++) {
				nParameters += 2;
				parameters.push_back(Form("N_{fdc,%d}",i+1));
				parameters.push_back(Form("#Delta#mu_{fdc,%d}",i+1));
			}
			break;
		case 4:
			for(int i=0; i<m_muFDC_omega.size(); i++) {
				nParameters += 2;
				parameters.push_back(Form("N_{fdc,#omega,%d}",i+1));
				parameters.push_back(Form("#mu_{fdc,#omega,%d}",i+1));
			}
			for(int i=0; i<m_muFDC_eta.size(); i++) {
				nParameters += 2;
				parameters.push_back(Form("N_{fdc,#eta,%d}",i+1));
				parameters.push_back(Form("#mu_{fdc,#eta,%d}",i+1));
			}
			break;
	}
	
	return nParameters;
}
