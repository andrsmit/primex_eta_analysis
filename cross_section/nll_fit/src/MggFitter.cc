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
	
	for(int i=0; i<2; i++) {
		h_full[i]      = nullptr;
		h_empty[i]     = nullptr;
		h_emptyWide[i] = nullptr;
	}
	
	// create a dummy histogram and a simple TF1
	TH1F hDummy("hDummy","dummy",10,0,1);
	TF1 fDummy("fDummy","x",0,1);
	ROOT::Fit::Fitter fitter;
	fitter.Config().MinimizerOptions().SetPrintLevel(0);
	fitter.Config().SetMinimizer("Minuit2", "Migrad");
	
	// a simple fit to trigger first-call initialization
	double pars[1] = {1.0};
	ROOT::Math::Functor fcn([&](const double* p){ return 0.0; }, 1);
	fitter.SetFCN(fcn);
	fitter.Config().SetParamsSettings(1, pars);
	fitter.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	fitter.FitFCN();
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
			- nominally: 1 normalization and 1 weight parameter for (quasi)elastic signal
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
	
	vector<TString> locEMParameters;
	int nParsEM = GetEMParameters(locEMParameters);
	
	vector<TString> locEtaPrimeParameters;
	int nParsEtaPrime = GetEtaPrimeParameters(locEtaPrimeParameters);
	
	vector<TString> locBeamlineParameters;
	int nParsBeamline = GetBeamlineParameters(locBeamlineParameters);
	
	// combine all parameters into a single vector:
	
	m_parametersFull.clear();
	m_parametersFull.insert(m_parametersFull.end(),   locSignalParameters.begin(),   locSignalParameters.end());
	m_parametersFull.insert(m_parametersFull.end(),    locOmegaParameters.begin(),    locOmegaParameters.end());
	m_parametersFull.insert(m_parametersFull.end(),       locEMParameters.begin(),       locEMParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locEtaPrimeParameters.begin(), locEtaPrimeParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locBeamlineParameters.begin(), locBeamlineParameters.end());
	
	m_parametersFull.push_back("A_{empty}");           // (density of cold gas + target walls) / (density of liquid + target walls)
	m_parametersFull.push_back("#alpha_{flux}");       // flux_full / flux_empty
	m_parametersFull.push_back("#alpha_{acc}");        // 1.0 / (# of sidebands)
	m_parametersFull.push_back("#alpha_{acc,switch}"); // just used for drawing purposes
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
		
		// all but two parameters are used when fitting the empty target mgg spectrum:
		if((m_parametersFull[ipar]!="#alpha_{acc,switch}") && (m_parametersFull[ipar]!="binSize")) m_parIndexEmpty.push_back(ipar);
	}
}

struct MggFitter::CombinedNLL {
	
	MggFitter *fitter;
	
	double x1, x2;
	TH1F *hFull, *hEmpty;
	vector<int> indexFull, indexEmpty;
	
	bool debugPrint = false;
	
	CombinedNLL(MggFitter *f, const vector<int> &index1, const vector<int> &index2) 
		: fitter(f), indexFull(index1), indexEmpty(index2) {
		x1 = fitter->minFitRange;
		x2 = fitter->maxFitRange;
	}
	CombinedNLL(MggFitter *f, double a, double b, const vector<int> &index1, const vector<int> &index2) 
		: fitter(f), x1(a), x2(b), indexFull(index1), indexEmpty(index2) {}
	CombinedNLL() {};
	
	double operator()(const double* p) const {
		
		static int callCounter = 0;
		++callCounter;
		if(callCounter==1 && debugPrint) std::cout << "[DEBUG] First call to CombinedNLL" << std::endl;
		
		vector<double> p1(indexFull.size()), p2(indexEmpty.size());
		for(size_t ipar=0; ipar<indexFull.size(); ipar++)
			p1[ipar] = p[indexFull[ipar]];
		for(size_t jpar=0; jpar<indexEmpty.size(); jpar++)
			p2[jpar] = p[indexEmpty[jpar]];
		
		double nll = 0.0;
		
		int minBin = fitter->h_full[0]->FindBin(x1);
		int maxBin = fitter->h_full[0]->FindBin(x2);
		
		bool fitFull = false, fitEmpty = false;
		switch(fitter->combinedFit) {
			case 0:
				fitFull  = true;
				fitEmpty = false;
				//printf("fitting full only\n");
				break;
			case 1:
				fitFull  = true;
				fitEmpty = true;
				break;
			case 2:
				fitFull  = false;
				fitEmpty = true;
				//printf("fitting empty only\n");
				break;
		}
		
		// Get some parameters we'll need for computing the likelihood:
		
		vector<TString> locParNames;
		
		//--------------------------------------------//
		// Signal + Hadronic Background:
		
		locParNames.clear();
		int nParsEta = fitter->GetSignalParameters(locParNames);
		
		vector<double> parsEta(nParsEta+1);
		for(int ipar=0; ipar<nParsEta; ipar++) {
			int locParIndex = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), locParNames[ipar].Data()) - 
				fitter->m_parametersFull.begin());
			parsEta[ipar] = p[locParIndex];
		}
		
		//--------------------------------------------//
		// Omega:
		
		locParNames.clear();
		int nParsOmega = fitter->GetOmegaParameters(locParNames);
		
		vector<double> parsOmega(nParsOmega+1);
		for(int ipar=0; ipar<nParsOmega; ipar++) {
			int locParIndex = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), locParNames[ipar].Data()) - 
				fitter->m_parametersFull.begin());
			parsOmega[ipar] = p[locParIndex];
		}
		
		//--------------------------------------------//
		// e+e- bkgd:
		
		locParNames.clear();
		int nParsEM = fitter->GetEMParameters(locParNames);
		
		vector<double> parsEM(nParsEM+1);
		for(int ipar=0; ipar<nParsEM; ipar++) {
			int locParIndex = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), locParNames[ipar].Data()) - 
				fitter->m_parametersFull.begin());
			parsEM[ipar] = p[locParIndex];
		}
		
		//--------------------------------------------//
		// Beamline:
		
		locParNames.clear();
		int nParsBeamline = fitter->GetBeamlineParameters(locParNames);
		
		vector<double> parsBeamline(nParsBeamline+1);
		for(int ipar=0; ipar<nParsBeamline; ipar++) {
			int locParIndex = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), locParNames[ipar].Data()) - 
				fitter->m_parametersFull.begin());
			parsBeamline[ipar] = p[locParIndex];
		}
		
		//--------------------------------------------//
		// Accidentals:
		
		int par_alpha_acc = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), "#alpha_{acc}") - 
			fitter->m_parametersFull.begin());
		double alpha_acc = p[par_alpha_acc];
		
		//--------------------------------------------//
		// Empty scaling factor (from photon flux):
		
		int par_alpha_flux = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), "#alpha_{flux}") - 
			fitter->m_parametersFull.begin());
		double alpha_flux = p[par_alpha_flux];
		
		//--------------------------------------------//
		// Cold gas scaling factor:
		
		int par_A_empty = (int)(find(fitter->m_parametersFull.begin(), fitter->m_parametersFull.end(), "A_{empty}") - 
			fitter->m_parametersFull.begin());
		double A_empty = p[par_A_empty];
		
		//------------------------------------------------------//
		// Likelihood contribution from full target:
		
		// last parameter of each array should be the histogram bin index:
		
		parsEta[nParsEta]           = fitter->h_full[0]->GetBinWidth(1);
		parsOmega[nParsOmega]       = fitter->h_full[0]->GetBinWidth(1);
		parsEM[nParsEM]             = fitter->h_full[0]->GetBinWidth(1);
		parsBeamline[nParsBeamline] = fitter->h_full[0]->GetBinWidth(1);
		
		if(debugPrint) {
			printf("Eta parameters:\n");
			for(int i=0; i<nParsEta; i++) {
				printf("  %f\n", parsEta[i]);
			}
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
				
				double n_i = fitter->h_full[0]->GetBinContent(ibin); // observed counts in prompt-timing distribution
				double m_i = fitter->h_full[1]->GetBinContent(ibin); // observed counts in accidental-timing distribution
				
				// expected counts (from model) in prompt-timing distribution (minus the accidental contribution):
				
				double s_i = fitter->model_eta(mgg, parsEta.data()) + fitter->model_omega(mgg, parsOmega.data())
					+ fitter->model_em(mgg, parsEM.data()) + alpha_flux * fitter->model_beamline(mgg, parsBeamline.data());
				
				// calculated (profiled) expected counts in the accidental-timing distribution:
				
				double A = alpha_acc * (alpha_acc + 1.0);
				double B = alpha_acc * (s_i - n_i - m_i) + s_i;
				double C = -m_i * s_i;
				double D = pow(B,2.0) - 4.0*A*C;
				if(D<0.0) D = 0.0;
				double eps = 1e-12;
				double a_i_1 = (-B + sqrt(D + eps))/(2.0*A);
				double a_i_2 = (-B - sqrt(D + eps))/(2.0*A);
				
				// choose the root that is positive (if both are positive, choose the larger value):
				
				double a_i;
				if(a_i_1<0.0 && a_i_2<0.0) {
					a_i = 0.0;
				}
				else if(a_i_1>0.0 && a_i_2>0.0) {
					a_i = a_i_1 > a_i_2 ? a_i_1 : a_i_2;
				}
				else {
					a_i = a_i_1 > 0.0 ? a_i_1 : a_i_2;
				}
				
				// likelihood contribution from prompt-timing distribution:
				
				double mu_FP = s_i + alpha_acc * a_i;
				if(mu_FP<=0.0) continue;
				
				nll += (mu_FP - n_i*log(mu_FP));
				
				// likelihood contribution from accidental-timing distribution:
				
				if(a_i>0.0) nll += (a_i - m_i*log(a_i));
				
				if(debugPrint) {
            		// print a summary for each bin
					printf("Bin %d: mgg=%f n_i=%f m_i=%f s_i=%f a_i=%f NLL=%f\n",
						minBin, fitter->h_full[0]->GetBinCenter(minBin),
						fitter->h_full[0]->GetBinContent(minBin),
						fitter->h_full[1]->GetBinContent(minBin),
						s_i, // replace with s_i if you store it
						a_i, // replace with b_i
						nll);
        		}
				
				fitter->h_acc_fit_full->SetBinContent(ibin, a_i);
			}
		}
		
		//------------------------------------------------------//
		// Likelihood contribution from empty target:
		
		// last parameter of each array should be the histogram bin index:
		
		parsEta[nParsEta]           = fitter->h_emptyWide[0]->GetBinWidth(1);
		parsOmega[nParsOmega]       = fitter->h_emptyWide[0]->GetBinWidth(1);
		parsEM[nParsEM]             = fitter->h_emptyWide[0]->GetBinWidth(1);
		parsBeamline[nParsBeamline] = fitter->h_emptyWide[0]->GetBinWidth(1);
		
		if(fitEmpty) 
		{
			minBin = fitter->h_emptyWide[0]->FindBin(x1);
			maxBin = fitter->h_emptyWide[0]->FindBin(x2);
			
			for(int ibin=minBin; ibin<=maxBin; ibin++) {
				double mgg = fitter->h_emptyWide[0]->GetBinCenter(ibin);
				
				int skipBin = 0;
				for(int iexc=0; iexc<(int)fitter->excludeRegions.size(); iexc++) {
					if(fitter->excludeRegions[iexc].first < mgg && mgg < fitter->excludeRegions[iexc].second) {
						skipBin = 1;
						break;
					}
				}
				if(skipBin) continue;
				
				double n_i = fitter->h_emptyWide[0]->GetBinContent(ibin); // observed counts in prompt-timing distribution
				double m_i = fitter->h_emptyWide[1]->GetBinContent(ibin); // observed counts in accidental-timing distribution
				
				// expected counts (from model) in prompt-timing distribution (minus the accidental contribution):
				
				double s_i = (A_empty/alpha_flux) * (fitter->model_eta(mgg, parsEta.data()) + fitter->model_omega(mgg, parsOmega.data())
					+ fitter->model_em(mgg, parsEM.data())) + fitter->model_beamline(mgg, parsBeamline.data());
				
				// calculated (profiled) expected counts in the accidental-timing distribution:
				
				double D = alpha_acc * (alpha_acc + 1.0);
				double E = alpha_acc * (s_i - n_i - m_i) + s_i;
				double F = -m_i * s_i;
				double G = pow(E,2.0) - 4.0*D*F;
				if(G<0.0) G = 0.0;
				double eps = 1e-12;
				double b_i_1 = (-E + sqrt(G+eps))/(2.0*D);
				double b_i_2 = (-E - sqrt(G+eps))/(2.0*D);
				
				// choose the root that is positive (if both are positive, choose the larger value):
				
				double b_i;
				if(b_i_1<0.0 && b_i_2<0.0) {
					b_i = 0.0;
				}
				else if(b_i_1>0.0 && b_i_2>0.0) {
					b_i = b_i_1 > b_i_2 ? b_i_1 : b_i_2;
				}
				else {
					b_i = b_i_1 > 0.0 ? b_i_1 : b_i_2;
				}
				
				// likelihood contribution from prompt-timing distribution:
				
				double mu_EP = s_i + alpha_acc * b_i;
				if(mu_EP<=0.0) continue;
				
				nll += (mu_EP - n_i*log(mu_EP));
				
				// likelihood contribution from accidental-timing distribution:
				
				if(b_i>0.0) nll += (b_i - m_i*log(b_i));
				
				//fitter->h_acc_fit_empty->SetBinContent(ibin, b_i);
				
				if(debugPrint) {
            		// print a summary for each bin
					printf("Bin %d: mgg=%f n_i=%f m_i=%f s_i=%f b_i=%f NLL=%f\n",
						minBin, fitter->h_full[0]->GetBinCenter(minBin),
						fitter->h_full[0]->GetBinContent(minBin),
						fitter->h_full[1]->GetBinContent(minBin),
						s_i, // replace with s_i if you store it
						b_i, // replace with b_i
						nll);
        		}
				
				fitter->h_acc_fit_empty->SetBinContent(ibin, b_i);
			}
		}
		
		//------------------------------------------------------//
		// Likelihood contribution from Gaussian constraints:
		
		// add constraint on empty target flux ratio:
		double constr1 = 0.5*pow((alpha_flux - fitter->m_emptyFluxRatio)/fitter->m_emptyFluxRatioErr, 2.0);
		
		// add constraint on accidental scaling factor:
		double constr2 = 0.5*pow((alpha_acc-0.1)/0.005, 2.0);
		
		nll += (constr1 + constr2);
		
		return nll;
	}
};



//================================================================================================================//
// Main fit routine:

void MggFitter::FitData()
{
	if(TVirtualFitter::GetFitter()) {
		delete TVirtualFitter::GetFitter();
		TVirtualFitter::SetFitter(nullptr);
	}
	
	int debug = 0;
	
	if(f_emptyWide) delete f_emptyWide;
	if(f_empty)     delete f_empty;
	if(f_full)      delete f_full;
	
	if(h_acc_fit_full) {
		delete h_acc_fit_full;
		h_acc_fit_full = (TH1F*)h_full[1]->Clone("acc_fit_full");
	}
	if(h_acc_fit_empty) {
		delete h_acc_fit_empty;
		h_acc_fit_empty = (TH1F*)h_empty[1]->Clone("acc_fit_empty");
	}
	
	excludeRegions.clear();
	m_parIndexFull.clear();
	m_parIndexEmpty.clear();
	
	/*-----------------------------*/
	// First thing to do is set up vector's storing parameters for full target and empty target fits:
	
	if(debug) printf("\nInitializing parameter arrays...\n");
	InitializeParameterArrays();
	
	// define some useful parameters we will play with later:
	
	// parameter to account for uncertainty in full/empty flux ratio:
	int alpha_flux_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{flux}") - m_parametersFull.begin());
	
	// parameter to account for density of cold gas and target walls:
	int A_empty_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{empty}") - m_parametersFull.begin());
	
	// parameter to account for offset of simulated lineshape compared to dta.
	int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
	
	/*-----------------------------*/
	// Next, initialize the fit functions themselves:
	
	f_etaLineshape = new TF1("f_etaLineshape", this, &MggFitter::EtaLineshape, minFitRange, maxFitRange, 1);
	f_etaLineshape->SetParameter(0, 1.0);
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
	vector<double> dummyPars(locNPars, 0.0);
	
	//=======================================================================================================//
	// For the first fit we only let the smooth, in-target background float:
	
	if(debug) printf("  fit 1: In-target EM background floating...\n");
	
	CombinedNLL locNLL1(this, m_parIndexFull, m_parIndexEmpty);
	
	ROOT::Math::Functor locFCN1(locNLL1, locNPars);
	
	ROOT::Fit::Fitter locFitter1;
	locFitter1.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter1.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter1.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter1.SetFCN(locFCN1);
	locFitter1.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	InitializeFitParameters(locFitter1);
	
	vector<double> locTestPars(locNPars, 0.0);
	for(int ipar=0; ipar<locNPars; ipar++) locTestPars[ipar] = locFitter1.Config().ParSettings(ipar).Value();
	locNLL1(locTestPars.data());
	
	ReleaseBeamlineParameters(locFitter1);
	combinedFit = 2;
	locFitter1.FitFCN();
	
	/*
	auto result0 = locFitter1.Result();
	excludeRegions.clear();
	UpdateFitFunctions(result0);
	return;
	*/
	
	excludeRegions.push_back({0.50,0.85});
	
	FixBeamlineParameters(locFitter1);
	ReleaseEMParameters(locFitter1);
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
	
	if(debug) printf("  fit 2: Empty target floating as well...\n");
	
	CombinedNLL locNLL2(this, m_parIndexFull, m_parIndexEmpty);
	
	ROOT::Math::Functor locFCN2(locNLL2, locNPars);
	
	ROOT::Fit::Fitter locFitter2;
	locFitter2.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter2.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter2.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter2.SetFCN(locFCN2);
	locFitter2.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	SetFitParameters(locFitter2, locFitter1); // set from previous fit
	
	// Only do the combined fit in the region where empty background is large:
	//if(angle < 2.0) {
		
		ReleaseBeamlineParameters(locFitter2);
		
		locFitter2.Config().ParSettings(alpha_flux_par).Release();
		locFitter2.Config().ParSettings(alpha_flux_par).SetLimits(0.9*m_emptyFluxRatio, 1.1*m_emptyFluxRatio);
		combinedFit = 1;
		locFitter2.FitFCN();
		
		auto result2 = locFitter2.Result();
		if(debug) {
			printf("\n\nFit Results (2nd fit):\n");
			result2.Print(std::cout);
		}
	//}
	
	//=======================================================================================================//
	// Next, fit the region to the right of the peak, allowing the omega parameters to float:
	
	if(debug) printf("  fit 3: omega parameters floating...\n");
	
	excludeRegions.clear();
	excludeRegions.push_back({minFitRange,0.65});
	
	CombinedNLL locNLL3(this, m_parIndexFull, m_parIndexEmpty);
	locNLL3.x1 = 0.65;
	locNLL3.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN3(locNLL3, locNPars);
	
	ROOT::Fit::Fitter locFitter3;
	locFitter3.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter3.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter3.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter3.SetFCN(locFCN3);
	locFitter3.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	SetFitParameters(locFitter3, locFitter2); // set from previous fit:
	
	// Fix parameters that were floating in previous step:
	FixEMParameters(locFitter3);
	FixBeamlineParameters(locFitter3);
	ReleaseOmegaParameters(locFitter3);
	
	combinedFit = 0;
	locFitter3.FitFCN();
	
	auto result3 = locFitter3.Result();
	if(debug) {
		printf("\n\nFit Results (3rd fit):\n");
		result3.Print(std::cout);
	}
	
	//=======================================================================================================//
	// Next, fit the eta mass region with everything else fixed:
	
	if(debug) printf("  fit 4: eta parameters floating...\n");
	
	excludeRegions.clear();
	excludeRegions.push_back({minFitRange,0.5});
	excludeRegions.push_back({0.6,maxFitRange});
	
	CombinedNLL locNLL4(this, m_parIndexFull, m_parIndexEmpty);
	locNLL4.x1 = 0.5;
	locNLL4.x2 = 0.6;
	
	ROOT::Math::Functor locFCN4(locNLL4, locNPars);
	
	ROOT::Fit::Fitter locFitter4;
	locFitter4.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter4.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter4.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter4.SetFCN(locFCN4);
	locFitter4.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	SetFitParameters(locFitter4, locFitter3); // set from previous fit:
	
	FixOmegaParameters(locFitter4);
	ReleaseFDCParameters(locFitter4);
	ReleaseEtaParameters(locFitter4);
	/*
	locFitter4.Config().ParSettings(A_empty_par).Release();
	locFitter4.Config().ParSettings(A_empty_par).SetLimits(0.02, 0.08);
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
	//excludeRegions.clear();
	//UpdateFitFunctions(result4);
	//return;
	//=======================================================================================================//
	// Fix Eta parameters, and re-fit the empty target background coming from the 2nd FDC package:
	/*
	if(debug) printf("  fit 5: eta parameters fixed. Empty parameters floating again...\n");
	
	excludeRegions.clear();
	
	CombinedNLL locNLL5(this, m_parIndexFull, m_parIndexEmpty);
	locNLL5.x1 = minFitRange;
	locNLL5.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN5(locNLL5, locNPars);
	
	ROOT::Fit::Fitter locFitter5;
	locFitter5.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter5.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter5.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter5.SetFCN(locFCN5);
	locFitter5.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	SetFitParameters(locFitter5, locFitter4); // set from previous fit
	
	FixEtaParameters(locFitter5);
	ReleaseBeamlineParameters(locFitter5);
	
	combinedFit = 1;
	locFitter5.FitFCN();
	
	auto result5 = locFitter5.Result();
	if(debug) {
		printf("\n\nFit Results (5th fit):\n");
		result5.Print(std::cout);
	}
	FixBeamlineParameters(locFitter5);
	*/
	//=======================================================================================================//
	// Finally, let everything float:
	
	if(debug) printf("  fit 6: all parameters floating...\n");
	
	excludeRegions.clear();
	
	CombinedNLL locNLL6(this, m_parIndexFull, m_parIndexEmpty);
	locNLL6.x1 = minFitRange;
	locNLL6.x2 = maxFitRange;
	
	ROOT::Math::Functor locFCN6(locNLL6, locNPars);
	
	ROOT::Fit::Fitter locFitter6;
	locFitter6.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter6.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter6.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter6.SetFCN(locFCN6);
	locFitter6.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	SetFitParameters(locFitter6, locFitter4); // set from previous fit:
	
	ReleaseEtaParameters(locFitter6);
	ReleaseEMParameters(locFitter6);
	ReleaseOmegaParameters(locFitter6);
	
	combinedFit = 0;
	//if(angle<2.0) {
		combinedFit = 1;
		ReleaseBeamlineParameters(locFitter6);
	//}
	
	// which parameters do we want to run Minos errors on:
	unsigned int yieldPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
	unsigned int etaPiPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
	
	//locFitter6.Config().SetMinosErrors({yieldPar, etaPiPar});
	
	bool ok = locFitter6.FitFCN();
	
	/*
	FixBeamlineParameters(locFitter6);
	FixOmegaParameters(locFitter6);
	FixEMParameters(locFitter6);
	combinedFit = 1;
	
	bool ok = locFitter6.FitFCN();
	*/
	
	if(!ok) {
		std::cerr << "Migrad failed, retrying with all empty target parameters fixed...\n";
		FixBeamlineParameters(locFitter6);
		ok = locFitter6.FitFCN();
	}
	if(!ok) {
		std::cerr << "Migrad failed, retrying with all other background parameters fixed...\n";
		FixOmegaParameters(locFitter6);
		FixEMParameters(locFitter6);
		ok = locFitter6.FitFCN();
	}
	
	result = locFitter6.Result();
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
	
	
	double etaYield   = result.Parameter(yieldPar);
	double etaStatErr = result.ParError(yieldPar);
	
	double etaPiYield   = result.Parameter(etaPiPar);
	double etaPiStatErr = result.ParError(etaPiPar);
	
	double covExcInc = result.CovMatrix(yieldPar, etaPiPar);
	double incStatErr = sqrt(pow(etaStatErr,2.0) + pow(etaPiStatErr,2.0) + 2.0*covExcInc);
	
	printf("\n\n\n");
	printf("Eta Yield: %f +/- %f\n", etaYield, etaStatErr);
	printf("Eta+Pi Yield: %f +/- %f\n", etaPiYield, etaPiStatErr);
	printf("Covariance: %f\n", covExcInc);
	printf("Inclusive Yield: %f +/- %f\n", etaYield+etaPiYield, incStatErr);
	printf("\n\n\n");
	
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


void MggFitter::ReleaseEMParameters(ROOT::Fit::Fitter &fitter) {
	
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
			fitter.Config().ParSettings(p0Par).SetLimits(0.00, 1.e7);
			
			//fitter.Config().ParSettings(p1Par).Release();
			//fitter.Config().ParSettings(p1Par).SetLimits(0.00, 1.e6);
			
			//fitter.Config().ParSettings(p2Par).Release();
			//fitter.Config().ParSettings(p2Par).SetLimits(-1.e2, 1.e2);
			
			//fitter.Config().ParSettings(p3Par).Release();
			//fitter.Config().ParSettings(p3Par).SetLimits(-1.e2, 1.e2);
			//fitter.Config().ParSettings(p3Par).SetLimits(-1.e2, 1.e2);
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

void MggFitter::FixEMParameters(ROOT::Fit::Fitter &fitter)
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

void MggFitter::ReleaseBeamlineParameters(ROOT::Fit::Fitter &fitter) {
	
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
	
	ReleaseFDCParameters(fitter);
}

void MggFitter::ReleaseFDCParameters(ROOT::Fit::Fitter &fitter) {
	
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
					if((locMin<m_muFDC[i]) && (m_muFDC[i]<locMax)) {
						skip = true;
						break;
					}
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
				fitter.Config().ParSettings(locMuPar).SetLimits(m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				
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

void MggFitter::FixBeamlineParameters(ROOT::Fit::Fitter &fitter) {
	
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
	
	int alpha_flux_par = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#alpha_{flux}") - m_parametersFull.begin());
	fitter.Config().ParSettings(alpha_flux_par).Fix();
}


void MggFitter::ReleaseOmegaParameters(ROOT::Fit::Fitter &fitter)
{
	switch(fitOption_omega) {
		default:
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
	
	switch(fitOption_rho) {
		default:
			break;
		case 1:
		{
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#rho}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			if((fitOption_omega!=2) && (fitOption_omega!=3)) {
				int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
				fitter.Config().ParSettings(dMuPar).Release();
				fitter.Config().ParSettings(dMuPar).SetLimits(-0.02, 0.02);
			}
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
	
	switch(fitOption_rho) {
		default:
			break;
		case 1:
		{
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#rho}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Fix();
			
			if((fitOption_omega!=2) && (fitOption_omega!=3)) {
				int dMuPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#omega}") - m_parametersFull.begin());
				fitter.Config().ParSettings(dMuPar).Fix();
			}
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
			//if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Release();
				fitter.Config().ParSettings(NHadBgPar).SetLimits(0.0, 10.0);
			//}
			if(vetoOption>6) {
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
			/*
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			fitter.Config().ParSettings(zPar).Release();
			
			double minLimit = 0.0, maxLimit = 1.0;
			if(incFraction_theory>0.1) minLimit = incFraction_theory - 0.1;
			if(incFraction_theory<0.9) maxLimit = incFraction_theory + 0.1;
			
			fitter.Config().ParSettings(zPar).SetLimits(minLimit, maxLimit);
			*/
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NEtaPiPar).Release();
			fitter.Config().ParSettings(NEtaPiPar).SetLimits(0.0, 1.e6);
			
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			//if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Release();
				fitter.Config().ParSettings(NHadBgPar).SetLimits(0.0, 10.0);
			//}
			if(vetoOption>6) {
				fitter.Config().ParSettings(NHadBgPar).SetValue(0.0);
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
		case 3:
		{
			//int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			//fitter.Config().ParSettings(offsetPar).Release();
			//fitter.Config().ParSettings(offsetPar).SetLimits(-0.002, 0.002);
			
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			/*
			int zPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			fitter.Config().ParSettings(zPar).Release();
			
			double minLimit = 0.0, maxLimit = 1.0;
			if(incFraction_theory>0.1) minLimit = incFraction_theory - 0.1;
			if(incFraction_theory<0.9) maxLimit = incFraction_theory + 0.1;
			
			fitter.Config().ParSettings(zPar).SetLimits(minLimit, maxLimit);
			*/
			int fEtaXPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "f_{#etaX}") - m_parametersFull.begin());
			fitter.Config().ParSettings(fEtaXPar).Release();
			fitter.Config().ParSettings(fEtaXPar).SetLimits(0.0, 1.0);
			
			int rEtaXPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "r_{#etaX}") - m_parametersFull.begin());
			fitter.Config().ParSettings(rEtaXPar).Release();
			fitter.Config().ParSettings(rEtaXPar).SetLimits(0.0, 1.0);
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
			int NPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			int zPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(offsetPar).Fix();
			fitter.Config().ParSettings(zPar).Fix();
			fitter.Config().ParSettings(NEtaPiPar).Fix();
			fitter.Config().ParSettings(NHadBgPar).Fix();
			break;
		}
		case 2:
		{
			int NPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			int zPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") - m_parametersFull.begin());
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(offsetPar).Fix();
			fitter.Config().ParSettings(zPar).Fix();
			fitter.Config().ParSettings(NEtaPiPar).Fix();
			fitter.Config().ParSettings(NHadBgPar).Fix();
			break;
		}
		case 3:
		{
			int NPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			int offsetPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "#Delta#mu_{#eta}") - m_parametersFull.begin());
			int zPar      = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "z_{qf}") - m_parametersFull.begin());
			int fEtaXPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "f_{#etaX}") - m_parametersFull.begin());
			int rEtaXPar  = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "r_{#etaX}") - m_parametersFull.begin());
			
			fitter.Config().ParSettings(NPar).Fix();
			fitter.Config().ParSettings(offsetPar).Fix();
			fitter.Config().ParSettings(zPar).Fix();
			fitter.Config().ParSettings(fEtaXPar).Fix();
			fitter.Config().ParSettings(rEtaXPar).Fix();
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
		
		h_cohLineshape->GetXaxis()->SetRangeUser(0.5, 0.6);
		h_qfLineshape->GetXaxis()->SetRangeUser(0.5, 0.6);
		f_cohLineshape->SetParameter(12, h_cohLineshape->GetXaxis()->GetBinWidth(1));
		
		h_qfLineshape->SetLineColor(kRed);
		h_qfLineshape->SetMarkerColor(kRed);
		f_qfLineshape->SetLineColor(kRed);
		f_cohLineshape->SetLineColor(kBlue);
		f_etaLineshape->SetLineColor(kBlack);
		
		f_qfLineshape->SetLineStyle(7);
		f_cohLineshape->SetLineStyle(7);
		
		h_cohLineshape->Draw();
		h_qfLineshape->Draw("same");
		f_qfLineshape->SetRange(0.4,0.7);
		f_qfLineshape->SetNpx(1000);
		f_qfLineshape->Draw("same");
		f_cohLineshape->Draw("same");
		f_etaLineshape->Draw("same");
		
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
		/*
		TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9);
		leg1->AddEntry(h_cohLineshape, "Coh Simulation", "PE");
		leg1->AddEntry(h_qfLineshape, "QF Simulation", "PE");
		leg1->AddEntry(f_cohLineshape, "f_{Coh}",    "l");
		leg1->AddEntry(f_qfLineshape, "f_{QF}", "l");
		leg1->AddEntry(f_etaLineshape, "f_{#eta} (z=0.9)", "l");
		leg1->Draw();
		*/
		cQFLS->Update();
		//cQFLS->SaveAs("lineshape_fit.pdf");
		getchar();
		
		delete cQFLS;
	}
	f_qfLineshape->FixParameter(12, 1.0);
	f_cohLineshape->SetParameter(12, 1.0);
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
	
	if(fitOption_rho==1) {
		f_rhoLineshape = new TF1("f_rhoLineshape", DoubleCrystalBallPDF, minFitRange, maxFitRange, 10);
		f_rhoLineshape->SetParameters(fOmega2->GetParameters());
		f_rhoLineshape->SetParameter(0, fOmega2->GetParameter(0)-0.01);
		f_rhoLineshape->SetParameter(1, 2.75*fOmega2->GetParameter(1));
		f_rhoLineshape->SetParameter(5, 2.75*fOmega2->GetParameter(5));
		f_rhoLineshape->FixParameter(9, 1.0);
	}
	
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
					//yieldErr = locYieldIncErr;
					
					unsigned int etaPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") 
						- m_parametersFull.begin());
					unsigned int etaPiPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") 
						- m_parametersFull.begin());
					
					double etaStatErr   = result.ParError(etaPar);
					double etaPiStatErr = result.ParError(etaPiPar);
					double covExcInc    = result.CovMatrix(etaPar, etaPiPar);
					yieldErr = sqrt(pow(etaStatErr,2.0) + pow(etaPiStatErr,2.0) + 2.0*covExcInc);
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
					//yieldErr = locYieldIncErr;
					
					unsigned int etaPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") 
						- m_parametersFull.begin());
					unsigned int etaPiPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta#pi}") 
						- m_parametersFull.begin());
					
					double etaStatErr   = result.ParError(etaPar);
					double etaPiStatErr = result.ParError(etaPiPar);
					double covExcInc    = result.CovMatrix(etaPar, etaPiPar);
					yieldErr = sqrt(pow(etaStatErr,2.0) + pow(etaPiStatErr,2.0) + 2.0*covExcInc);
				}
				break;
			}
			case 3:
			{
				// Eta (parameterized lineshape) + Eta+Pion Background (histogram) + Other Hadronic Bkgd (histogram):
				
				double lsShift = f_full->GetParameter("#Delta#mu_{#eta}");
				
				int yieldPar = f_full->GetParNumber("N_{#eta}");
				int fEtaXPar = f_full->GetParNumber("f_{#etaX}");
				int rEtaXPar = f_full->GetParNumber("r_{#etaX}");
				
				double N_eta     = f_full->GetParameter(yieldPar);
				double N_etaErr  = f_full->GetParError(yieldPar);
				
				double f_etaX    = f_full->GetParameter(fEtaXPar);
				double f_etaXErr = f_full->GetParError(fEtaXPar);
				
				double r_etaX    = f_full->GetParameter(rEtaXPar);
				double r_etaXErr = f_full->GetParError(rEtaXPar);
				
				//----------------------------------------------------------------------------------------//
				// Overall yields without cut on mgg:
				
				double N_signal      = N_eta * (1.0 - f_etaX);
				double N_signal_err  = sqrt(pow(N_etaErr*(1.0-f_etaX),2.0) + pow(f_etaXErr*N_eta, 2.0));
				
				double N_etapi       = N_eta * f_etaX * r_etaX;
				double N_etapi_err   = sqrt(pow(N_etaErr*f_etaX*r_etaX,2.0) + pow(N_eta*f_etaXErr*r_etaX, 2.0)
					+ pow(N_eta*f_etaX*r_etaXErr,2.0));
				
				double N_etapipi     = N_eta * f_etaX * (1.0-r_etaX);
				double N_etapipi_err = sqrt(pow(N_etaErr*f_etaX*(1.0-r_etaX),2.0) + pow(N_eta*f_etaXErr*(1.0-r_etaX), 2.0)
					+ pow(N_eta*f_etaX*r_etaXErr,2.0));
				
				//----------------------------------------------------------------------------------------//
				// Corrections based on lineshape integrations over mgg cut range:
				
				double locCorrection_signal = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				
				double locCorrection_etapi = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift-0.002),
					h_etaPionLineshape->FindBin(maxMggCut-lsShift-0.002));
				
				double locCorrection_etapipi = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift-0.002), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift-0.002));
				
				//----------------------------------------------------------------------------------------//
				// Corrected yields:
				
				double locYield_signal      = N_signal      * locCorrection_signal;
				double locYield_signal_err  = N_signal_err  * locCorrection_signal;
				
				double locYield_etapi       = N_etapi       * locCorrection_etapi;
				double locYield_etapi_err   = N_etapi_err   * locCorrection_etapi;
				
				double locYield_etapipi     = N_etapipi     * locCorrection_etapipi;
				double locYield_etapipi_err = N_etapipi_err * locCorrection_etapipi;
				
				//----------------------------------------------------------------------------------------//
				// Return the results based on request:
				
				double locYieldInc    = locYield_signal + locYield_etapi + locYield_etapipi;
				double locYieldIncErr = sqrt(pow(locYield_signal_err,2.0) + pow(locYield_etapi_err,2.0) + pow(locYield_etapipi_err,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYield_signal;
					yieldErr = locYield_signal_err;
				} else {
					yield    = locYieldInc;
					//yieldErr = locYieldIncErr;
					
					unsigned int etaPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") 
						- m_parametersFull.begin());
					unsigned int etaXPar = (unsigned int)(find(m_parametersFull.begin(), m_parametersFull.end(), "f_{#etaX}") 
						- m_parametersFull.begin());
					
					double etaStatErr  = result.ParError(etaPar);
					double etaXStatErr = result.ParError(etaXPar);
					double covExcInc   = result.CovMatrix(etaPar, etaXPar);
					yieldErr = sqrt(pow(etaStatErr,2.0) + pow(etaXStatErr,2.0) + 2.0*covExcInc);
				}
				break;
			}
		}
	}
	
	delete locfBkgd;
	return;
}

void MggFitter::GetEmptyYield(double &yield, double &yieldErr) {
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the empty target background:
	
	ZeroSignalPars(locfEta);
	ZeroOmegaPars(locfEta);
	ZeroRhoPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroAccPars(locfEta);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locEtas = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locEtas;
	}
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
	ZeroRhoPars(locfEta);
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

void MggFitter::GetRhoYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	if(fitOption_rho==0) return;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_full->GetParameters());
	
	// zero the parameters not associated with the omega background:
	
	ZeroSignalPars(locfEta);
	ZeroOmegaPars(locfEta);
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
	
	int rhoYieldPar = f_full->GetParNumber("N_{#rho}");
	
	double locRelErr = f_full->GetParError(rhoYieldPar) / f_full->GetParameter(rhoYieldPar);
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
	
	ZeroOmegaPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroEmptyPars(locfEta);
	ZeroAccPars(locfEta);
	
	double relError = 0.0;
	
	switch(fitOption_signal) {
		case 1:
		{
			ZeroSignalPars(locfEta, 1);
			locfEta->SetParameter("A_{#eta#pi}", 0.0);
			int nPar = f_full->GetParNumber("A_{#eta#pi#pi}");
			relError = f_full->GetParError(nPar) / f_full->GetParameter(nPar);
			break;
		}
		case 2:
		{
			ZeroSignalPars(locfEta, 1);
			locfEta->SetParameter("A_{#eta#pi}", 0.0);
			int nPar = f_full->GetParNumber("A_{#eta#pi#pi}");
			relError = f_full->GetParError(nPar) / f_full->GetParameter(nPar);
			break;
		}
		case 3:
		{
			int nPar = f_full->GetParNumber("N_{#eta}");
			int fPar = f_full->GetParNumber("f_{#etaX}");
			int rPar = f_full->GetParNumber("r_{#etaX}");
			
			double N  = f_full->GetParameter(nPar);
			double f  = f_full->GetParameter(fPar);
			double r  = f_full->GetParameter(rPar);
			double Ne = f_full->GetParError(nPar);
			double fe = f_full->GetParError(fPar);
			double re = f_full->GetParError(rPar);
			
			locfEta->SetParameter(nPar, N*f*(1.0-r));
			locfEta->SetParameter(fPar, 1.0);
			locfEta->SetParameter(rPar, 0.0);
			
			relError = sqrt(pow(Ne*f*(1.0-r),2.0) + pow(N*fe*(1.0-r),2.0) + pow(N*f*re,2.0)) / (N*f*r);
			locfEta->SetParError(nPar, relError);
			break;
		}
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	if(relError>2.0) relError = 2.0;
	yieldErr = yield * relError;
	
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
	
	ZeroOmegaPars(locfEta);
	ZeroBkgdPars(locfEta);
	ZeroEtaPrimePars(locfEta);
	ZeroEmptyPars(locfEta);
	ZeroAccPars(locfEta);
	
	double relError = 0.0;
	
	switch(fitOption_signal) {
		case 1:
		{
			ZeroSignalPars(locfEta, 1);
			locfEta->SetParameter("A_{#eta#pi#pi}", 0.0);
			int nPar = f_full->GetParNumber("A_{#eta#pi}");
			relError = f_full->GetParError(nPar) / f_full->GetParameter(nPar);
			break;
		}
		case 2:
		{
			ZeroSignalPars(locfEta, 1);
			locfEta->SetParameter("A_{#eta#pi#pi}", 0.0);
			int nPar = f_full->GetParNumber("A_{#eta#pi}");
			relError = f_full->GetParError(nPar) / f_full->GetParameter(nPar);
			break;
		}
		case 3:
		{
			int nPar = f_full->GetParNumber("N_{#eta}");
			int fPar = f_full->GetParNumber("f_{#etaX}");
			int rPar = f_full->GetParNumber("r_{#etaX}");
			
			double N  = f_full->GetParameter(nPar);
			double f  = f_full->GetParameter(fPar);
			double r  = f_full->GetParameter(rPar);
			double Ne = f_full->GetParError(nPar);
			double fe = f_full->GetParError(fPar);
			double re = f_full->GetParError(rPar);
			
			locfEta->SetParameter(nPar, N*f*r);
			locfEta->SetParameter(fPar, 1.0);
			locfEta->SetParameter(rPar, 1.0);
			
			relError = sqrt(pow(Ne*f*r,2.0) + pow(N*fe*r,2.0) + pow(N*f*re,2.0)) / (N*f*r);
			locfEta->SetParError(nPar, relError);
			break;
		}
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_full[0]->FindBin(minMggCut);
	int maxMggBin = h_full[0]->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full[0]->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	if(relError>2.0) relError = 2.0;
	yieldErr = yield * relError;
	
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
		case 3:
			if(excludeHadronicBkgd) {
				double N = f1->GetParameter("N_{#eta}");
				double f = f1->GetParameter("f_{#etaX}");
				f1->SetParameter("f_{#etaX}", 1.0); // zeros the signal
				f1->SetParameter("N_{#eta}", N*f);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			break;
	}
	return;
}

void MggFitter::ZeroHadronicBkgdPars(TF1 *f1)
{
	switch(fitOption_signal) {
		case 1:
			f1->SetParameter("N_{#eta#pi}", 0.0);
			break;
		case 2:
			f1->SetParameter("N_{#eta#pi}", 0.0);
			break;
		case 3:
		{
			double N = f1->GetParameter("N_{#eta}");
			double f = f1->GetParameter("f_{#etaX}");
			f1->SetParameter("N_{#eta}", N*(1.0-f));
			f1->SetParameter("f_{#etaX}", 0.0);
			break;
		}
	}
	return;
}

void MggFitter::ZeroOmegaPars(TF1 *f1)
{
	f1->SetParameter("N_{#omega}", 0.0);
	return;
}

void MggFitter::ZeroRhoPars(TF1 *f1)
{
	f1->SetParameter("N_{#rho}", 0.0);
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
	f1->SetParameter("#alpha_{flux}",0.0);
	return;
}

void MggFitter::ZeroAccPars(TF1 *f1)
{
	f1->SetParameter("#alpha_{acc}",0.0);
	return;
}

void MggFitter::ZeroEmptyBkgdPars(TF1 *f1)
{
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
	if(h_empty[0]->GetXaxis()->GetNbins() != h_empty[0]->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of empty pull histogram.\n" << endl;
	}
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_empty[0]->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_empty[0]->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_empty[0]->GetBinContent(ibin) - f_empty->Eval(loc_mgg))/loc_unc);
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
	
	if(fitOption_signal==3)
	{
		int fEtaXPar = f_full->GetParNumber("f_{#etaX}");
		int rEtaXPar = f_full->GetParNumber("r_{#etaX}");
		
		double fEtaX     = f_full->GetParameter(fEtaXPar);
		double fEtaX_err = f_full->GetParError(fEtaXPar);
		
		double rEtaX     = f_full->GetParameter(rEtaXPar);
		double rEtaX_err = f_full->GetParError(rEtaXPar);
		
		fraction    = fEtaX * (1.0-rEtaX);
		fractionErr = sqrt(pow(fEtaX_err*(1.0-rEtaX),2.0) + pow(fEtaX*rEtaX_err,2.0));
		return;
	}
	
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
	
	if(fitOption_signal==3)
	{
		int fEtaXPar = f_full->GetParNumber("f_{#etaX}");
		int rEtaXPar = f_full->GetParNumber("r_{#etaX}");
		
		double fEtaX     = f_full->GetParameter(fEtaXPar);
		double fEtaX_err = f_full->GetParError(fEtaXPar);
		
		double rEtaX     = f_full->GetParameter(rEtaXPar);
		double rEtaX_err = f_full->GetParError(rEtaXPar);
		
		fraction    = fEtaX * rEtaX;
		fractionErr = sqrt(pow(fEtaX_err*rEtaX,2.0) + pow(fEtaX*rEtaX_err,2.0));
		return;
	}
	
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
		case 3:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("z_{qf}");
			parameters.push_back("f_{#etaX}");
			parameters.push_back("r_{#etaX}");
			nParameters += 5;
			break;
		}
	}
	return nParameters;
}

int MggFitter::GetOmegaParameters(vector<TString> &parameters)
{
	// Rho0 and Omega backgrounds are treated together
	
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
	
	switch(fitOption_rho) {
		case 0:
			break;
		case 1:
		{
			parameters.push_back("N_{#rho}");
			nParameters++;
			if((fitOption_omega!=2) && (fitOption_omega!=3)) {
				parameters.push_back("#Delta#mu_{#omega}");
				nParameters++;
			}
			break;
		}
	}
	
	return nParameters;
}

int MggFitter::GetEMParameters(vector<TString> &parameters)
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

int MggFitter::GetBeamlineParameters(vector<TString> &parameters)
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
