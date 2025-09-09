#include "MggFitter.h"
#include "EtaAnalyzer.h"
#include "CrossSection.h"

//================================================================================================================//
// Constructor:

MggFitter::MggFitter() :
	h_full(nullptr), h_empty(nullptr), h_emptyWide(nullptr), 
	f_full(nullptr), f_empty(nullptr), f_emptyWide(nullptr),
	h_hadronicBkgdLineshape(nullptr), f_hadronicBkgdLineshape(nullptr), 
	h_etaLineshape(nullptr),          f_etaLineshape(nullptr),
	h_omegaLineshape(nullptr),        f_omegaLineshape(nullptr),
	h_rhoLineshape(nullptr), 
	h_etaPionLineshape(nullptr),      f_etaPionLineshape(nullptr), 
	h_fdcOmegaLineshape(nullptr),     f_fdcOmegaLineshape(nullptr)
{
	f_chebyshev = new TF1("chebyshev", "cheb5", 0.0, 1.0);
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
			- nominally: 1 normalization and 1 shift parameter for elastic signal
			- nominally: 2 normalization (eta+pi and eta+2pi) and 1 shift parameter for hadronic bkgd
		
		- Omega(->pi0+gamma)::
			- nominally: 5 for crystal ball function
		
		- Smooth (non-peaking) Background:
			- nominally: 4 parameters for exponential background
		
		- Empty target Background:
			- nominally: 5x3 for FDC peaking structures (Gaussians), 4 for smooth background
	
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
	m_parametersFull.push_back("angle");
	m_parametersFull.insert(m_parametersFull.end(), locSignalParameters.begin(),   locSignalParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locOmegaParameters.begin(),    locOmegaParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locBkgdParameters.begin(),     locBkgdParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locEtaPrimeParameters.begin(), locEtaPrimeParameters.end());
	m_parametersFull.insert(m_parametersFull.end(), locEmptyParameters.begin(),    locEmptyParameters.end());
	m_parametersFull.push_back("binSize");
	m_parametersFull.push_back("binSize_empty");
	
	int nTotalParameters = (int)m_parametersFull.size();
	
	// Now we need to initial vectors which store the index of each parameter that is used across both functions:
	
	m_parIndexFull.clear();
	m_parIndexFull.push_back(0);
	
	m_parIndexEmpty.clear();
	m_parIndexEmpty.push_back(0);
	
	for(int ipar=1; ipar<nTotalParameters; ipar++) {
		
		// all but the last parameter are used when fitting the full target mgg spectrum:
		if(ipar<(nTotalParameters-1)) m_parIndexFull.push_back(ipar);
		
		// only the parameters within the 'locEmptyParameters' vector are used in the empty target fit
		if(find(locEmptyParameters.begin(), locEmptyParameters.end(), m_parametersFull[ipar]) != locEmptyParameters.end()) {
			m_parIndexEmpty.push_back(ipar);
		}
		// as well as the last parameter:
		if(ipar==(nTotalParameters-1)) m_parIndexEmpty.push_back(ipar);
	}
}

//================================================================================================================//
// Main fit routine:

void MggFitter::FitData()
{
	int debug = 0;
	
	/*-----------------------------*/
	// First thing to do is set up vector's storing parameters for full target and empty target fits:
	
	if(debug) printf("\nInitializing parameter arrays...\n");
	InitializeParameterArrays();
	
	/*-----------------------------*/
	// Next, initialize the fit functions themselves:
	
	if(debug) printf("Initializing fit functions...\n");
	InitializeFitFunction(&f_full);
	InitializeEmptyFitFunction(&f_empty);
	
	//=======================================================================================================//
	// Define the struct needed to compute global chi2 from combined fit:
	
	struct GlobalChi2 {
		const ROOT::Math::IMultiGenFunction *chiSqFull;
		const ROOT::Math::IMultiGenFunction *chiSqEmpty;
		vector<int> indexFull, indexEmpty;
		
		GlobalChi2(const ROOT::Math::IMultiGenFunction &c1,
			const ROOT::Math::IMultiGenFunction &c2,
			const vector<int> &index1, const vector<int> &index2) 
			: chiSqFull(&c1), chiSqEmpty(&c2), indexFull(index1), indexEmpty(index2) {}
		GlobalChi2() {};
		
		double operator()(const double* p) const {
			vector<double> p1(indexFull.size()), p2(indexEmpty.size());
			for(size_t ipar=0; ipar<indexFull.size(); ipar++)
				p1[ipar] = p[indexFull[ipar]];
			for(size_t jpar=0; jpar<indexEmpty.size(); jpar++)
				p2[jpar] = p[indexEmpty[jpar]];
			return (*chiSqFull)(p1.data()) + (*chiSqEmpty)(p2.data());
		}
	};
	
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
	
	ROOT::Fit::DataOptions locOpts;
	
	//=======================================================================================================//
	// For the first fit we only let the smooth, in-target background float:
	
	ROOT::Fit::Fitter locFitter1;
	locFitter1.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter1.Config().SetMinimizer("Minuit", "Migrad");
	
	InitializeFitParameters(locFitter1);
	//UpdateFitFunctions(locFitter1);
	//return;
	
	ROOT::Fit::DataRange fittingRange1(minFitRange, maxFitRange);
	
	excludeRegions.push_back({0.51,0.60});
	excludeRegions.push_back({0.70,0.85});
	
	ROOT::Fit::BinData dataFull1(locOpts, fittingRange1);
	ROOT::Fit::FillData(dataFull1,  h_full);
	
	ROOT::Fit::BinData dataEmpty1(locOpts, fittingRange1);
	ROOT::Fit::FillData(dataEmpty1, h_emptyWide);
	
	ROOT::Math::WrappedMultiTF1 wfFull1(*f_full, f_full->GetNdim());
	ROOT::Math::WrappedMultiTF1 wfEmpty1(*f_empty, f_empty->GetNdim());
	
	ROOT::Fit::Chi2Function locChi2Full1(dataFull1,  wfFull1);
	ROOT::Fit::Chi2Function locChi2Empty1(dataEmpty1, wfEmpty1);
	
	GlobalChi2 locChi2FCN1(locChi2Full1, locChi2Empty1, m_parIndexFull, m_parIndexEmpty);
	
	ReleaseBkgdParameters(locFitter1);
	
	locFitter1.SetFCN(locNPars, locChi2FCN1, nullptr, dataFull1.Size() + dataEmpty1.Size(), 1);
	locFitter1.FitFCN();
	
	auto result1 = locFitter1.Result();
	if(debug) {
		printf("\n\nFit Results (1st fit):\n");
		result1.Print(std::cout);
	}
	//UpdateFitFunctions(result1);
	//return;
	
	//=======================================================================================================//
	// Now, let the empty target background parameters float if they are within this fitting range:
	
	ROOT::Fit::Fitter locFitter2;
	locFitter2.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter2.Config().SetMinimizer("Minuit", "Migrad");
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter2, locFitter1);
	
	ROOT::Fit::DataRange fittingRange2(minFitRange, maxFitRange);
	
	ROOT::Fit::BinData dataFull2(locOpts, fittingRange2);
	ROOT::Fit::FillData(dataFull2,  h_full);
	
	ROOT::Fit::BinData dataEmpty2(locOpts, fittingRange2);
	ROOT::Fit::FillData(dataEmpty2, h_emptyWide);
	
	ROOT::Math::WrappedMultiTF1 wfFull2(*f_full, f_full->GetNdim());
	ROOT::Math::WrappedMultiTF1 wfEmpty2(*f_empty, f_empty->GetNdim());
	
	ROOT::Fit::Chi2Function locChi2Full2(dataFull2,  wfFull2);
	ROOT::Fit::Chi2Function locChi2Empty2(dataEmpty2, wfEmpty2);
	
	GlobalChi2 locChi2FCN2(locChi2Full2, locChi2Empty2, m_parIndexFull, m_parIndexEmpty);
	
	// Only do the combined fit in the region where empty background is large:
	if(angle<2.0) ReleaseEmptyParameters(locFitter2);
	
	locFitter2.SetFCN(locNPars, locChi2FCN2, nullptr, dataFull2.Size() + dataEmpty2.Size(), 1);
	locFitter2.FitFCN();
	
	auto result2 = locFitter2.Result();
	if(debug) {
		printf("\n\nFit Results (2nd fit):\n");
		result2.Print(std::cout);
	}
	//UpdateFitFunctions(result2);
	//return;
	
	//=======================================================================================================//
	// Next, fit the region to the right of the peak, allowing the omega parameters to float:
	
	excludeRegions.clear();
	excludeRegions.push_back({0.51,0.60});
	
	ROOT::Fit::Fitter locFitter3;
	locFitter3.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter3.Config().SetMinimizer("Minuit", "Migrad");
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter3, locFitter2);
	
	// Fix parameters that were floating in previous step:
	
	ROOT::Fit::DataRange fittingRange3(minFitRange, maxFitRange);
	
	ROOT::Fit::BinData dataFull3(locOpts, fittingRange3);
	ROOT::Fit::FillData(dataFull3,  h_full);
	
	ROOT::Fit::BinData dataEmpty3(locOpts, fittingRange3);
	ROOT::Fit::FillData(dataEmpty3, h_emptyWide);
	
	ROOT::Math::WrappedMultiTF1 wfFull3(*f_full, f_full->GetNdim());
	ROOT::Math::WrappedMultiTF1 wfEmpty3(*f_empty, f_empty->GetNdim());
	
	ROOT::Fit::Chi2Function locChi2Full3(dataFull3,  wfFull3);
	ROOT::Fit::Chi2Function locChi2Empty3(dataEmpty3, wfEmpty3);
	
	GlobalChi2 locChi2FCN3(locChi2Full3, locChi2Empty3, m_parIndexFull, m_parIndexEmpty);
	
	FixBkgdParameters(locFitter3);
	FixEmptyParameters(locFitter3);
	ReleaseOmegaParameters(locFitter3);
	
	locFitter3.SetFCN(locNPars, locChi2FCN3, nullptr, dataFull3.Size() + dataEmpty3.Size(), 1);
	locFitter3.FitFCN();
	
	auto result3 = locFitter3.Result();
	if(debug) {
		printf("\n\nFit Results (3rd fit):\n");
		result3.Print(std::cout);
	}
	//UpdateFitFunctions(result3);
	//return;
	
	//=======================================================================================================//
	// Next, fit the eta mass region with everything else fixed:
	
	excludeRegions.clear();
	
	ROOT::Fit::Fitter locFitter4;
	locFitter4.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter4.Config().SetMinimizer("Minuit", "Migrad");
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter4, locFitter3);
	
	// Reconfigure fitter:
	
	ROOT::Fit::DataRange fittingRange4(minFitRange, maxFitRange);
	
	ROOT::Fit::BinData dataFull4(locOpts, fittingRange4);
	ROOT::Fit::FillData(dataFull4,  h_full);
	
	ROOT::Fit::BinData dataEmpty4(locOpts, fittingRange4);
	ROOT::Fit::FillData(dataEmpty4, h_emptyWide);
	
	ROOT::Math::WrappedMultiTF1 wfFull4(*f_full, f_full->GetNdim());
	ROOT::Math::WrappedMultiTF1 wfEmpty4(*f_empty, f_empty->GetNdim());
	
	ROOT::Fit::Chi2Function locChi2Full4(dataFull4,  wfFull4);
	ROOT::Fit::Chi2Function locChi2Empty4(dataEmpty4, wfEmpty4);
	
	GlobalChi2 locChi2FCN4(locChi2Full4, locChi2Empty4, m_parIndexFull, m_parIndexEmpty);
	
	FixOmegaParameters(locFitter4);
	ReleaseEtaParameters(locFitter4);
	
	locFitter4.SetFCN(locNPars, locChi2FCN4, nullptr, dataFull4.Size() + dataEmpty4.Size(), 1);
	locFitter4.FitFCN();
	
	auto result4 = locFitter4.Result();
	if(debug) {
		printf("\n\nFit Results (4th fit):\n");
		result4.Print(std::cout);
	}
	//UpdateFitFunctions(result4);
	//return;
	
	//=======================================================================================================//
	// Finally, let everything float:
	
	ROOT::Fit::Fitter locFitter5;
	locFitter5.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter5.Config().SetMinimizer("Minuit", "Migrad");
	
	// setup parameter settings based on previous fit:
	SetFitParameters(locFitter5, locFitter4);
	
	// Reconfigure fitter:
	
	ROOT::Fit::DataRange fittingRange5(minFitRange, maxFitRange);
	
	ROOT::Fit::BinData dataFull5(locOpts, fittingRange5);
	ROOT::Fit::FillData(dataFull5,  h_full);
	
	ROOT::Fit::BinData dataEmpty5(locOpts, fittingRange5);
	ROOT::Fit::FillData(dataEmpty5, h_emptyWide);
	
	ROOT::Math::WrappedMultiTF1 wfFull5(*f_full, f_full->GetNdim());
	ROOT::Math::WrappedMultiTF1 wfEmpty5(*f_empty, f_empty->GetNdim());
	
	ROOT::Fit::Chi2Function locChi2Full5(dataFull5,  wfFull5);
	ROOT::Fit::Chi2Function locChi2Empty5(dataEmpty5, wfEmpty5);
	
	GlobalChi2 locChi2FCN5(locChi2Full5, locChi2Empty5, m_parIndexFull, m_parIndexEmpty);
	
	ReleaseOmegaParameters(locFitter5);
	ReleaseBkgdParameters(locFitter5);
	
	// Only do the combined fit in the region where empty background is large:
	if(angle<2.0) ReleaseEmptyParameters(locFitter5);
	
	locFitter5.SetFCN(locNPars, locChi2FCN5, nullptr, dataFull5.Size() + dataEmpty5.Size(), 3);
	//locFitter5.MinFCN();
	
	if(debug) printf("Performing fit...\n");
	bool ok = locFitter5.FitFCN();
	
	auto result = locFitter5.Result();
	result.Print(std::cout);
	if (!ok) std::cerr << "Fit did not converge!" << std::endl;
	
	UpdateFitFunctions(result);
	
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
			f_full->SetParError(locParNumber, locParErr);
			
			// check if f_full->GetParameter() agrees:
			double locDiff = locParVal - f_full->GetParameter(locParNumber);
			if(fabs(locDiff-0.0)>1.e-6) {
				printf("\n\n");
				printf("=================================================================\n");
				printf("WARNING!!! f_full parameter values do not align with fit results!\n");
				printf("   par: %s, f_full->GetParameter() = %f, result.Parameter() = %f\n", 
					m_parametersFull[ipar].Data(), f_full->GetParameter(locParNumber), locParVal);
				printf("=================================================================\n");
				printf("\n\n");
			}
		}
		
		// Do the same for f_empty (even though I don't think it's needed anywhere):
		if(find(m_parIndexEmpty.begin(), m_parIndexEmpty.end(), ipar) != m_parIndexEmpty.end()) {
			int locParNumber = f_empty->GetParNumber(m_parametersFull[ipar].Data());
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
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int pPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), Form("p%d",ipar)) - m_parametersFull.begin());
				fitter.Config().ParSettings(pPar).Release();
				fitter.Config().ParSettings(pPar).SetLimits(-1.e6, 1.e6);
			}
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
				fitter.Config().ParSettings(locMuPar).SetLimits(m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				
				fitter.Config().ParSettings(locSigPar).Release();
				fitter.Config().ParSettings(locSigPar).SetLimits(0.0075, 0.025);
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
		case 12:
		{
			int NPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "N_{#eta}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NPar).Release();
			fitter.Config().ParSettings(NPar).SetLimits(0.0, 1.e6);
			
			int NEtaPiPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi}") - m_parametersFull.begin());
			fitter.Config().ParSettings(NEtaPiPar).Release();
			fitter.Config().ParSettings(NEtaPiPar).SetLimits(0.0, 10.0);
			
			int NHadBgPar = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), "A_{#eta#pi#pi}") - m_parametersFull.begin());
			if(m_hadronicBkgdYieldBGGEN>(0.02*fitter.Config().ParSettings(NPar).Value())) {
				fitter.Config().ParSettings(NHadBgPar).Release();
				fitter.Config().ParSettings(NHadBgPar).SetLimits(0.0, 10.0);
			}
			if(vetoOption==6) {
				fitter.Config().ParSettings(NHadBgPar).SetValue(1.0);
				fitter.Config().ParSettings(NHadBgPar).Fix();
			}
			break;
		}
	}
}


//==============================================================//
// Lineshape Fits:

void MggFitter::FitEtaLineshape(int drawOption)
{
	if(h_etaLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_etaLineshape, "Eta Lineshape");
	/*
	for(int ibin=1; ibin<=h_etaLineshape->GetXaxis()->GetNbins(); ibin++) {
		double locX = h_etaLineshape->GetBinCenter(ibin);
		if((locX<minMggCut) || (locX>maxMggCut)) {
			h_etaLineshape->SetBinContent(ibin, 0.0);
			h_etaLineshape->SetBinError(ibin, 0.0);
		}
	}
	h_etaLineshape->Scale(1.0/h_etaLineshape->Integral());
	*/
	
	h_etaLineshape->SetMarkerStyle(8);
	h_etaLineshape->SetMarkerSize(0.7);
	h_etaLineshape->SetMarkerColor(kBlue);
	h_etaLineshape->SetLineColor(kBlue);
	h_etaLineshape->GetYaxis()->SetTitle(Form("Normalized Counts / %d MeV/c^{2}", (int)(1.e3 * h_etaLineshape->GetBinWidth(1))));
	h_etaLineshape->GetYaxis()->SetTitleSize(0.05);
	h_etaLineshape->GetYaxis()->SetTitleOffset(1.2);
	h_etaLineshape->GetYaxis()->CenterTitle(true);
	h_etaLineshape->GetXaxis()->SetTitleSize(0.05);
	h_etaLineshape->GetXaxis()->SetTitleOffset(1.0);
	h_etaLineshape->GetXaxis()->CenterTitle(true);
	h_etaLineshape->SetTitle("");
	
	switch(useRawMass) {
		case 0:
		{
			// Two Crystal Ball Functions + One Gaussian:
			
			TF1 *locfEta = new TF1("locfEta", DoubleCrystalBallPlusGausPDF, 0.30, 0.80, 13);
			locfEta->SetParameter(0, 0.540);  // mu1
			locfEta->SetParameter(1, 0.0075); // sigma1
			locfEta->SetParameter(2, 1.400);  // alpha1
			locfEta->SetParameter(3, 9.000);  // n1
			locfEta->SetParameter(4, 0.015);  // mu2-mu1
			locfEta->SetParameter(5, 0.010);  // sigma2
			locfEta->SetParameter(6, 1.400);  // alpha2
			locfEta->SetParameter(7, 13.50);  // n2
			locfEta->SetParameter(8, 0.535);  // mu3
			locfEta->SetParameter(9, 0.007);  // sigma3
			locfEta->SetParameter(10, 0.63);   // fraction1
			locfEta->SetParameter(11, 0.20);   // fraction2
			locfEta->FixParameter(12, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParName(0, "#mu_{1}");
			locfEta->SetParName(1, "#sigma_{1}");
			locfEta->SetParName(2, "#alpha_{1}");
			locfEta->SetParName(3, "n_{1}");
			locfEta->SetParName(4, "#mu_{2}-#mu_{1}");
			locfEta->SetParName(5, "#sigma_{2}");
			locfEta->SetParName(6, "#alpha_{2}");
			locfEta->SetParName(7, "n_{2}");
			locfEta->SetParName(8, "#mu_{3}-#mu_{1}");
			locfEta->SetParName(9, "#sigma_{3}");
			locfEta->SetParName(10, "frac1");
			locfEta->SetParName(11, "frac2");
			
			locfEta->SetParLimits(0, 0.530,  0.560);
			locfEta->SetParLimits(1, 0.004,  0.050);
			locfEta->SetParLimits(2, 0.200,  9.999);
			locfEta->SetParLimits(3, 1.100, 49.999);
			
			locfEta->SetParLimits(4,-0.050,  0.050);
			locfEta->SetParLimits(5, 0.004,  0.050);
			locfEta->SetParLimits(6, 0.200,  9.999);
			locfEta->SetParLimits(7, 1.100, 49.999);
			
			locfEta->SetParLimits(10, 0.0, 1.0);
			locfEta->SetParLimits(11, 0.0, 1.0);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(4, 0.0);
			//locfEta->FixParameter(8, 0.0);
			
			locfEta->FixParameter(8, 0.0);
			locfEta->FixParameter(9, 0.007);
			locfEta->FixParameter(11, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0Q");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				cEtaLS->SetLeftMargin(0.13); cEtaLS->SetRightMargin(0.07);
				cEtaLS->SetBottomMargin(0.13); cEtaLS->SetTopMargin(0.07);
				//gStyle->SetOptFit(0);
				//h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.5,0.6);
				h_etaLineshape->Draw();
				locfEta->SetRange(0.4,0.7);
				locfEta->SetNpx(1000);
				locfEta->Draw("same");
				//h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
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
				locLat.DrawLatexNDC(0.165, 0.725, Form("#int_{%.2f}^{%.2f}#color[632]{f_{#eta}}dm_{#gamma#gamma} = %.4f", 
					minMggCut, maxMggCut, fracAccepted));
				
				cEtaLS->Update();
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
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(12, 1.0);
			
			// Widen the lineshape by 10% (test):
			//f_etaLineshape->SetParameter(1, 1.1*f_etaLineshape->GetParameter(1));
			//f_etaLineshape->SetParameter(5, 1.1*f_etaLineshape->GetParameter(5));
			
			locfEta->Delete();
			break;
		}
		case 2:
		{
			// Two Gaussian Functions:
			
			TF1 *locfEta = new TF1("locfEta", DoubleGausPDF, 0.52, 0.59, 6);
			locfEta->SetParameters(
				0.547, //      mu1
				0.010, //   sigma1
				0.000, //  mu2-mu1
				0.015, //   sigma2
				0.300  // fraction
			);
			locfEta->FixParameter(5, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParLimits(0, 0.530,  0.570);
			locfEta->SetParLimits(1, 0.005,  0.050);
			locfEta->SetParLimits(2,-0.050,  0.050);
			locfEta->SetParLimits(3, 0.010,  0.050);
			locfEta->SetParLimits(4, 0.000,  1.000);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(2, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0QL");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				locfEta->SetRange(0.4,0.7);
				h_etaLineshape->Draw();
				locfEta->Draw("same");
				h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
				printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
				
				cEtaLS->Update();
				getchar();
				cEtaLS->SetLogy();
				cEtaLS->Update();
				getchar();
				delete cEtaLS;
			}
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleGausPDF, minFitRange, maxFitRange, 6);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(5, 1.0);
			
			locfEta->Delete();
			break;
		}
		case 1:
		{
			// Two Gaussian Functions:
			
			TF1 *locfEta = new TF1("locfEta", DoubleGausPDF, 0.40, 0.70, 6);
			locfEta->SetParameters(
				0.547, //      mu1
				0.015, //   sigma1
				0.000, //  mu2-mu1
				0.025, //   sigma2
				0.300  // fraction
			);
			locfEta->FixParameter(5, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParLimits(0, 0.530,  0.570);
			locfEta->SetParLimits(1, 0.005,  0.050);
			locfEta->SetParLimits(2,-0.050,  0.050);
			locfEta->SetParLimits(3, 0.010,  0.050);
			locfEta->SetParLimits(4, 0.000,  1.000);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(2, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0QL");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				locfEta->SetRange(0.4,0.7);
				h_etaLineshape->Draw();
				locfEta->Draw("same");
				h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
				printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
				
				cEtaLS->Update();
				getchar();
				cEtaLS->SetLogy();
				cEtaLS->Update();
				getchar();
				delete cEtaLS;
			}
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleGausPDF, minFitRange, maxFitRange, 6);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(5, 1.0);
			
			locfEta->Delete();
			break;
		}
	}
	return;
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
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_omegaLineshape, "Omega Lineshape");
	
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

void MggFitter::FitFDCOmegaLineshape(int drawOption)
{
	if(h_fdcOmegaLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_fdcOmegaLineshape, "FDC Omega Lineshape");
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBallPDF, 0.2, 0.7, 5);
	
	double    muGuess = h_fdcOmegaLineshape->GetBinCenter(h_fdcOmegaLineshape->GetMaximumBin());
	double sigmaGuess = 0.025;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fOmega1->SetParameter(0,    muGuess);
	fOmega1->SetParameter(1, sigmaGuess);
	fOmega1->SetParameter(2, alphaGuess);
	fOmega1->SetParameter(3,     nGuess);
	
	fOmega1->SetParLimits(0, muGuess-0.05, muGuess+0.05);
	fOmega1->SetParLimits(1, 0.005,  0.050);
	fOmega1->SetParLimits(2, 0.200,  9.999);
	fOmega1->SetParLimits(3, 1.100, 49.999);
	
	fOmega1->FixParameter(4, h_fdcOmegaLineshape->GetBinWidth(1));
	
	h_fdcOmegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", DoubleCrystalBallPDF, 0.2, 0.7, 10);
	fOmega2->SetParameter(0, fOmega1->GetParameter(0));
	fOmega2->SetParameter(1, fOmega1->GetParameter(1));
	fOmega2->SetParameter(2, fOmega1->GetParameter(2));
	fOmega2->SetParameter(3, fOmega1->GetParameter(3));
	fOmega2->SetParameter(4, 0.0);
	fOmega2->SetParameter(5, fOmega1->GetParameter(1)*2.0);
	fOmega2->SetParameter(6, fOmega1->GetParameter(2));
	fOmega2->SetParameter(7, fOmega1->GetParameter(3));
	fOmega2->SetParameter(8, 0.0);
	
	fOmega2->SetParLimits(0,  muGuess-0.05, muGuess+0.05);
	fOmega2->SetParLimits(1,  0.005,  0.050);
	fOmega2->SetParLimits(2,  0.200,  9.999);
	fOmega2->SetParLimits(3,  1.100, 49.999);
	fOmega2->SetParLimits(4, -0.050,  0.050);
	fOmega2->SetParLimits(5,  0.005,  0.150);
	fOmega2->SetParLimits(6,  0.200,  9.999);
	fOmega2->SetParLimits(7,  1.100, 49.999);
	fOmega2->SetParLimits(8,  0.000,  1.000);
	
	fOmega2->FixParameter(9, h_fdcOmegaLineshape->GetBinWidth(1));
	
	h_fdcOmegaLineshape->Fit(fOmega2, "R0QL");
	
	if(drawOption) {
		TCanvas *cFDCOmegaLS = new TCanvas("cFDCOmegaLS", "cFDCOmegaLS", 950, 700);
		//cFDCOmegaLS->SetLogy();
		h_fdcOmegaLineshape->Draw();
		fOmega2->Draw("same");
		
		cFDCOmegaLS->Update();
		getchar();
		delete cFDCOmegaLS;
	}
	
	f_fdcOmegaLineshape = new TF1("f_fdcOmegaLineshape", DoubleCrystalBallPDF, minFitRange, maxFitRange, 10);
	f_fdcOmegaLineshape->SetParameters(fOmega2->GetParameters());
	f_fdcOmegaLineshape->FixParameter(9, 1.0);
	
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
	
	int minMggBin = h_full->FindBin(minMggCut+0.001*binSize);
	int maxMggBin = h_full->FindBin(maxMggCut-0.001*binSize);
	
	double locMinMggCut = h_full->GetBinCenter(minMggBin) - 0.5*binSize;
	double locMaxMggCut = h_full->GetBinCenter(maxMggBin) + 0.5*binSize;
	
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
			double locData = h_full->GetBinContent(ibin);
			double locBkgd = locfBkgd->Eval(h_full->GetBinCenter(ibin));
			double locFit  = f_full->Eval(h_full->GetBinCenter(ibin));
			fitCounts  += locFit;
			dataCounts += locData;
			yield    += locData - locBkgd;
			yieldErr += pow(h_full->GetBinError(ibin),2.0) + locBkgd;
			//printf("  mgg, data, fit: %f %f %f\n", h_full->GetBinCenter(ibin), locData, locFit);
		}
		yieldErr = sqrt(yieldErr);
		
		if(verbose) printf("\nDifference between data and fit: %f sigma\n", (dataCounts-fitCounts)/yieldErr);
	}
	else {
		switch(fitOption_signal) {
			case 1:
			{
				// Single Gaussian:
				int yieldPar    = f_full->GetParNumber("N_{#eta}");
				
				double locMu    = f_full->GetParameter("#mu_{#eta}");
				double locSigma = f_full->GetParameter("#sigma_{#eta}");
				double locInt   = IntegrateGaussian(locMu, locSigma, locMinMggCut, locMaxMggCut);
				
				double locN     = f_full->GetParameter(yieldPar);
				double locNErr  = f_full->GetParError(yieldPar);
				
				yield    = locN    * locInt;
				yieldErr = locNErr * locInt;
				break;
			}
			case 2:
			{
				// Double Gaussian:
				int yieldPar     = f_full->GetParNumber("N_{#eta}");
				int fractionPar  = f_full->GetParNumber("fraction_{#eta}");
				
				double locYield       = f_full->GetParameter(yieldPar);
				double locYieldErr    = f_full->GetParError(yieldPar);
				double locFraction    = f_full->GetParameter(fractionPar);
				double locFractionErr = f_full->GetParError(fractionPar);
				
				double locMu1    = f_full->GetParameter("#mu_{#eta,1}");
				double locSigma1 = f_full->GetParameter("#sigma_{#eta,1}");
				double locInt1   = IntegrateGaussian(locMu1, locSigma1, locMinMggCut, locMaxMggCut);
				
				double locMu2    = f_full->GetParameter("#mu_{#eta,2}-#mu_{#eta,1}") + locMu1;
				double locSigma2 = f_full->GetParameter("#sigma_{#eta,2}");
				double locInt2   = IntegrateGaussian(locMu2, locSigma2, locMinMggCut, locMaxMggCut);
				
				double locInt    = (1.0 - locFraction)*locInt1 + locFraction*locInt2;
				double locIntErr = locFractionErr * fabs(locInt2 - locInt1);
				
				printf("Int1, Int2, locInt = %f, %f, %f\n", locInt1, locInt2, locInt);
				
				yield    = locYield * locInt;
				yieldErr = sqrt(pow(locYieldErr*locInt,2.0) + pow(locYield*locIntErr,2.0));
				break;
			}
			case 3:
			{
				// Crystal Ball:
				TF1 *locPDF = new TF1("locCrystalBallPDF", CrystalBallPDF_flip, 0.5, 0.6, 5);
				locPDF->SetParameters(
					f_full->GetParameter("#mu_{#eta}"),
					f_full->GetParameter("#sigma_{#eta}"),
					f_full->GetParameter("#alpha_{#eta}"),
					f_full->GetParameter("n_{#eta}"),
					1.0
				);
				double locInt = locPDF->Integral(locMinMggCut, locMaxMggCut);
				
				int yieldPar = f_full->GetParNumber("N_{#eta}");
				
				yield    = locInt * f_full->GetParameter(yieldPar);
				yieldErr = locInt * f_full->GetParError(yieldPar);
				delete locPDF;
				break;
			}
			case 4:
			{
				// Crystal Ball + Gaussian (not yet fully implemented):
				yield    = 0.0;
				yieldErr = 0.0;
				break;
			}
			case 5:
			{
				// Lineshape fit (with histogram):
				int yieldPar = f_full->GetParNumber("N_{#eta}");
				
				int locMinMggBin = h_etaLineshape->FindBin(locMinMggCut+0.5*h_etaLineshape->GetBinWidth(1));
				int locMaxMggBin = h_etaLineshape->FindBin(locMaxMggCut-0.5*h_etaLineshape->GetBinWidth(1));
				
				double lsMinMggCut = h_etaLineshape->GetBinCenter(locMinMggBin) - 0.5*h_etaLineshape->GetBinWidth(1);
				double lsMaxMggCut = h_etaLineshape->GetBinCenter(locMaxMggBin) + 0.5*h_etaLineshape->GetBinWidth(1);
				
				// check that lsMinMggCut and lsMaxMggCut align with desired cut range:
				if((fabs(lsMinMggCut-locMinMggCut)>1.e-6) || (fabs(lsMaxMggCut-locMaxMggCut)>1.e-6)) {
					printf("\n\nWarning: Mgg binning of eta lineshape does not overlap with data.\n");
					printf("  Lineshape cut range: %f GeV - %f GeV\n", lsMinMggCut, lsMaxMggCut);
					printf("  Data cut range: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
				}
				
				yield    = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_full->GetParameter(yieldPar);
				yieldErr = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_full->GetParError(yieldPar);
				break;
			}
			case 6:
			{
				// Eta Lineshape (PDF) + Eta+Pion Lineshape (PDF):
				
				int yieldPar          = f_full->GetParNumber("N_{#eta}");
				double N_eta          = f_full->GetParameter(yieldPar);
				double N_etaErr       = f_full->GetParError(yieldPar);
				double locCorrection  = f_etaLineshape->Integral(minMggCut, maxMggCut);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				int bkgdPar           = f_full->GetParNumber("N_{#eta,bkgd}");
				double N_etaBkgd      = f_full->GetParameter(bkgdPar);
				double N_etaBkgdErr   = f_full->GetParError(bkgdPar);
				
				TF1 *flocClone;
				InitializeFitFunction(&flocClone,"locf1");
				flocClone->SetParameters(f_full->GetParameters());
				ZeroSignalPars(flocClone, 1);
				ZeroOmegaPars(flocClone);
				ZeroBkgdPars(flocClone);
				ZeroEtaPrimePars(flocClone);
				ZeroEmptyPars(flocClone);
				
				double locYieldBkgd    = flocClone->Integral(minMggCut,maxMggCut);
				double locYieldBkgdErr = locYieldBkgd * (N_etaBkgdErr/N_etaBkgd);
				
				double locYieldInc    = locYieldEta + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 7:
			{
				// Eta Lineshape (PDF) + Hadronic Lineshape (PDF):
				
				int yieldPar    = f_full->GetParNumber("N_{#eta}");
				int fractionPar = f_full->GetParNumber("frac_{bkgd}");
				
				double lsShift  = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double frac_bkgd    = f_full->GetParameter(fractionPar);
				double frac_bkgdErr = f_full->GetParError(fractionPar);
				
				double N_bkgd    = N_eta * frac_bkgd;
				double N_bkgdErr = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta' and 'N_etapi' above represent the yield of exclusive eta's and eta+pions integrated over all mgg.
				But to be consistent with our efficiency correction that will be applied later, we need to correct
				this for the small fraction of events that fall outside our mgg cut range.
				
				Recall that at this point f_etaLineshape and f_etaPionLineshape are normalized to have unit-integrals.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrection_bkgd = f_hadronicBkgdLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldHadronicBkgd    = locCorrection_bkgd * N_bkgd;
				double locYieldHadronicBkgdErr = locCorrection_bkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldHadronicBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldHadronicBkgdErr,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 8:
			{
				// Eta Lineshape (PDF) + Hadronic Lineshape (PDF):
				
				int yieldPar    = f_full->GetParNumber("N_{#eta}");
				int fractionPar = f_full->GetParNumber("frac_{bkgd}");
				
				double lsShift  = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double frac_bkgd    = f_full->GetParameter(fractionPar);
				double frac_bkgdErr = f_full->GetParError(fractionPar);
				
				double N_bkgd    = N_eta * frac_bkgd;
				double N_bkgdErr = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta' and 'N_etapi' above represent the yield of exclusive eta's and eta+pions integrated over all mgg.
				But to be consistent with our efficiency correction that will be applied later, we need to correct
				this for the small fraction of events that fall outside our mgg cut range.
				
				Recall that at this point f_etaLineshape and f_etaPionLineshape are normalized to have unit-integrals.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrection_bkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldHadronicBkgd    = locCorrection_bkgd * N_bkgd;
				double locYieldHadronicBkgdErr = locCorrection_bkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldHadronicBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldHadronicBkgdErr,2.0));
				
				if(subtractHadBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 9:
			{
				// Eta (Histogram) + Eta+Pion Background (Histogram) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar         = f_full->GetParNumber("N_{#eta}");
				int fractionParEtaPi = f_full->GetParNumber("frac_{#eta#pi}");
				int fractionParBkgd  = f_full->GetParNumber("frac_{bkgd}");
				
				double lsShift    = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double frac_etapi    = f_full->GetParameter(fractionParEtaPi);
				double frac_etapiErr = f_full->GetParError(fractionParEtaPi);
				
				double frac_bkgd     = f_full->GetParameter(fractionParBkgd);
				double frac_bkgdErr  = f_full->GetParError(fractionParBkgd);
				
				double N_etapi    = N_eta * frac_etapi;
				double N_etapiErr = sqrt(pow(N_etaErr*frac_etapi,2.0) + pow(N_eta*frac_etapiErr,2.0));
				
				double N_bkgd     = N_eta * frac_bkgd;
				double N_bkgdErr  = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = h_etaLineshape->Integral(h_etaLineshape->FindBin(minMggCut-lsShift), h_etaLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = f_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift), 
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion    = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
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
			case 10:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar         = f_full->GetParNumber("N_{#eta}");
				int fractionParEtaPi = f_full->GetParNumber("frac_{#eta#pi}");
				int fractionParBkgd  = f_full->GetParNumber("frac_{bkgd}");
				
				double lsShift    = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double frac_etapi    = f_full->GetParameter(fractionParEtaPi);
				double frac_etapiErr = f_full->GetParError(fractionParEtaPi);
				
				double frac_bkgd     = f_full->GetParameter(fractionParBkgd);
				double frac_bkgdErr  = f_full->GetParError(fractionParBkgd);
				
				double N_etapi    = N_eta * frac_etapi;
				double N_etapiErr = sqrt(pow(N_etaErr*frac_etapi,2.0) + pow(N_eta*frac_etapiErr,2.0));
				
				double N_bkgd     = N_eta * frac_bkgd;
				double N_bkgdErr  = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				//double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locCorrection  = h_etaLineshape->Integral(h_etaLineshape->FindBin(minMggCut-lsShift), 
					h_etaLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift), 
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
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
			case 11:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar      = f_full->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_full->GetParNumber("A_{#eta#pi}");
				int yieldParBkgd  = f_full->GetParNumber("A_{#eta#pi#pi}");
				
				double lsShift    = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double N_etapi    = f_full->GetParameter(yieldParEtaPi) * m_etaPionYieldBGGEN;
				double N_etapiErr = f_full->GetParError(yieldParEtaPi) * m_etaPionYieldBGGEN;
				
				double N_bkgd     = f_full->GetParameter(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				double N_bkgdErr  = f_full->GetParError(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = f_etaPionLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
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
			case 12:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar      = f_full->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_full->GetParNumber("A_{#eta#pi}");
				int yieldParBkgd  = f_full->GetParNumber("A_{#eta#pi#pi}");
				
				double lsShift    = f_full->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_full->GetParameter(yieldPar);
				double N_etaErr = f_full->GetParError(yieldPar);
				
				double N_etapi    = f_full->GetParameter(yieldParEtaPi) * m_etaPionYieldBGGEN;
				double N_etapiErr = f_full->GetParError(yieldParEtaPi) * m_etaPionYieldBGGEN;
				
				double N_bkgd     = f_full->GetParameter(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				double N_bkgdErr  = f_full->GetParError(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift),
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
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
		
		ZeroEmptyOmegaPars(locfEta);
		ZeroEmptyBkgdPars(locfEta);
		ZeroEmptyFDCPars(locfEta);
		
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
	
	//-----------------------------------------------//
	
	int minMggBin = h_full->FindBin(minMggCut);
	int maxMggBin = h_full->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full->GetBinCenter(ibin));
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
	
	//-----------------------------------------------//
	
	int minMggBin = h_full->FindBin(minMggCut);
	int maxMggBin = h_full->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfBkgd->Eval(h_full->GetBinCenter(ibin));
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
	
	if(fitOption_signal<6) return;
	
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
	
	if(fitOption_signal>=9) {
		if(fitOption_signal>10) locfEta->SetParameter("A_{#eta#pi}", 0.0);
		else locfEta->SetParameter("frac_{#eta#pi}", 0.0);
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_full->FindBin(minMggCut);
	int maxMggBin = h_full->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int bkgdYieldPar;
	if(fitOption_signal>10) bkgdYieldPar = f_full->GetParNumber("A_{#eta#pi#pi}");
	else bkgdYieldPar = f_full->GetParNumber("frac_{bkgd}");
	
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
	
	if(fitOption_signal<9) return;
	
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
	
	if(fitOption_signal>10) locfEta->SetParameter("A_{#eta#pi#pi}", 0.0);
	else locfEta->SetParameter("frac_{bkgd}", 0.0);
	
	//-----------------------------------------------//
	
	int minMggBin = h_full->FindBin(minMggCut);
	int maxMggBin = h_full->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_full->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int etaPionYieldPar;
	if(fitOption_signal>10) etaPionYieldPar = f_full->GetParNumber("A_{#eta#pi}");
	else etaPionYieldPar = f_full->GetParNumber("frac_{#eta#pi}");
	
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
	
	if(fitOption_signal<6) return;
	
	int ipar = f_full->GetParNumber("#Delta#mu_{#eta}");
	shift    = f_full->GetParameter(ipar);
	shiftErr = f_full->GetParError(ipar);
	return;
}

void MggFitter::ZeroSignalPars(TF1 *f1, int excludeHadronicBkgd)
{
	switch(fitOption_signal) {
		case 1:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 2:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 3:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 4:
			f1->SetParameter("N_{#eta,1}", 0.0);
			f1->SetParameter("N_{#eta,2}", 0.0);
			break;
		case 5:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 6:
			f1->SetParameter("N_{#eta}", 0.0);
			if(!excludeHadronicBkgd) f1->SetParameter("N_{#eta,bkgd}", 0.0);
			break;
		case 7:
			if(excludeHadronicBkgd) {
				f1->SetParameter("frac_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			break;
		case 8:
			if(excludeHadronicBkgd) {
				f1->SetParameter("frac_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			break;
		case 9:
			if(excludeHadronicBkgd) {
				f1->SetParameter("frac_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			break;
		case 10:
			if(excludeHadronicBkgd) {
				f1->SetParameter("frac_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			break;
		case 11:
			if(excludeHadronicBkgd) {
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("A_{#eta#pi}", 0.0);
				f1->SetParameter("A_{#eta#pi#pi}",    0.0);
			}
			break;
		case 12:
			if(excludeHadronicBkgd) {
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("A_{#eta#pi}", 0.0);
				f1->SetParameter("A_{#eta#pi#pi}",    0.0);
			}
			break;
	}
	return;
}

void MggFitter::ZeroHadronicBkgdPars(TF1 *f1)
{
	switch(fitOption_signal) {
		case 6:
			f1->SetParameter("N_{#eta,bkgd}",  0.0);
			break;
		case 7:
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 8:
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 9:
			f1->SetParameter("frac_{#eta#pi}", 0.0);
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 10:
			f1->SetParameter("frac_{#eta#pi}", 0.0);
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 11:
			f1->SetParameter("A_{#eta#pi}",    0.0);
			f1->SetParameter("A_{#eta#pi#pi}", 0.0);
			break;
		case 12:
			f1->SetParameter("A_{#eta#pi}",    0.0);
			f1->SetParameter("A_{#eta#pi#pi}", 0.0);
			break;
	}
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
	ZeroEmptyBkgdPars(f1);
	ZeroEmptyFDCPars(f1);
	ZeroEmptyEtaPars(f1);
	ZeroEmptyOmegaPars(f1);	
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

void MggFitter::ZeroEmptyEtaPars(TF1 *f1)
{
	if(fitOption_empty==0) return;
	if(emptyFitOption_eta==0) return;
	f1->SetParameter("N_{#eta,empty}",0.0);
	return;
}

void MggFitter::ZeroEmptyOmegaPars(TF1 *f1)
{
	if(fitOption_empty==0) return;
	if(emptyFitOption_omega==0) return;
	f1->SetParameter("N_{#omega,empty}",0.0);
	return;
}

void MggFitter::FillPull(TH1F *h_pull)
{
	if(h_pull==NULL) {
		return;
	}
	
	// Just check that the binning of supplied histogram and h_full are consistent:
	if(h_pull->GetXaxis()->GetNbins() != h_full->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of pull histogram.\n" << endl;
	}
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_full->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_full->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_full->GetBinContent(ibin) - f_full->Eval(loc_mgg))/loc_unc);
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
	
	// redundant check:
	if(fitOption_signal<7) return;
	
	int fractionPar;
	if(fitOption_signal<11) {
		fractionPar = f_full->GetParNumber("frac_{bkgd}");
		fraction    = f_full->GetParameter(fractionPar);
		fractionErr = f_full->GetParError(fractionPar);
		return;
	}
	
	// For fit options 11 and 12, we don't have the fraction as a direct parameter of the fit function.
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
	
	// redundant check:
	if(fitOption_signal<9) return;
	
	int fractionPar;
	if(fitOption_signal<11) {
		fractionPar = f_full->GetParNumber("frac_{#eta#pi}");
		fraction    = f_full->GetParameter(fractionPar);
		fractionErr = f_full->GetParError(fractionPar);
		return;
	}
	
	// For fit options 11 and 12, we don't have the fraction as a direct parameter of the fit function.
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
			parameters.push_back("#mu_{#eta}");
			parameters.push_back("#sigma_{#eta}");
			nParameters += 3;
			break;
		}
		case 2:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("fraction_{#eta}");
			parameters.push_back("#mu_{#eta,1}");
			parameters.push_back("#mu_{#eta,2}-#mu_{#eta,1}");
			parameters.push_back("#sigma_{#eta,1}");
			parameters.push_back("#sigma_{#eta,2}");
			nParameters += 6;
			break;
		}
		case 3:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#mu_{#eta}");
			parameters.push_back("#sigma_{#eta}");
			parameters.push_back("#alpha_{#eta}");
			parameters.push_back("n_{#eta}");
			nParameters += 5;
			break;
		}
		case 4:
		{
			parameters.push_back("N_{#eta,1}");
			parameters.push_back("#mu_{#eta,1}");
			parameters.push_back("#sigma_{#eta,1}");
			parameters.push_back("#alpha_{#eta,1}");
			parameters.push_back("n_{#eta,1}");
			parameters.push_back("N_{#eta,2}");
			parameters.push_back("#mu_{#eta,2}");
			parameters.push_back("#sigma_{#eta,2}");
			nParameters += 8;
			break;
		}
		case 5:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			nParameters += 2;
			break;
		}
		case 6:
		{
			parameters.push_back(        "N_{#eta}"     );
			parameters.push_back("#Delta#mu_{#eta}"     );
			parameters.push_back(        "N_{#eta,bkgd}");
			parameters.push_back(      "#mu_{#eta,bkgd}");
			parameters.push_back(   "#sigma_{#eta,bkgd}");
			parameters.push_back(   "#alpha_{#eta,bkgd}");
			parameters.push_back(        "n_{#eta,bkgd}");
			nParameters += 7;
			break;
		}
		case 7:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("frac_{#eta}");
			parameters.push_back("frac_{bkgd}");
			nParameters += 4;
			break;
		}
		case 8:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("frac_{#eta}");
			parameters.push_back("frac_{bkgd}");
			nParameters += 4;
			break;
		}
		case 9:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("frac_{#eta}");
			parameters.push_back("frac_{#eta#pi}");
			parameters.push_back("frac_{bkgd}");
			nParameters += 5;
			break;
		}
		case 10:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("frac_{#eta}");
			parameters.push_back("frac_{#eta#pi}");
			parameters.push_back("frac_{bkgd}");
			nParameters += 5;
			break;
		}
		case 11:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("A_{#eta#pi}");
			parameters.push_back("A_{#eta#pi#pi}");
			nParameters += 4;
			break;
		}
		case 12:
		{
			parameters.push_back("N_{#eta}");
			parameters.push_back("#Delta#mu_{#eta}");
			parameters.push_back("A_{#eta#pi}");
			parameters.push_back("A_{#eta#pi#pi}");
			nParameters += 4;
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
			nParameters = 0;
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
	
	// Eta mesons produced from residual gas and target walls:
	switch(emptyFitOption_eta) {
		case 0:
			break;
		case 1:
			nParameters = 3;
			parameters.push_back("N_{#eta,empty}");
			parameters.push_back("#mu_{#eta,empty}");
			parameters.push_back("#sigma_{#eta,empty}");
			break;
		case 2:
			nParameters += 2;
			parameters.push_back("N_{#eta,empty}");
			parameters.push_back("#Delta#mu_{#eta,empty}");
			break;
	}
	
	// Omega mesons produced from residual gas and target walls:
	switch(emptyFitOption_omega) {
		case 0:
			break;
		case 1:
			nParameters += 5;
			parameters.push_back("N_{#omega,empty}");
			parameters.push_back("#mu_{#omega,empty}");
			parameters.push_back("#sigma_{#omega,empty}");
			parameters.push_back("#alpha_{#omega,empty}");
			parameters.push_back("n_{#omega,empty}");
			break;
		case 2:
			nParameters += 2;
			parameters.push_back("N_{#omega,empty}");
			parameters.push_back("#Delta#mu_{#omega,empty}");
			break;
	}
	
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
