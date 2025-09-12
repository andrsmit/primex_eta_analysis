#include "CrossSection.h"
#include "MggFitter.h"
#include "EtaAnalyzer.h"

int MggFitter::InitializeFitFunction(TF1 **f1, TString funcName)
{
	vector<TString> parNames; parNames.clear();
	
	int nParameters = 0;
	
	for(int ipar=0; ipar<(int)m_parIndexFull.size(); ipar++) {
		TString locParName = m_parametersFull[m_parIndexFull[ipar]];
		parNames.push_back(locParName);
		nParameters++;
	}
	
	//------------------------------------------//
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &MggFitter::MggFitFunction, 0.3, 1.20, nParameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < nParameters; ipar++) (*f1)->SetParName(ipar, parNames[ipar]);
	return nParameters;
}

void MggFitter::InitializeFitParameters(ROOT::Fit::Fitter &fitter)
{
	int debug = 0;
	
	// Initialize parameter values to something that will not produce nans:
	
	int nParameters      = 0;
	int nTotalParameters = (int)m_parametersFull.size();
	// (by the end of this function nParameters should be equal to nTotalParameters)
	
	if(debug) printf("  Total number of parameters: %d\n", nTotalParameters);
	
	// Start by setting all parameters to zero:
	double *dummyPars = new double[nTotalParameters];
	for(int ipar=0; ipar<(int)m_parametersFull.size(); ipar++) dummyPars[ipar] = 0.0;
	
	fitter.Config().SetParamsSettings(nTotalParameters, dummyPars);
	
	if(debug) printf("  Total number of parameter settings: %d\n", fitter.Config().NPar());
	if(debug) {
		for(int ipar=0; ipar<nTotalParameters; ipar++) {
			printf("    par %d: %s\n", ipar, m_parametersFull[ipar].Data());
		}
	}
	
	// The next set of parameters are for the signal + hadronic background:
	if(debug) printf("  guessing eta parameters...\n");
	
	vector<double> locSignalPars;
	GuessEtaParameters(locSignalPars);
	
	for(int ipar=0; ipar<(int)locSignalPars.size(); ipar++) {
		fitter.Config().ParSettings(nParameters+ipar).Set(m_parametersFull[nParameters+ipar].Data(), locSignalPars[ipar]);
	}
	nParameters += (int)locSignalPars.size();
	
	// Next are the parameters associated with the omega->pi0+gamma peak:
	if(debug) printf("  guessing omega parameters...\n");
	
	vector<double> locOmegaPars;
	GuessOmegaParameters(locOmegaPars);
	
	for(int ipar=0; ipar<(int)locOmegaPars.size(); ipar++) {
		fitter.Config().ParSettings(nParameters+ipar).Set(m_parametersFull[nParameters+ipar].Data(), locOmegaPars[ipar]);
	}
	nParameters += (int)locOmegaPars.size();
	
	// Next are the parameters for the smooth (non-peaking) background:
	if(debug) printf("  guessing bkgd parameters...\n");
	
	vector<double> locBkgdPars;
	GuessBkgdParameters(locBkgdPars);
	
	for(int ipar=0; ipar<(int)locBkgdPars.size(); ipar++) {
		fitter.Config().ParSettings(nParameters+ipar).Set(m_parametersFull[nParameters+ipar].Data(), locBkgdPars[ipar]);
	}
	nParameters += (int)locBkgdPars.size();
	
	// Next are the parameters for the eta-prime peak (usually this is ignored):
	if(debug) printf("  guessing eta-prime parameters...\n");
	
	vector<double> locEtaPrimePars;
	GuessEtaPrimeParameters(locEtaPrimePars);
	
	for(int ipar=0; ipar<(int)locEtaPrimePars.size(); ipar++) {
		fitter.Config().ParSettings(nParameters+ipar).Set(m_parametersFull[nParameters+ipar].Data(), locEtaPrimePars[ipar]);
	}
	nParameters += (int)locEtaPrimePars.size();
	
	// Finally are the empty target parameters (which at the time this function is called, should come from the initial fit to h_emptyWide)::
	if(debug) printf("  copying over empty parameters...\n");
	
	vector<double> locEmptyPars;
	GuessEmptyParameters(locEmptyPars);
	
	for(int ipar=0; ipar<(int)locEmptyPars.size(); ipar++) {
		fitter.Config().ParSettings(nParameters+ipar).Set(m_parametersFull[nParameters+ipar].Data(), locEmptyPars[ipar]);
		if(debug) printf("    got parameter %d (%s)\n", nParameters+ipar, m_parametersFull[nParameters+ipar].Data());
	}
	nParameters += (int)locEmptyPars.size();
	
	// The next parameter is for "A_{empty}":
	fitter.Config().ParSettings(nParameters).Set("A_{empty}", 0.04);
	nParameters++;
	
	// The next parameter is for "#alpha_{empty}":
	fitter.Config().ParSettings(nParameters).Set("#alpha_{empty}", 1.0);
	nParameters++;
	
	// The next parameter is for "#alpha_{acc}":
	fitter.Config().ParSettings(nParameters).Set("#alpha_{acc}", 1.0);
	nParameters++;
	
	// The next parameter is for "#alpha_{acc,switch}" 
	// (used to turn off accidental contribution from full without turning it off for empty as well):
	fitter.Config().ParSettings(nParameters).Set("#alpha_{acc,switch}", 1.0);
	nParameters++;
	
	fitter.Config().ParSettings(nParameters+0).Set(m_parametersFull[nParameters+0].Data(), binSize);
	fitter.Config().ParSettings(nParameters+1).Set(m_parametersFull[nParameters+1].Data(), emptyBinSize);
	
	if(debug) printf("    done.\n");
	return;
}

void MggFitter::SetFitParameters(ROOT::Fit::Fitter &fitterNew, ROOT::Fit::Fitter &fitterOld)
{
	
	int nTotalParameters = (int)m_parametersFull.size();
	
	// Start by setting all parameters to zero:
	double *dummyPars = new double[nTotalParameters];
	for(int ipar=0; ipar<(int)m_parametersFull.size(); ipar++) dummyPars[ipar] = 0.0;
	
	fitterNew.Config().SetParamsSettings(nTotalParameters, dummyPars);
	
	for(int ipar=0; ipar<nTotalParameters; ipar++) {
		fitterNew.Config().ParSettings(ipar).Set(fitterOld.Config().ParSettings(ipar).Name(), 
			fitterOld.Config().ParSettings(ipar).Value());
		
		if(!fitterOld.Config().ParSettings(ipar).IsFixed()) {
			fitterNew.Config().ParSettings(ipar).Release();
			fitterNew.Config().ParSettings(ipar).SetLimits(fitterOld.Config().ParSettings(ipar).LowerLimit(), 
				fitterOld.Config().ParSettings(ipar).UpperLimit());
			fitterNew.Config().ParSettings(ipar).SetStepSize(fitterOld.Config().ParSettings(ipar).StepSize());
		}
	}
	return;
}

void MggFitter::GuessEtaParameters(vector<double> &parGuesses)
{
	parGuesses.clear();
	
	// guess the number of etas by integrating histogram and subtracting background:
	
	double minEtaFit = 0.50;
	double maxEtaFit = 0.60;
	
	double NGuess = 0.0;
	for(int ihist=0; ihist<2; ihist++) {
		double weight = ihist==0 ? 1.0 : -0.1;
		NGuess += (weight * h_full[ihist]->Integral(h_full[ihist]->FindBin(minEtaFit), h_full[ihist]->FindBin(maxEtaFit)));
	}
	double MuGuess    = EtaAnalyzer::m_massEta;
	double SigmaGuess = 0.013;
	
	// subtract empty-target background based on fit to h_emptyWide:
	
	double nEmptyBkgd = h_emptyWide->Integral(h_emptyWide->FindBin(minEtaFit), h_emptyWide->FindBin(maxEtaFit));
	NGuess -= nEmptyBkgd;
	if(NGuess < 0.0) NGuess = 0.0;
	
	// Use peak position between 0.54-0.56 as 'mu' guess:
	double locMax = 0.0;
	for(int ibin=h_full[0]->FindBin(0.54); ibin<=h_full[0]->FindBin(0.56); ibin++) {
		if(h_full[0]->GetBinContent(ibin) > locMax) {
			locMax  = h_full[0]->GetBinContent(ibin);
			MuGuess = h_full[0]->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_signal) {
		case 11:
		{
			// Lineshape  (Signal fit + EtaPi fit + other hadronic bkgd hist): N_exc, deltaMu, A_etapi, A_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(lineshapeOffset);
			parGuesses.push_back(1.0); // fraction of quasifree lineshape
			parGuesses.push_back(0.0); // eta+pi yield
			parGuesses.push_back(1.0); // dummy parameter (should always be 1)
			parGuesses.push_back(0.0); // eta+pi+pi yield
			break;
		}
		case 12:
		{
			// Lineshape  (Signal fit + EtaPi hist + other hadronic bkgd hist): N_exc, deltaMu, A_etapi, A_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(lineshapeOffset);
			parGuesses.push_back(1.0); // fraction of quasifree lineshape
			parGuesses.push_back(0.0); // eta+pi yield
			parGuesses.push_back(1.0); // dummy parameter (should always be 1)
			parGuesses.push_back(0.0); // eta+pi+pi yield
			break;
		}
	}
	return;
}

void MggFitter::GuessOmegaParameters(vector<double> &parGuesses)
{
	parGuesses.clear();
	
	// guess the number of etas by integrating histogram and subtracting background:
	
	double minOmegaFit = 0.68;
	double maxOmegaFit = 0.85;
	
	double NGuess = 0.0;
	for(int ihist=0; ihist<2; ihist++) {
		double weight = ihist==0 ? 1.0 : -0.1;
		NGuess += (weight * h_full[ihist]->Integral(h_full[ihist]->FindBin(minOmegaFit), h_full[ihist]->FindBin(maxOmegaFit)));
	}
	double MuGuess    = EtaAnalyzer::m_massOmega;
	double SigmaGuess = 0.032;
	
	// subtract background from h_emptyWide:
	
	double nEmptyBkgd = h_emptyWide->Integral(h_emptyWide->FindBin(minOmegaFit), h_emptyWide->FindBin(maxOmegaFit));
	NGuess -= nEmptyBkgd;
	if(NGuess < 0.0) NGuess = 0.0;
	
	// Use peak position between 0.75-0.85 as 'mu' guess:
	/*
	double locMax = 0.0;
	for(int ibin=h_full[0]->FindBin(0.75); ibin<=h_full[0]->FindBin(0.85); ibin++) {
		if(h_full[0]->GetBinContent(ibin) > locMax) {
			locMax  = h_full[0]->GetBinContent(ibin);
			MuGuess = h_full[0]->GetBinCenter(ibin);
		}
	}
	*/
	switch(fitOption_omega) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// Crystal Ball: N, mu, sigma, alpha, n
			parGuesses.push_back(NGuess);
			parGuesses.push_back(MuGuess);
			parGuesses.push_back(SigmaGuess);
			parGuesses.push_back( 1.0);
			parGuesses.push_back(10.0);
			break;
		}
		case 2:
		{
			// Lineshape fit: N, deltaMu
			parGuesses.push_back(NGuess);
			parGuesses.push_back(lineshapeOffset);
			break;
		}
		case 3:
		{
			// Lineshape hist: N, deltaMu
			parGuesses.push_back(NGuess);
			parGuesses.push_back(lineshapeOffset);
			break;
		}
		case 4:
		{
			// Lineshape hist (with rho as separate fit parameter): N_omega, N_rho, deltaMu
			parGuesses.push_back(0.88*NGuess);
			parGuesses.push_back(0.12*NGuess);
			parGuesses.push_back(lineshapeOffset);
			break;
		}
		case 5:
		{
			// Double Crystal Ball: N, mu1, sigma1, alpha1, n1, (mu2-mu1), sigma2, alpha2, n2, frac
			parGuesses.push_back(NGuess);
			parGuesses.push_back(f_omegaLineshape->GetParameter(0));
			parGuesses.push_back(f_omegaLineshape->GetParameter(1));
			parGuesses.push_back(f_omegaLineshape->GetParameter(2));
			parGuesses.push_back(f_omegaLineshape->GetParameter(3));
			parGuesses.push_back(f_omegaLineshape->GetParameter(4));
			parGuesses.push_back(f_omegaLineshape->GetParameter(5));
			parGuesses.push_back(f_omegaLineshape->GetParameter(6));
			parGuesses.push_back(f_omegaLineshape->GetParameter(7));
			parGuesses.push_back(f_omegaLineshape->GetParameter(8));
			break;
		}
	}
	return;
}

void MggFitter::GuessBkgdParameters(vector<double> &parGuesses)
{
	parGuesses.clear();
	
	switch(fitOption_bkgd) {
		case 0:
			return;
		case 1:
		{
			// Polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) parGuesses.push_back(0.0);
			break;
		}
		case 2:
		{
			// Exponential:
			double p0Guess =  0.0;
			for(int ihist=0; ihist<2; ihist++) {
				double weight = ihist==0 ? 1.0 : -0.1;
				p0Guess += (weight * h_full[ihist]->GetBinContent(h_full[ihist]->FindBin(minFitRange)));
			}
			p0Guess -= f_emptyWide->Eval(minFitRange) * (h_full[0]->GetBinWidth(1) / h_empty->GetBinWidth(1));
			double p1Guess =  minFitRange;
			double p2Guess = -10.0;
			double p3Guess =  10.0;
			if(p0Guess < 0.0) p0Guess = 0.0;
			
			parGuesses.push_back(p0Guess);
			parGuesses.push_back(p1Guess);
			parGuesses.push_back(p2Guess);
			parGuesses.push_back(p3Guess);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) parGuesses.push_back(0.0001);
			break;
		}
		case 4:
		{
			// Parametric approximation to dsigma/dM for e+e- pair production
			parGuesses.push_back(0.0);
			break;
		}
	}
	return;
}

void MggFitter::GuessEtaPrimeParameters(vector<double> &parGuesses)
{
	parGuesses.clear();
	
	double NGuess     = 0.0;
	double MuGuess    = EtaAnalyzer::m_massEtap;
	double SigmaGuess = 0.035;
	
	switch(fitOption_etap) {
		case 0:
			return;
		case 1:
		{
			// Gaussian: N, mu, sigma
			parGuesses.push_back(NGuess);
			parGuesses.push_back(MuGuess);
			parGuesses.push_back(SigmaGuess);
			return;
		}
	}
	return;
}

void MggFitter::GuessEmptyParameters(vector<double> &parGuesses)
{
	parGuesses.clear();
	
	if(fitOption_empty==0) return;
	
	// get list of empty target parameters:
	vector<TString> locEmptyParameters;
	int nParsEmpty = GetEmptyParameters(locEmptyParameters);
	
	for(int ipar=0; ipar<nParsEmpty; ipar++) {
		// get the index of this parameter within the 'm_parametersFull' vector:
		int locIndex1 = (int)(find(m_parametersFull.begin(), m_parametersFull.end(), locEmptyParameters[ipar]) - m_parametersFull.begin());
		
		// check if this index is found within 'm_parIndexEmpty' and use its index as the parameter number for f_emptyWide:
		int locIndex2 = (int)(find(m_parIndexEmpty.begin(), m_parIndexEmpty.end(), locIndex1) - m_parIndexEmpty.begin());
		if(locIndex2 < (int)m_parIndexEmpty.size()) {
			
			// Set initial guess for empty target background from f_emptyWide fit result:
			parGuesses.push_back(f_emptyWide->GetParameter(locIndex2));
		}
	}
	return;
}


double MggFitter::MggFitFunction(double *x, double *par)
{
	double locMgg = x[0];
	
	// for excluding sub-regions of the fit:
	
	for(int iexc = 0; iexc < (int)excludeRegions.size(); iexc++) {
		if(excludeRegions[iexc].first < locMgg && locMgg < excludeRegions[iexc].second) {
			TF1::RejectPoint();
			return 0;
		}
	}
	
	vector<TString> dummyVec;
	
	int nParsEta    = GetSignalParameters(dummyVec);
	int nParsOmega  = GetOmegaParameters(dummyVec);
	int nParsBkgd   = GetBkgdParameters(dummyVec);
	int nParsEtap   = GetEtaPrimeParameters(dummyVec);
	int nParsEmpty  = GetEmptyParameters(dummyVec);
	
	int nParameters = nParsEta + nParsOmega + nParsBkgd + nParsEtap + nParsEmpty;
	
	double alpha_empty      = par[nParameters+0];
	double alpha_acc        = par[nParameters+1];
	double alpha_acc_switch = par[nParameters+2]; 
	// should be 1 when fitting, but can set to 0 to turn off the contribution from accidentals in the full target data, 
	// without removing accidentals from the empty target 
	
	nParameters += 3;
	
	//==================================================================================//
	// Eta:
	
	vector<double> parsEta(nParsEta+1);
	for(int ipar=0; ipar<nParsEta; ipar++) {
		parsEta[ipar] = par[ipar];
	}
	parsEta[nParsEta] = par[nParameters];
	
	double fEta = model_eta(locMgg, parsEta.data());
	
	//==================================================================================//
	// omega->pi0+gamma:
	
	vector<double> parsOmega(nParsOmega+1);
	for(int ipar=0; ipar<nParsOmega; ipar++) {
		parsOmega[ipar] = par[nParsEta + ipar];
	}
	parsOmega[nParsOmega] = par[nParameters];
	
	double fOmega = model_omega(locMgg, parsOmega.data());
	
	//==================================================================================//
	// Electromagnetic:
	
	vector<double> parsBkgd(nParsBkgd+1);
	for(int ipar=0; ipar<nParsBkgd; ipar++) {
		parsBkgd[ipar] = par[nParsEta + nParsOmega + ipar];
	}
	parsBkgd[nParsBkgd] = par[nParameters];
	
	double fBkgd = model_em(locMgg, parsBkgd.data());
	
	//==================================================================================//
	// Eta-prime:
	
	vector<double> parsEtap(nParsEtap+1);
	for(int ipar=0; ipar<nParsEtap; ipar++) {
		parsEtap[ipar] = par[nParsEta + nParsOmega + nParsBkgd + ipar];
	}
	parsEtap[nParsEtap] = par[nParameters];
	
	double fEtaPrime = model_etap(locMgg, parsEtap.data());
	
	//==================================================================================//
	// Empty target:
	
	vector<double> parsEmpty(nParsEmpty+1);
	for(int ipar=0; ipar<nParsEmpty; ipar++) {
		parsEmpty[ipar] = par[nParsEta + nParsOmega + nParsBkgd + nParsEtap + ipar];
	}
	parsEmpty[nParsEmpty] = par[nParameters];
	
	double fEmpty = alpha_empty * model_empty(locMgg, parsEmpty.data());
	
	//==================================================================================//
	// Accidentals:
	
	double fAcc = alpha_acc * alpha_acc_switch * h_full[1]->GetBinContent(h_full[1]->FindBin(locMgg));
	
	//==================================================================================//
	
	double fMgg = fEta + fOmega + fBkgd + fEtaPrime + fEmpty + fAcc;
	return fMgg;
}
