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
	// The first parameter represents the central value of the angular bin we're fitting (should always remain fixed):
	
	fitter.Config().ParSettings(0).Set(m_parametersFull[0].Data(), angle);
	nParameters++;
	
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
	
	fitter.Config().ParSettings(nParameters+0).Set(m_parametersFull[nParameters+0].Data(), binSize);
	fitter.Config().ParSettings(nParameters+1).Set(m_parametersFull[nParameters+1].Data(), emptyBinSize);
	
	if(debug) printf("    done.\n");
	nParameters += (int)locEmptyPars.size();
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
	
	double NGuess     = h_full->Integral(h_full->FindBin(minEtaFit), h_full->FindBin(maxEtaFit));
	double MuGuess    = EtaAnalyzer::m_massEta;
	double SigmaGuess = 0.013;
	
	// subtract background from h_emptyWide:
	
	double nEmptyBkgd = h_emptyWide->Integral(h_emptyWide->FindBin(minEtaFit), h_emptyWide->FindBin(maxEtaFit));
	NGuess -= nEmptyBkgd;
	if(NGuess < 0.0) NGuess = 0.0;
	
	// Use peak position between 0.54-0.56 as 'mu' guess:
	double locMax = 0.0;
	for(int ibin=h_full->FindBin(0.54); ibin<=h_full->FindBin(0.56); ibin++) {
		if(h_full->GetBinContent(ibin) > locMax) {
			locMax  = h_full->GetBinContent(ibin);
			MuGuess = h_full->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_signal) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// Gaussian: N, mu, sigma
			parGuesses.push_back(NGuess);
			parGuesses.push_back(MuGuess);
			parGuesses.push_back(SigmaGuess);
			break;
		}
		case 2:
		{
			// Double Gaussian: N, frac, mu1, dmu, sigma1, sigma2
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.2);
			parGuesses.push_back(MuGuess);
			parGuesses.push_back(0.0);
			parGuesses.push_back(0.9*SigmaGuess);
			parGuesses.push_back(1.5*SigmaGuess);
			break;
		}
		case 3:
		{
			// Crsytal ball: N, mu, sigma, alpha, n
			parGuesses.push_back(NGuess);
			parGuesses.push_back(MuGuess);
			parGuesses.push_back(SigmaGuess);
			parGuesses.push_back(1.0);
			parGuesses.push_back(2.0);
			break;
		}
		case 4:
		{
			// Not implemented
			break;
		}
		case 5:
		{
			// Lineshape fit: N, deltaMu
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			break;
		}
		case 6:
		{
			// Lineshape + Crystal Ball: N, deltaMu, N_inc, mu_inc, sigma_inc, alpha_inc, n_inc
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(0.0);
			parGuesses.push_back(MuGuess+0.02);
			parGuesses.push_back(2.0*SigmaGuess);
			parGuesses.push_back(1.0);
			parGuesses.push_back(2.0);
			break;
		}
		case 7:
		{
			// Lineshape (Signal fit + hadronic bkgd fit): N, deltaMu, frac_exc, frac_inc
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(1.0);
			parGuesses.push_back(0.0);
			break;
		}
		case 8:
		{
			// Lineshape  (Signal fit + hadronic bkgd hist): N, deltaMu, frac_exc, frac_inc
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(1.0);
			parGuesses.push_back(0.0);
			break;
		}
		case 9:
		{
			// Lineshape  (Signal fit + EtaPi fit + other hadronic bkgd hist): N, deltaMu, frac_exc, frac_etapi, frac_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(1.0);
			parGuesses.push_back(0.0);
			parGuesses.push_back(0.0);
			break;
		}
		case 10:
		{
			// Lineshape  (Signal hist + EtaPi hist + other hadronic bkgd hist): N, deltaMu, frac_exc, frac_etapi, frac_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(1.0);
			parGuesses.push_back(0.0);
			parGuesses.push_back(0.0);
			break;
		}
		case 11:
		{
			// Lineshape  (Signal fit + EtaPi fit + other hadronic bkgd hist): N_exc, deltaMu, A_etapi, A_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			parGuesses.push_back(0.0);
			parGuesses.push_back(0.0);
			break;
		}
		case 12:
		{
			double locShift = 0.0025;
			//if(angle>1.5) locShift = 0.0035;
			
			// Lineshape  (Signal fit + EtaPi hist + other hadronic bkgd hist): N_exc, deltaMu, A_etapi, A_other
			parGuesses.push_back(NGuess);
			parGuesses.push_back(locShift);
			parGuesses.push_back(0.0);
			parGuesses.push_back(0.0);
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
	
	double NGuess     = h_full->Integral(h_full->FindBin(minOmegaFit), h_full->FindBin(maxOmegaFit));
	double MuGuess    = EtaAnalyzer::m_massEta;
	double SigmaGuess = 0.032;
	
	// subtract background from h_emptyWide:
	
	double nEmptyBkgd = h_emptyWide->Integral(h_emptyWide->FindBin(minOmegaFit), h_emptyWide->FindBin(maxOmegaFit));
	NGuess -= nEmptyBkgd;
	if(NGuess < 0.0) NGuess = 0.0;
	
	// Use peak position between 0.75-0.85 as 'mu' guess:
	double locMax = 0.0;
	for(int ibin=h_full->FindBin(0.75); ibin<=h_full->FindBin(0.85); ibin++) {
		if(h_full->GetBinContent(ibin) > locMax) {
			locMax  = h_full->GetBinContent(ibin);
			MuGuess = h_full->GetBinCenter(ibin);
		}
	}
	
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
			parGuesses.push_back(0.0025);
			break;
		}
		case 3:
		{
			// Lineshape hist: N, deltaMu
			parGuesses.push_back(NGuess);
			parGuesses.push_back(0.0025);
			break;
		}
		case 4:
		{
			// Lineshape hist (with rho as separate fit parameter): N_omega, N_rho, deltaMu
			parGuesses.push_back(0.88*NGuess);
			parGuesses.push_back(0.12*NGuess);
			parGuesses.push_back(0.0025);
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
			double p0Guess =  h_full->GetBinContent(h_full->FindBin(minFitRange)) - f_emptyWide->Eval(minFitRange);
			double p1Guess =  minFitRange;
			double p2Guess = -8.0;
			double p3Guess =  6.0;
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
			for(int ipar=0; ipar<=fitOption_poly; ipar++) parGuesses.push_back(0.0);
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
	
	for(int ipar=1; ipar<f_emptyWide->GetNpar()-1; ipar++) {
		parGuesses.push_back(f_emptyWide->GetParameter(ipar));
	}
	return;
}


double MggFitter::MggFitFunction(double *x, double *par)
{
	double locMgg = x[0];
	
	// for excluding sub-regions of the fit:
	
	for(int iexc = 0; iexc < excludeRegions.size(); iexc++) {
		if(excludeRegions[iexc].first < locMgg && locMgg < excludeRegions[iexc].second) {
			TF1::RejectPoint();
			return 0;
		}
	}
	
	int nParameters = 0;
	double locAngle = par[nParameters];
	nParameters += 1;
	
	//==================================================================================//
	// Signal:
	
	double fEta = 0.;
	switch(fitOption_signal) {
		case 1:
		{
			// single Gaussian:
			
			double N     = par[nParameters+0];
			double mu    = par[nParameters+1];
			double sigma = par[nParameters+2];
			nParameters += 3;
			
			fEta = N * NormGaus(locMgg, mu, sigma);
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			double N        = par[nParameters+0];
			double fraction = par[nParameters+1];
			double mu1      = par[nParameters+2];
			double mu2      = par[nParameters+3] + mu1;
			double sigma1   = par[nParameters+4];
			double sigma2   = par[nParameters+5];
			nParameters    += 6;
			
			fEta = N * ((1.0-fraction)*NormGaus(locMgg, mu1, sigma1)
				 + fraction * NormGaus(locMgg, mu2, sigma2));
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			
			double N     = par[nParameters+0];
			double mu    = par[nParameters+1];
			double sigma = par[nParameters+2];
			double alpha = par[nParameters+3];
			double n     = par[nParameters+4];
			nParameters += 5;
			
			fEta = N * NormCrystalBall(locMgg, mu, sigma, alpha, n, 1);
			break;
		}
		case 4:
		{
			// Crystal Ball + Gaussian function:
			
			nParameters += 8;
			
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			
			double N     = par[nParameters+0];
			double dmu   = par[nParameters+1];
			nParameters += 2;
			
			fEta = N * (h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu))
				/ h_etaLineshape->GetBinWidth(1));
			break;
		}
		case 6:
		{
			// Line shape from simulation + Crystal Ball for hadronic background:
			
			double N          = par[nParameters+0];
			double dmu        = par[nParameters+1];
			double fExclusive = N * f_etaLineshape->Eval(locMgg-dmu);
			
			double N_bkgd     = par[nParameters+2];
			double mu_bkgd    = par[nParameters+3];
			double sigma_bkgd = par[nParameters+4];
			double alpha_bkgd = par[nParameters+5];
			double n_bkgd     = par[nParameters+6];
			double fInclusive = N_bkgd * NormCrystalBall(locMgg, mu_bkgd, sigma_bkgd, alpha_bkgd, n_bkgd, 1);
			
			fEta = fExclusive + fInclusive;
			nParameters += 7;
			break;
		}
		case 7:
		{
			// Line shape from simulation:
			
			double N         = par[nParameters+0];
			double dmu       = par[nParameters+1];
			double frac_eta  = par[nParameters+2]; // only exists to remove signal contribution when drawing
			double frac_bkgd = par[nParameters+3];
			nParameters      += 4;
			
			fEta = N * (frac_eta * f_etaLineshape->Eval(locMgg-dmu) 
				+ frac_bkgd * f_hadronicBkgdLineshape->Eval(locMgg-dmu));
			break;
		}
		case 8:
		{
			// Line shape from simulation:
			
			double N         = par[nParameters+0];
			double dmu       = par[nParameters+1];
			double frac_eta  = par[nParameters+2]; // only exists to remove signal contribution when drawing
			double frac_bkgd = par[nParameters+3];
			nParameters      += 4;
			
			double corrRatioBkgd = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			fEta = N * (frac_eta * f_etaLineshape->Eval(locMgg-dmu) 
				+ frac_bkgd * h_hadronicBkgdLineshape->GetBinContent(h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd);
			break;
		}
		case 9:
		{
			// Line shape from simulation:
			
			double N          = par[nParameters+0];
			double dmu        = par[nParameters+1];
			double frac_eta   = par[nParameters+2]; // only exists to remove signal contribution when drawing
			double frac_etapi = par[nParameters+3];
			double frac_bkgd  = par[nParameters+4];
			nParameters       += 5;
			
			double corrRatioBkgd = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			fEta = N * (
				frac_eta   * f_etaLineshape->Eval(locMgg-dmu) + 
				frac_etapi * f_etaPionLineshape->Eval(locMgg-dmu) + 
				frac_bkgd  * h_hadronicBkgdLineshape->GetBinContent(h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd
			);
			break;
		}
		case 10:
		{
			// Line shape from simulation:
			
			double N          = par[nParameters+0];
			double dmu        = par[nParameters+1];
			double frac_eta   = par[nParameters+2]; // only exists to remove signal contribution when drawing
			double frac_etapi = par[nParameters+3];
			double frac_bkgd  = par[nParameters+4];
			nParameters       += 5;
			
			double corrRatioEta   = 1.0 / h_etaLineshape->GetBinWidth(1);
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			/*
			fEta = N * (
				frac_eta   * f_etaLineshape->Eval(locMgg-dmu) + 
				frac_etapi * h_etaPionLineshape->GetBinContent(h_etaPionLineshape->FindBin(locMgg-dmu)) * corrRatioEtaPi + 
				frac_bkgd  * h_hadronicBkgdLineshape->GetBinContent(h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd
			);
			*/
			fEta = N * (
				frac_eta   * h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu)) * corrRatioEta + 
				frac_etapi * h_etaPionLineshape->GetBinContent(h_etaPionLineshape->FindBin(locMgg-dmu)) * corrRatioEtaPi + 
				frac_bkgd  * h_hadronicBkgdLineshape->GetBinContent(h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd
			);
			break;
		}
		case 11:
		{
			// Line shape from simulation:
			
			double N_eta   = par[nParameters+0];
			double dmu     = par[nParameters+1];
			double A_etapi = par[nParameters+2];
			double A_bkgd  = par[nParameters+3];
			nParameters   += 4;
			
			double corrRatioBkgd = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			double fEta_exc = 0., fEta_pi = 0., fEta_bkgd = 0.;
			
			fEta_exc = N_eta * f_etaLineshape->Eval(locMgg-dmu);
			
			if(m_etaPionYieldBGGEN>0.1) fEta_pi = A_etapi * m_etaPionYieldBGGEN * f_etaPionLineshape->Eval(locMgg-dmu);
			
			if(m_hadronicBkgdYieldBGGEN>0.1) fEta_bkgd = A_bkgd * m_hadronicBkgdYieldBGGEN * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEta_pi + fEta_bkgd;
			break;
		}
		case 12:
		{
			// Line shape from simulation:
			
			double N_eta   = par[nParameters+0];
			double dmu     = par[nParameters+1];
			double A_etapi = par[nParameters+2];
			double A_bkgd  = par[nParameters+3];
			nParameters   += 4;
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			double fEta_exc = 0., fEta_pi = 0., fEta_bkgd = 0.;
			
			fEta_exc = N_eta * f_etaLineshape->Eval(locMgg-dmu);
			
			if(m_etaPionYieldBGGEN>0.1) fEta_pi = A_etapi * m_etaPionYieldBGGEN * h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-0.001)) * corrRatioEtaPi;
			
			if(m_hadronicBkgdYieldBGGEN>0.1) fEta_bkgd = A_bkgd * m_hadronicBkgdYieldBGGEN * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-0.001)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEta_pi + fEta_bkgd;
			break;
		}
	}
	
	//==================================================================================//
	// Background:
	
	//-------------------------//
	// omega->pi0+gamma:
	
	double fOmega = 0.;
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
		{
			double N     = par[nParameters+0];
			double mu    = par[nParameters+1];
			double sigma = par[nParameters+2];
			double alpha = par[nParameters+3];
			double n     = par[nParameters+4];
			nParameters += 5;
			
			if(N==0) fOmega = 0;
			else fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			double N     = par[nParameters+0];
			double dmu   = par[nParameters+1];
			nParameters += 2;
			
			if(N==0) fOmega = 0;
			else fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
		case 3:
		{
			double N     = par[nParameters+0];
			double dmu   = par[nParameters+1];
			nParameters += 2;
			
			if(N==0) fOmega = 0;
			else fOmega = N * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu))
				/ h_omegaLineshape->GetBinWidth(1);
			break;
		}
		case 4:
		{
			double N_omega = par[nParameters+0];
			double N_rho   = par[nParameters+1];
			double dmu     = par[nParameters+2];
			nParameters += 3;
			
			//if(N_omega==0) fOmega = 0;
			fOmega = N_omega * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu))
				/ h_omegaLineshape->GetBinWidth(1) +
					N_rho * h_rhoLineshape->GetBinContent(h_rhoLineshape->FindBin(locMgg-dmu)) / h_rhoLineshape->GetBinWidth(1);
			break;
		}
		case 5:
		{
			double N       = par[nParameters+0];
			double mu1     = par[nParameters+1];
			double sigma1  = par[nParameters+2];
			double alpha1  = par[nParameters+3];
			double n1      = par[nParameters+4];
			double mu2     = par[nParameters+5] + mu1;
			double sigma2  = par[nParameters+6];
			double alpha2  = par[nParameters+7];
			double n2      = par[nParameters+8];
			double frac    = par[nParameters+9];
			nParameters += 10;
			
			if(N==0) fOmega = 0;
			else fOmega = N * ((1.0-frac)*NormCrystalBall(locMgg, mu1, sigma1, alpha1, n1) 
				+ frac*NormCrystalBall(locMgg, mu2, sigma2, alpha2, n2));
			break;
		}
	}
	
	//-------------------------//
	// electromagnetic:
	
	double fBkgd = 0.;
	switch(fitOption_bkgd) {
		case 1:
		{
			// polynomial background:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[nParameters+ipar];
				double   dPar = (double)(ipar);
				fBkgd += locPar*pow(locMgg,dPar);
			}
			nParameters += (fitOption_poly+1);
			break;
		}
		case 2:
		{
			// exponential background:
			double p0    = par[nParameters+0];
			double p1    = par[nParameters+1];
			double p2    = par[nParameters+2];
			double p3    = par[nParameters+3];
			nParameters += 4;
			
			fBkgd = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			
			/*
			// exponential background (normalized to have unit integral within fit range):
			
			double N     = par[nParameters+0];
			double p1    = par[nParameters+1];
			double p2    = par[nParameters+2];
			double p3    = par[nParameters+3];
			nParameters += 4;
			
			double erfdiff = erf((2.0*p3*maxFitRange + (p2 - 2.0*p1*p3))/(2.0*sqrt(p3)))
				- erf((2.0*p3*minFitRange + (p2 - 2.0*p1*p3))/(2.0*sqrt(p3)));
			
			double p0 = 2.0*sqrt(p3/TMath::Pi())*exp(-(pow(p2 - 2*p1*p3,2.0)/(4.0*p3) - (p1*p1*p3 - p1*p2))) * (1.0/erfdiff);
			
			fBkgd = N * p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			*/
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[nParameters+ipar];
				f_chebyshev->SetParameter(ipar, locPar);
			}
			nParameters += (fitOption_poly+1);
			
			fBkgd = f_chebyshev->Eval(locMgg);
			break;
		}
	}
	if(fBkgd<0.0) return -1.e6;
	
	//-------------------------//
	// eta-prime:
	
	double fEtaPrime = 0.0;
	if(fitOption_etap==1) {
		double N     = par[nParameters+0];
		double mu    = par[nParameters+1];
		double sigma = par[nParameters+2];
		nParameters += 3;
		
		fEtaPrime = N * NormGaus(locMgg, mu, sigma);
	}
	
	//==================================================================================//
	// Empty target:
	
	//::::::::::::::::::::::::::::::://
	// Eta Mass Region:
	
	double fEtaEmpty = 0.;
	switch(emptyFitOption_eta) {
		case 0:
		{
			// Don't try to fit the eta peak:
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			double N     = par[nParameters+0];
			double mu    = par[nParameters+1];
			double sigma = par[nParameters+2];
			nParameters += 3;
			
			fEtaEmpty = N * NormGaus(locMgg, mu, sigma);
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			double N   = par[nParameters+0];
			double dmu = par[nParameters+1];
			nParameters += 2;
			
			fEtaEmpty = N * f_etaLineshape->Eval(locMgg-dmu);
			break;
		}
	}
	
	//::::::::::::::::::::::::::::::://
	// Omega Mass Region:
	
	double fOmegaEmpty = 0.;
	switch(emptyFitOption_omega) {
		case 0:
		{
			// Don't try to fit the omega peak:
			break;
		}
		case 1:
		{
			double N     = par[nParameters+0];
			double mu    = par[nParameters+1];
			double sigma = par[nParameters+2];
			double alpha = par[nParameters+3];
			double n     = par[nParameters+4];
			nParameters += 5;
			
			fOmegaEmpty = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			double N     = par[nParameters+0];
			double dmu   = par[nParameters+1];
			nParameters += 2;
			
			fOmegaEmpty = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
	}
	
	//::::::::::::::::::::::::::::::://
	// Smooth background:
	
	double fBkgdEmpty = 0.;
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// polynomial background:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				double locPar = par[nParameters+ipar];
				double   dPar = (double)(ipar);
				fBkgdEmpty += locPar*pow(locMgg,dPar);
			}
			nParameters += emptyFitOption_poly+1;
			break;
		}
		case 2:
		{
			// exponential background:
			double p0 = par[nParameters+0];
			double p1 = par[nParameters+1];
			double p2 = par[nParameters+2];
			double p3 = par[nParameters+3];
			nParameters += 4;
			
			fBkgdEmpty = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				double locPar = par[nParameters+ipar];
				f_chebyshev->SetParameter(ipar, locPar);
			}
			nParameters += emptyFitOption_poly+1;
			fBkgdEmpty = f_chebyshev->Eval(locMgg);
			break;
		}
		case 4:
		{
			break;
		}
	}
	//if(fBkgdEmpty<0.0) return -1.e6;
	
	//==================================================================================//
	// Peaking structures from FDC:
	//   1. omega from first package (0.61 GeV)
	//   2. omega from second package (0.54 GeV)
	//   3. eta from first package+omega from second package (0.44 GeV)
	//   
	double fFDC = 0.;
	int nFDCParameters = 0;
	switch(emptyFitOption_fdc) {
		case 0:
		{
			// Don't try to fit peaking structures:
			break;
		}
		case 1:
		{
			// Gaussian functions for each structure:
			
			for(int i=0; i<m_muFDC.size(); i++) {
				double N     = par[nParameters+0];
				double mu    = par[nParameters+1];
				double sigma = par[nParameters+2];
				fFDC += (N * NormGaus(locMgg, mu, sigma));
				nParameters += 3;
			}
			break;
		}
		case 2:
		{
			// Use lineshape from simulation of omegas from 1st FDC package:
			// In this option, only one normalization parameter is used for all peaking structures
			
			double N = par[nParameters+0];
			nParameters += 1;
			
			for(int i=0; i<m_muFDC.size(); i++) {
				double dmu = par[nParameters];
				fFDC += (N * f_fdcOmegaLineshape->Eval(locMgg-dmu));
				nParameters++;
			}
			break;
		}
		case 3:
		{
			// Use lineshape from simulation of omegas from 1st FDC package:
			// In this option, separate normalization parameters are used for all peaking structures
			
			for(int i=0; i<m_muFDC.size(); i++) {
				double N   = par[nParameters+0];
				double dmu = par[nParameters+1];
				fFDC += (N * f_fdcOmegaLineshape->Eval(locMgg-dmu));
				nParameters += 2;
			}
			break;
		}
		case 4:
		{
			break;
		}
	}
	
	double fEmpty = fEtaEmpty + fOmegaEmpty + fBkgdEmpty + fFDC;
	
	//==================================================================================//
	
	double locBinWidth = par[nParameters];
	double fMgg = (fEta + fOmega + fBkgd + fEtaPrime + fEmpty) * locBinWidth;
	return fMgg;
}
