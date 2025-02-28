#include "CrossSection.h"
#include "MggFitter.h"
#include "EtaAnalyzer.h"

int MggFitter::InitializeFitFunction(TF1 **f1, TString funcName)
{
	vector<TString> parNames; parNames.clear();
	
	int nEtaPars = 0;
	switch(fitOption_signal) {
		case 1:
			nEtaPars = 3;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#mu_{#eta}");
			parNames.push_back("#sigma_{#eta}");
			break;
		case 2:
			nEtaPars = 6;
			parNames.push_back("N_{#eta}");
			parNames.push_back("fraction_{#eta}");
			parNames.push_back("#mu_{#eta,1}");
			parNames.push_back("#mu_{#eta,2}-#mu_{#eta,1}");
			parNames.push_back("#sigma_{#eta,1}");
			parNames.push_back("#sigma_{#eta,2}");
			break;
		case 3:
			nEtaPars = 5;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#mu_{#eta}");
			parNames.push_back("#sigma_{#eta}");
			parNames.push_back("#alpha_{#eta}");
			parNames.push_back("n_{#eta}");
			break;
		case 4:
			nEtaPars = 8;
			parNames.push_back("N_{#eta,1}");
			parNames.push_back("#mu_{#eta,1}");
			parNames.push_back("#sigma_{#eta,1}");
			parNames.push_back("#alpha_{#eta,1}");
			parNames.push_back("n_{#eta,1}");
			parNames.push_back("N_{#eta,2}");
			parNames.push_back("#mu_{#eta,2}");
			parNames.push_back("#sigma_{#eta,2}");
			break;
		case 5:
			nEtaPars = 2;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			break;
		case 6:
			nEtaPars = 5;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			parNames.push_back("fraction_{#eta#pi}");
			parNames.push_back("#mu_{#eta#pi}");
			parNames.push_back("#sigma_{#eta#pi}");
			break;
		case 7:
			nEtaPars = 5;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			parNames.push_back("#Delta#mu_{#eta#pi}");
			parNames.push_back("frac_{#eta}");
			parNames.push_back("frac_{#eta#pi}");
			break;
	}
	
	int nOmegaPars = 0;
	switch(fitOption_omega) {
		case 1:
			nOmegaPars = 5;
			parNames.push_back("N_{#omega}");
			parNames.push_back("#mu_{#omega}");
			parNames.push_back("#sigma_{#omega}");
			parNames.push_back("#alpha_{#omega}");
			parNames.push_back("n_{#omega}");
			break;
		case 2:
			nOmegaPars = 2;
			parNames.push_back("N_{#omega}");
			parNames.push_back("#Delta#mu_{#omega}");
			break;
		case 3:
			nOmegaPars = 2;
			parNames.push_back("N_{#omega}");
			parNames.push_back("#Delta#mu_{#omega}");
			break;
	}
	
	int nBkgdPars = 0;
	switch(fitOption_bkgd) {
		case 1:
			nBkgdPars = fitOption_poly + 1;
			break;
		case 2:
			nBkgdPars = 4;
			break;
		case 3:
			nBkgdPars = fitOption_poly + 1;
			break;
		case 4:
			nBkgdPars = 0;
			break;
	}
	for(int ipar=0; ipar<nBkgdPars; ipar++) parNames.push_back(Form("p%d",ipar));
	
	int nEtaPrimePars = 0;
	if(fitOption_etap==1) {
		nEtaPrimePars = 3;
		parNames.push_back("N_{#eta'}");
		parNames.push_back("#mu_{#eta'}");
		parNames.push_back("#sigma_{#eta'}");
	}
	
	int nParameters = nEtaPars + nOmegaPars + nBkgdPars + nEtaPrimePars;
	
	// Empty Target PDF:
	if(fitOption_empty==1) {
		parNames.push_back("N_{empty}");
		nParameters++;
	}
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &MggFitter::MggFitFunction, minFitRange, maxFitRange, nParameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < parNames.size(); ipar++) (*f1)->SetParName(ipar, parNames[ipar]);
	return nParameters;
}

int MggFitter::InitializeFitParameters()
{
	// Initialize parameter values to something that will not produce nans:
	
	excludeRegions.clear();
	
	// eta fit parameters:
	int nParameters = 0;
	switch(fitOption_signal) {
		case 1:
			nParameters = 3;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(2, 0.02);
			break;
		case 2:
			nParameters = 6;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, 0.0);
			f_fit->FixParameter(2, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(3, 0.0);
			f_fit->FixParameter(4, 0.02);
			f_fit->FixParameter(5, 0.02);
			break;
		case 3:
			nParameters = 5;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(2, 0.02);
			f_fit->FixParameter(3, 1.0);
			f_fit->FixParameter(4, 2.0);
			break;
		case 4:
			nParameters = 8;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(2, 0.02);
			f_fit->FixParameter(3, 1.0);
			f_fit->FixParameter(4, 2.0);
			f_fit->FixParameter(5, 0.0);
			f_fit->FixParameter(6, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(7, 0.02);
			break;
		case 5:
			nParameters = 2;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, 0.0);
			break;
		case 6:
			nParameters = 5;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, 0.0);
			f_fit->FixParameter(2, 0.0);
			f_fit->FixParameter(3, 0.56);
			f_fit->FixParameter(4, 0.02);
			break;
		case 7:
			nParameters = 5;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, 0.0);
			f_fit->FixParameter(2, 0.0);
			f_fit->FixParameter(3, 1.0);
			f_fit->FixParameter(4, m_etaPionFraction);
			break;
	}
	
	switch(fitOption_omega) {
		case 1:
			f_fit->FixParameter(nParameters+0, 0.0);
			f_fit->FixParameter(nParameters+1, EtaAnalyzer::m_massOmega);
			f_fit->FixParameter(nParameters+2, 0.02);
			f_fit->FixParameter(nParameters+3, 1.0);
			f_fit->FixParameter(nParameters+4, 1.0);
			nParameters += 5;
			break;
		case 2:
			f_fit->FixParameter(nParameters+0, 0.0);
			f_fit->FixParameter(nParameters+1, 0.0);
			nParameters += 2;
			break;
		case 3:
			f_fit->FixParameter(nParameters+0, 0.0);
			f_fit->FixParameter(nParameters+1, 0.0);
			nParameters += 2;
			break;
	}
	
	// background fit parameters:
	switch(fitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0);
			nParameters += (fitOption_poly+1);
			break;
		case 2:
			for(int ipar=0; ipar<4; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0);
			nParameters += 4;
			break;
		case 3:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0001);
			nParameters += (fitOption_poly+1);
			break;
		case 4:
			nParameters += 0;
			break;
	}
	
	// eta-prime fit parameters:
	f_fit->FixParameter(nParameters+0, 0.0);
	f_fit->FixParameter(nParameters+1, EtaAnalyzer::m_massEtap);
	f_fit->FixParameter(nParameters+2, 0.02);
	nParameters += 3;
	
	// A free parameter for the normalization of the empty target background:
	if(fitOption_empty==1) {
		f_fit->FixParameter(nParameters, 1.0);
		
		// Change the 'binWidth' parameter in the empty target fit function to 1, as the bin size correction
		// will be applied separately in the full fit funciton.
		// Just need to remember to re-set this parameter when plotting the empty target background.
		f_empty->SetParameter(f_empty->GetNpar()-1, 1.0);
	}
	else f_fit->FixParameter(nParameters, 0.0);
	nParameters += 1;
	
	return nParameters;
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
	
	//==================================================================================//
	// Signal:
	
	double fEta = 0.;
	
	int nSignalPars = 0;
	switch(fitOption_signal) {
		case 1:
		{
			// single Gaussian:
			
			nSignalPars = 3;
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			
			fEta = N * NormGaus(locMgg, mu, sigma);
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			nSignalPars = 6;
			double N        = par[0];
			double fraction = par[1];
			double mu1      = par[2];
			double mu2      = par[3] + mu1;
			double sigma1   = par[4];
			double sigma2   = par[5];
			
			fEta = N * ((1.0-fraction)*NormGaus(locMgg, mu1, sigma1)
				 + fraction * NormGaus(locMgg, mu2, sigma2));
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			
			nSignalPars = 5;
			
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			double alpha = par[3];
			double n     = par[4];
			
			fEta = N * NormCrystalBall(locMgg, mu, sigma, alpha, n, 1);
			break;
		}
		case 4:
		{
			// Crystal Ball + Gaussian function:
			
			nSignalPars = 8;
			
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			
			nSignalPars = 2;
			
			double N   = par[0];
			double dmu = par[1];
			fEta = N * (h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu))
				/ h_etaLineshape->GetBinWidth(1));
			break;
		}
		case 6:
		{
			// Line shape from simulation + Gaussian for eta+pion background:
			
			nSignalPars = 5;
			
			double N   = par[0];
			double dmu = par[1];
			
			double  frac_etapi = par[2];
			double    mu_etapi = par[3];
			double sigma_etapi = par[4];
			
			double f1 = h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu)) 
				/ h_etaLineshape->GetBinWidth(1);
			double f2 = NormGaus(locMgg, mu_etapi, sigma_etapi);
			
			fEta = N * (f1 + frac_etapi*f2);
			break;
		}
		case 7:
		{
			// Line shape from simulation:
			
			nSignalPars = 5;
			
			double N          = par[0];
			double dmu        = par[1];
			double dmu_etapi  = par[2] + dmu;
			double frac_eta   = par[3]; // only exists to remove signal contribution when drawing
			double frac_etapi = par[4];
			
			fEta = N * (frac_eta * f_etaLineshape->Eval(locMgg-dmu) 
				+ frac_etapi * f_etaPionLineshape->Eval(locMgg-dmu_etapi));
			break;
		}
	}
	
	//==================================================================================//
	// Background:
	
	//-------------------------//
	// omega->pi0+gamma:
	
	double fOmega = 0.;
	
	int nOmegaPars = 0;
	switch(fitOption_omega) {
		case 1:
		{
			nOmegaPars = 5;
			double N     = par[nSignalPars + 0];
			double mu    = par[nSignalPars + 1];
			double sigma = par[nSignalPars + 2];
			double alpha = par[nSignalPars + 3];
			double n     = par[nSignalPars + 4];
			
			fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			nOmegaPars = 2;
			double N   = par[nSignalPars + 0];
			double dmu = par[nSignalPars + 1];
			fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
		case 3:
		{
			nOmegaPars = 2;
			double N   = par[nSignalPars + 0];
			double dmu = par[nSignalPars + 1];
			fOmega = N * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu))
				/ h_omegaLineshape->GetBinWidth(1);
			break;
		}
	}
	
	//-------------------------//
	// electromagnetic:
	
	double fBkgd = 0.;
	
	int nBkgdPars = 0;
	switch(fitOption_bkgd) {
		case 1:
		{
			// polynomial background:
			nBkgdPars = fitOption_poly+1;
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[nSignalPars + nOmegaPars + ipar];
				double   dPar = (double)(ipar);
				fBkgd += locPar*pow(locMgg,dPar);
			}
			break;
		}
		case 2:
		{
			// exponential background:
			nBkgdPars = 4;
			
			double p0 = par[nSignalPars + nOmegaPars + 0];
			double p1 = par[nSignalPars + nOmegaPars + 1];
			double p2 = par[nSignalPars + nOmegaPars + 2];
			double p3 = par[nSignalPars + nOmegaPars + 3];
			
			fBkgd = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			nBkgdPars = fitOption_poly + 1;
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[nSignalPars + nOmegaPars + ipar];
				f_chebyshev->SetParameter(ipar, locPar);
			}
			fBkgd = f_chebyshev->Eval(locMgg);
			break;
		}
		case 4:
		{
			nBkgdPars = 0;
			break;
		}
	}
	
	//-------------------------//
	// eta-prime:
	
	double fEtaPrime = 0.0;
	
	int nEtapPars = 0;
	if(fitOption_etap==1) {
		nEtapPars = 3;
		double N     = par[nSignalPars + nOmegaPars + nBkgdPars + 0];
		double mu    = par[nSignalPars + nOmegaPars + nBkgdPars + 1];
		double sigma = par[nSignalPars + nOmegaPars + nBkgdPars + 2];
		fEtaPrime = N * NormGaus(locMgg, mu, sigma);
	}
	
	//==================================================================================//
	// Empty target:
	
	double fEmpty = 0.0;
	if(fitOption_empty==1) {
		double N = par[nSignalPars + nOmegaPars + nBkgdPars + nEtapPars];
		fEmpty = N * m_emptyRatio * f_empty->Eval(locMgg);
	}
	
	//==================================================================================//
	
	double fMgg = (fEta + fOmega + fBkgd + fEtaPrime + fEmpty) * binSize;
	return fMgg;
}
