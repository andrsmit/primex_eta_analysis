#include "CrossSection.h"
#include "MggFitter.h"
#include "EtaAnalyzer.h"

int MggFitter::InitializeFitFunction(TF1 **f1, TString funcName)
{
	vector<TString> parNames; parNames.clear();
	
	int nParameters = 0;
	
	parNames.push_back("angle");
	nParameters += 1;
	
	switch(fitOption_signal) {
		case 11:
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			parNames.push_back("A_{#eta#pi}");
			parNames.push_back("A_{bkgd}");
			nParameters += 4;
			break;
		default:
			break;
	}
	
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
			parNames.push_back("N_{#omega}");
			parNames.push_back("#mu_{#omega}");
			parNames.push_back("#sigma_{#omega}");
			parNames.push_back("#alpha_{#omega}");
			parNames.push_back("n_{#omega}");
			nParameters += 5;
			break;
		case 2:
			parNames.push_back("N_{#omega}");
			parNames.push_back("#Delta#mu_{#omega}");
			nParameters += 2;
			break;
		default:
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
	nParameters += nBkgdPars;
	for(int ipar=0; ipar<nBkgdPars; ipar++) parNames.push_back(Form("p%d",ipar));
	
	// Empty Target PDF:
	if(fitOption_empty==1) {
		parNames.push_back("N_{empty}");
		nParameters += 1;
	}
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &MggFitter::MggFitFunction, -0.3, 0.60, nParameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < parNames.size(); ipar++) (*f1)->SetParName(ipar, parNames[ipar]);
	return nParameters;
}

int MggFitter::InitializeFitParameters()
{
	// Initialize parameter values to something that will not produce nans:
	
	excludeRegions.clear();
	
	int nParameters = 0;
	
	// The first parameter represents the central value of the angular bin we're fitting (should always remain fixed):
	f_fit->FixParameter(nParameters, angle);
	nParameters += 1;
	
	// eta fit parameters:
	switch(fitOption_signal) {
		case 11:
			f_fit->FixParameter(nParameters+0, 0.0);
			f_fit->FixParameter(nParameters+1, 0.000);
			f_fit->FixParameter(nParameters+2, 1.0);
			f_fit->FixParameter(nParameters+3, 1.0);
			nParameters += 4;
			break;
		default:
			break;
	}
	
	switch(fitOption_omega) {
		case 0:
			break;
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
		default:
			break;
	}
	
	// background fit parameters:
	switch(fitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0);
			nParameters += (fitOption_poly+1);
			break;
		case 2:
			//for(int ipar=0; ipar<4; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0);
			f_fit->FixParameter(nParameters+0,  0.00);
			f_fit->FixParameter(nParameters+1,  minFitRange);
			f_fit->FixParameter(nParameters+2, -6.30);
			f_fit->FixParameter(nParameters+3,  0.00);
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
	
	int nParameters = 0;
	double angle = par[nParameters];
	nParameters += 1;
	//==================================================================================//
	// Signal:
	
	double fEta = 0.;
	switch(fitOption_signal) {
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
			if(m_etaPionYieldBGGEN>0.1) fEta_pi   = A_etapi * m_etaPionYieldBGGEN * f_etaPionLineshape->Eval(locMgg-dmu);
			if(m_hadronicBkgdYieldBGGEN>0.1) fEta_bkgd = A_bkgd  * m_hadronicBkgdYieldBGGEN * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEta_pi + fEta_bkgd;
			break;
		}
		default:
			break;
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
			
			fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			double N     = par[nParameters+0];
			double dmu   = par[nParameters+1];
			nParameters += 2;
			
			fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
		default:
			break;
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
		case 4:
		{
			break;
		}
	}
	if(fBkgd<0.0) return -1.e6;
	
	//==================================================================================//
	// Empty target:
	
	double fEmpty = 0.0;
	if(fitOption_empty==1) {
		double N = par[nParameters];
		nParameters += 1;
		fEmpty = N * m_emptyRatio * f_empty->Eval(locMgg);
	}
	
	//==================================================================================//
	
	double fMgg = (fEta + fOmega + fBkgd + fEmpty) * binSize;
	return fMgg;
}
