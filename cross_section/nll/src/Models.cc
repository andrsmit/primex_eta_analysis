#include "CrossSection.h"
#include "MggFitter.h"

double MggFitter::model_eta(double locMgg, double *par) {
	
	double fEta = 0.0;
	int nParameters = 0;
	switch(fitOption_signal) {
		case 1:
		{
			// single Gaussian:
			
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			nParameters  = 3;
			fEta = N * NormGaus(locMgg, mu, sigma);
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			double N        = par[0];
			double fraction = par[1];
			double mu1      = par[2];
			double mu2      = par[3] + mu1;
			double sigma1   = par[4];
			double sigma2   = par[5];
			nParameters     = 6;
			fEta = N * ((1.0-fraction)*NormGaus(locMgg, mu1, sigma1)
				 + fraction * NormGaus(locMgg, mu2, sigma2));
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			double alpha = par[3];
			double n     = par[4];
			nParameters  = 5;
			fEta = N * NormCrystalBall(locMgg, mu, sigma, alpha, n, 1);
			break;
		}
		case 4:
		{
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			
			double N     = par[0];
			double dmu   = par[1];
			nParameters  = 2;
			fEta = N * (h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu))
				/ h_etaLineshape->GetBinWidth(1));
			break;
		}
		case 6:
		{
			// Line shape from simulation + Crystal Ball for hadronic background:
			
			double N          = par[0];
			double dmu        = par[1];
			double fExclusive = N * f_etaLineshape->Eval(locMgg-dmu);
			
			double N_bkgd     = par[2];
			double mu_bkgd    = par[3];
			double sigma_bkgd = par[4];
			double alpha_bkgd = par[5];
			double n_bkgd     = par[6];
			nParameters       = 7;
			double fInclusive = N_bkgd * NormCrystalBall(locMgg, mu_bkgd, sigma_bkgd, alpha_bkgd, n_bkgd, 1);
			
			fEta = fExclusive + fInclusive;
			break;
		}
		case 7:
		{
			// Line shape from simulation:
			
			double N         = par[0];
			double dmu       = par[1];
			double frac_eta  = par[2]; // only exists to remove signal contribution when drawing
			double frac_bkgd = par[3];
			nParameters      = 4;
			fEta = N * (frac_eta * f_etaLineshape->Eval(locMgg-dmu) 
				+ frac_bkgd * f_hadronicBkgdLineshape->Eval(locMgg-dmu));
			break;
		}
		case 8:
		{
			// Line shape from simulation:
			
			double N         = par[0];
			double dmu       = par[1];
			double frac_eta  = par[2]; // only exists to remove signal contribution when drawing
			double frac_bkgd = par[3];
			nParameters      = 4;
			
			double corrRatioBkgd = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			fEta = N * (frac_eta * f_etaLineshape->Eval(locMgg-dmu) 
				+ frac_bkgd * h_hadronicBkgdLineshape->GetBinContent(h_hadronicBkgdLineshape->FindBin(locMgg-dmu)) * corrRatioBkgd);
			break;
		}
		case 9:
		{
			// Line shape from simulation:
			
			double N          = par[0];
			double dmu        = par[1];
			double frac_eta   = par[2]; // only exists to remove signal contribution when drawing
			double frac_etapi = par[3];
			double frac_bkgd  = par[4];
			nParameters       = 5;
			
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
			
			double N          = par[0];
			double dmu        = par[1];
			double frac_eta   = par[2]; // only exists to remove signal contribution when drawing
			double frac_etapi = par[3];
			double frac_bkgd  = par[4];
			nParameters       = 5;
			
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
			
			nParameters = 5;
			
			// (Quasi)elastic eta production:
			double N_eta   = par[0];
			double dmu     = par[1];
			double fEta_exc = N_eta * f_etaLineshape->Eval(locMgg-dmu);
			
			// Hadronic Background:
			double N_etapi   = par[2] * par[3];
			double N_etapipi = par[2] * par[4];
			
			double corrRatioBkgd = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double fEtaPi = N_etapi* f_etaPionLineshape->Eval(locMgg-0.0025);
			
			double fEtaPiPi = N_etapipi * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-0.0025)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEtaPi + fEtaPiPi;
			break;
		}
		case 12:
		{
			// Line shape from simulation:
			
			nParameters = 5;
			
			// (Quasi)elastic eta production:
			double N_eta   = par[0];
			double dmu     = par[1];
			double fEta_exc = N_eta * f_etaLineshape->Eval(locMgg-dmu);
			
			// Hadronic Background:
			double N_etapi   = par[2] * par[3];
			double N_etapipi = par[2] * par[4];
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double fEtaPi = N_etapi * h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-0.0035)) * corrRatioEtaPi;
			
			double fEtaPiPi = N_etapipi * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-0.0035)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEtaPi + fEtaPiPi;
			break;
		}
	}
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fEta*locBinWidth;
}

double MggFitter::model_omega(double locMgg, double *par) {
	
	double fOmega = 0.;
	int nParameters = 0;
	switch(fitOption_omega) {
		case 0:
			break;
		case 1:
		{
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			double alpha = par[3];
			double n     = par[4];
			nParameters  = 5;
			
			if(N==0) fOmega = 0;
			else fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			double N     = par[0];
			double dmu   = par[1];
			nParameters  = 2;
			
			if(N==0) fOmega = 0;
			else fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
		case 3:
		{
			double N     = par[0];
			double dmu   = par[1];
			nParameters  = 2;
			
			if(N==0) fOmega = 0;
			else fOmega = N * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu))
				/ h_omegaLineshape->GetBinWidth(1);
			break;
		}
		case 4:
		{
			double N_omega = par[0];
			double N_rho   = par[1];
			double dmu     = par[2];
			nParameters    = 3;
			
			//if(N_omega==0) fOmega = 0;
			fOmega = N_omega * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu))
				/ h_omegaLineshape->GetBinWidth(1) +
					N_rho * h_rhoLineshape->GetBinContent(h_rhoLineshape->FindBin(locMgg-dmu)) / h_rhoLineshape->GetBinWidth(1);
			break;
		}
		case 5:
		{
			double N       = par[0];
			double mu1     = par[1];
			double sigma1  = par[2];
			double alpha1  = par[3];
			double n1      = par[4];
			double mu2     = par[5] + mu1;
			double sigma2  = par[6];
			double alpha2  = par[7];
			double n2      = par[8];
			double frac    = par[9];
			nParameters    = 10;
			
			if(N==0) fOmega = 0;
			else fOmega = N * ((1.0-frac)*NormCrystalBall(locMgg, mu1, sigma1, alpha1, n1) 
				+ frac*NormCrystalBall(locMgg, mu2, sigma2, alpha2, n2));
			break;
		}
	}
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fOmega*locBinWidth;
}


double MggFitter::model_em(double locMgg, double *par) {
	
	double fBkgd = 0.;
	int nParameters = 0;
	switch(fitOption_bkgd) {
		case 0:
		{
			return 0.0;
		}
		case 1:
		{
			// polynomial background:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[ipar];
				double   dPar = (double)(ipar);
				fBkgd += locPar*pow(locMgg,dPar);
			}
			nParameters = (fitOption_poly+1);
			break;
		}
		case 2:
		{
			// exponential background:
			double p0 = par[0];
			double p1 = par[1];
			double p2 = par[2];
			double p3 = par[3];
			nParameters = 4;
			
			fBkgd = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[ipar];
				f_chebyshev->SetParameter(ipar, locPar);
			}
			nParameters = (fitOption_poly+1);
			
			fBkgd = f_chebyshev->Eval(locMgg);
			break;
		}
		case 4:
		{
			// Parametric approximation to dsigma/dM for e+e- pair production
			fBkgd = (par[0]/locMgg) * (1.0 + 2.0*pow(0.000511/locMgg,2.0)) * (6.0*log(10.0/locMgg) + 4.0*log(2.0));
			
			nParameters++;
			break;
		}
	}
	if(fBkgd<0.0) return -1.e6;
	
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fBkgd*locBinWidth;
}

double MggFitter::model_etap(double locMgg, double *par) {
	
	double fEtaPrime = 0.0;
	int nParameters = 0;
	switch(fitOption_etap) {
		case 0: 
			return 0.0;
		case 1:
		{
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			nParameters = 3;
			fEtaPrime = N * NormGaus(locMgg, mu, sigma);
			break;
		}
	}
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fEtaPrime*locBinWidth;
}

double MggFitter::model_empty(double locMgg, double *par) {
	
	int nParameters = 0;
	
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
	switch(emptyFitOption_fdc) {
		case 0:
		{
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
			break;
		}
		case 3:
		{
			break;
		}
		case 4:
		{
			break;
		}
	}
	double fEmpty = fBkgdEmpty + fFDC;
	
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fEmpty*locBinWidth;
}
