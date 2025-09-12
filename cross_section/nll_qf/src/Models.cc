#include "CrossSection.h"
#include "MggFitter.h"

double MggFitter::model_eta(double locMgg, double *par) {
	
	double fEta = 0.0;
	int nParameters = 0;
	switch(fitOption_signal) {
		case 12:
		{
			// Line shape from simulation:
			
			nParameters = 6;
			
			// (Quasi)elastic eta production:
			double N_eta   = par[0];
			double dmu     = par[1];
			double z_qf    = par[2];
			
			double fEta_exc = N_eta * (
				(1.0-z_qf)*f_cohLineshape->Eval(locMgg) +
				     z_qf * f_qfLineshape->Eval(locMgg));
			
			// Hadronic Background:
			double N_etapi   = par[3] * par[4];
			double N_etapipi = par[3] * par[5];
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double fEtaPi = N_etapi * h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-0.002)) * corrRatioEtaPi;
			
			double fEtaPiPi = N_etapipi * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-0.002)) * corrRatioBkgd;
			
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
