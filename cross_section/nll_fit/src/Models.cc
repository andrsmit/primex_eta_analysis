#include "CrossSection.h"
#include "MggFitter.h"

double MggFitter::model_eta(double locMgg, double *par) {
	
	double fEta = 0.0;
	int nParameters = 0;
	switch(fitOption_signal) {
		case 1:
		{
			// Use simulated histograms as lineshape:
			
			nParameters = 6;
			
			// (Quasi)elastic eta production:
			double N_eta   = par[0];
			double dmu     = par[1];
			double z_qf    = par[2];
			
			double corrRatioCoh = 1.0 / h_cohLineshape->GetBinWidth(1);
			double corrRatioQF  = 1.0 /  h_qfLineshape->GetBinWidth(1);
			
			double fEtaCoh = N_eta * (1.0 - z_qf) * h_cohLineshape->GetBinContent(
				h_cohLineshape->FindBin(locMgg-dmu)) * corrRatioCoh;
			
			double fEtaQF = N_eta * z_qf * h_qfLineshape->GetBinContent(
				h_qfLineshape->FindBin(locMgg-dmu)) * corrRatioQF;
			
			// Non-single-eta Background:
			double N_etapi   = par[3] * par[4];
			double N_etapipi = par[3] * par[5];
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double fEtaPi = N_etapi * h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioEtaPi;
			
			double fEtaPiPi = N_etapipi * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioBkgd;
			
			fEta = fEtaCoh + fEtaQF + fEtaPi + fEtaPiPi;
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			nParameters = 6;
			
			// (Quasi)elastic eta production:
			double N_eta   = par[0];
			double dmu     = par[1];
			double z_qf    = par[2];
			
			double fEta_exc = N_eta * (
				(1.0-z_qf)*f_cohLineshape->Eval(locMgg-dmu) +
				     z_qf * f_qfLineshape->Eval(locMgg-dmu));
			
			// Hadronic Background:
			double N_etapi   = par[3] * par[4];
			double N_etapipi = par[3] * par[5];
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double fEtaPi = N_etapi * h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioEtaPi;
			
			double fEtaPiPi = N_etapipi * h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioBkgd;
			
			fEta = fEta_exc + fEtaPi + fEtaPiPi;
			break;
		}
		case 3:
		{
			// Line shape from simulation:
			
			nParameters = 5;
			
			// (Quasi)elastic eta production:
			double N_eta = par[0];
			double dmu   = par[1];
			double z_qf  = par[2];
			double f_bg  = par[3];
			double r_bg  = par[4];
			
			double corrRatioEtaPi = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioBkgd  = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double shapeEta = (
				(1.0-z_qf)*f_cohLineshape->Eval(locMgg-dmu) +
				     z_qf * f_qfLineshape->Eval(locMgg-dmu)
			);
			double shapeEtaPi = h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioEtaPi;
			
			double shapeEtaPiPi = h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioBkgd;
			
			fEta = N_eta * (
				(1.0 - f_bg) * shapeEta +
				f_bg * (r_bg * shapeEtaPi + (1.0-r_bg)*shapeEtaPiPi)
			);
			break;
		}
		case 4:
		{
			// Line shape from simulation:
			
			nParameters = 6;
			
			// (Quasi)elastic eta production:
			double N_eta = par[0];
			double dmu   = par[1];
			double z_qf  = par[2];
			double f_bg  = par[3];
			double u_bg  = par[4];
			double v_bg  = par[5];
			
			double corrRatioEtaPi   = 1.0 / h_etaPionLineshape->GetBinWidth(1);
			double corrRatioEtaPiPi = 1.0 / h_etaPiPiLineshape->GetBinWidth(1);
			double corrRatioBkgd    = 1.0 / h_hadronicBkgdLineshape->GetBinWidth(1);
			
			double shapeEta = (
				(1.0-z_qf)*f_cohLineshape->Eval(locMgg-dmu) +
				     z_qf * f_qfLineshape->Eval(locMgg-dmu)
			);
			double shapeEtaPi = h_etaPionLineshape->GetBinContent(
				h_etaPionLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioEtaPi;
			
			double shapeEtaPiPi = h_etaPiPiLineshape->GetBinContent(
				h_etaPiPiLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioEtaPiPi;
			
			double shapeBkgd = h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-bggenLSOffset)) * corrRatioBkgd;
			
			double eu = exp(u_bg);
			double ev = exp(v_bg);
			double norm = 1.0 + eu + ev;
			
			double w1 = 1.0 / norm;
			double w2 =  eu / norm;
			double w3 =  ev / norm;
			
			fEta = N_eta * (
				(1.0 - f_bg) * shapeEta +
				f_bg * (w1*shapeEtaPi + w2*shapeEtaPiPi + w3*shapeBkgd)
			);
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			
			nParameters = 6;
			
			// (Quasi)elastic eta production:
			double N_eta = par[0];
			double dmu   = par[1];
			double z_qf  = par[2];
			double f_bg  = par[3];
			double u_bg  = par[4];
			double v_bg  = par[5];
			
			double shapeEta = (
				(1.0-z_qf)*f_cohLineshape->Eval(locMgg-dmu) +
				     z_qf * f_qfLineshape->Eval(locMgg-dmu)
			);
			double shapeEtaPi   = f_etaPionLineshape->Eval(locMgg-dmu-bggenLSOffset);
			double shapeEtaPiPi = f_etaPiPiLineshape->Eval(locMgg-dmu-bggenLSOffset);
			double shapeBkgd    = h_hadronicBkgdLineshape->GetBinContent(
				h_hadronicBkgdLineshape->FindBin(locMgg-dmu-bggenLSOffset));
			
			shapeBkgd *= (1.0/h_hadronicBkgdLineshape->GetBinWidth(1));
			
			double eu = exp(u_bg);
			double ev = exp(v_bg);
			double norm = 1.0 + eu + ev;
			
			double w1 = 1.0 / norm;
			double w2 =  eu / norm;
			double w3 =  ev / norm;
			
			fEta = N_eta * (
				(1.0 - f_bg) * shapeEta +
				f_bg * (w1*shapeEtaPi + w2*shapeEtaPiPi + w3*shapeBkgd)
			);
			break;
		}
	}
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fEta*locBinWidth;
}

double MggFitter::model_omega(double locMgg, double *par)
{
	double fOmega = 0., fRho = 0.;
	int nParameters = 0;
	
	if(fitOption_rho==2) {
		// Rho0 and Omega backgrounds combined into a single lineshape
		
		double N    = par[0];
		double dmu  = par[1];
		nParameters = 2;
		
		fOmega = N * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu)) / h_omegaLineshape->GetBinWidth(1);
		double locBinWidth = par[nParameters];
		return fOmega*locBinWidth;
	}
	
	switch(fitOption_omega) {
		default:
			break;
		case 1:
		{
			// Crystal Ball function:
			
			double N     = par[0];
			double mu    = par[1];
			double sigma = par[2];
			double alpha = par[3];
			double n     = par[4];
			nParameters  = 5;
			
			fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			// Fit using simulated histogram:
			
			double N    = par[0];
			double dmu  = par[1];
			nParameters = 2;
			fOmega = N * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu)) / h_omegaLineshape->GetBinWidth(1);
			break;
		}
		case 3:
		{
			// Fit using lineshape extracted from simulated histogram:
			
			double N    = par[0];
			double dmu  = par[1];
			nParameters = 2;
			fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
		case 4:
		{
			// Fit using lineshape extracted from simulated histogram, but allow width to float:
			
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
			
			fOmega = N * ((1.0-frac)*NormCrystalBall(locMgg, mu1, sigma1, alpha1, n1) 
				+ frac*NormCrystalBall(locMgg, mu2, sigma2, alpha2, n2));
			break;
		}
		case 5:
		{
			// Fit using lineshape extracted from simulated histogram, but allow width to float:
			
			double N      = par[0];
			double dmu    = par[1];
			double xsigma = par[2];
			nParameters   = 3;
			
			double    mu1 = f_omegaLineshape->GetParameter(0) - dmu;
			double sigma1 = f_omegaLineshape->GetParameter(1) * xsigma;
			double alpha1 = f_omegaLineshape->GetParameter(2);
			double     n1 = f_omegaLineshape->GetParameter(3);
			double    mu2 = f_omegaLineshape->GetParameter(4) + mu1;
			double sigma2 = f_omegaLineshape->GetParameter(5) * xsigma;
			double alpha2 = f_omegaLineshape->GetParameter(6);
			double     n2 = f_omegaLineshape->GetParameter(7);
			double   frac = f_omegaLineshape->GetParameter(8);
			
			fOmega = N * ((1.0-frac)*NormCrystalBall(locMgg, mu1, sigma1, alpha1, n1) 
				+ frac*NormCrystalBall(locMgg, mu2, sigma2, alpha2, n2));
			break;
		}
	}
	
	switch(fitOption_rho) {
		case 0:
			break;
		case 1:
		{
			double N   = par[nParameters+0];
			double dmu = par[nParameters+1];
			nParameters += 2;
			fRho = N * h_rhoLineshape->GetBinContent(h_rhoLineshape->FindBin(locMgg-dmu)) / h_rhoLineshape->GetBinWidth(1);
			break;
		}
	}
	
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return (fOmega+fRho)*locBinWidth;
}

double MggFitter::model_em(double locMgg, double *par) {
	
	double fEM = 0.;
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
				fEM += locPar*pow(locMgg,dPar);
			}
			nParameters = (fitOption_poly+1);
			fEM = pow(fEM,2.0);
			break;
		}
		case 2:
		case 3:
		{
			// exponential background:
			double p0 = par[0];
			double p1 = par[1];
			double p2 = par[2];
			double p3 = par[3];
			nParameters = 4;
			
			fEM = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
			break;
		}
		case 4:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				double locPar = par[ipar];
				f_chebyshev->SetParameter(ipar, locPar);
			}
			nParameters = (fitOption_poly+1);
			
			fEM = f_chebyshev->Eval(locMgg);
			break;
		}
		case 5:
		{
			// Parametric approximation to dsigma/dM for e+e- pair production
			fEM = (par[0]/locMgg) * (1.0 + 2.0*pow(0.000511/locMgg,2.0)) * (6.0*log(10.0/locMgg) + 4.0*log(2.0));
			
			nParameters++;
			break;
		}
	}
	if(fEM<0.0) return -1.e6;
	
	// multiply by the bin-width of the histogram being fit, which is stored in the parameter array:
	double locBinWidth = par[nParameters];
	return fEM*locBinWidth;
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

double MggFitter::model_beamline(double locMgg, double *par) {
	
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
