#include "CrossSection.h"
#include "MggFitter.h"

int MggFitter::InitializeEmptyFitFunction(TF1 **f1, TString funcName)
{
	vector<TString> parNames; parNames.clear();
	
	int nParameters = 0;
	
	for(int ipar=0; ipar<(int)m_parIndexEmpty.size(); ipar++) {
		TString locParName = m_parametersFull[m_parIndexEmpty[ipar]];
		parNames.push_back(locParName);
		nParameters++;
	}
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &MggFitter::EmptyMggFitFunction, minEmptyFitRange, maxEmptyFitRange, nParameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < nParameters; ipar++) (*f1)->SetParName(ipar, parNames[ipar]);
	return nParameters;
}

int MggFitter::InitializeEmptyWideFitParameters()
{
	// Initialize parameter values to something that will not produce nans:
	
	excludeRegions.clear();
	
	f_emptyWide->FixParameter(0, angle);
	
	// eta fit parameters:
	int nParameters = 1;
	
	// eta from gas:
	switch(emptyFitOption_eta) {
		case 0:
			break;
		case 1:
			f_emptyWide->FixParameter(nParameters+0, 0.00);
			f_emptyWide->FixParameter(nParameters+1, 0.54);
			f_emptyWide->FixParameter(nParameters+2, 0.02);
			nParameters += 3;
			break;
		case 2:
			f_emptyWide->FixParameter(nParameters+0, 0.0000);
			f_emptyWide->FixParameter(nParameters+1, 0.0025);
			nParameters += 2;
			break;
	}
	
	// omega from gas:
	switch(emptyFitOption_omega) {
		case 0:
			break;
		case 1:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.760);
			f_emptyWide->FixParameter(nParameters+2, 0.032);
			f_emptyWide->FixParameter(nParameters+3, 1.000);
			f_emptyWide->FixParameter(nParameters+4, 2.000);
			nParameters += 5;
			break;
		case 2:
			f_emptyWide->FixParameter(nParameters+0, 0.00);
			f_emptyWide->FixParameter(nParameters+1,-0.01);
			nParameters += 2;
			break;
	}
	
	// smooth background fit parameters:
	switch(emptyFitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0);
			nParameters += (emptyFitOption_poly+1);
			break;
		case 2:
			for(int ipar=0; ipar<4; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0);
			nParameters += 4;
			break;
		case 3:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0001);
			nParameters += (emptyFitOption_poly+1);
			break;
		case 4:
			nParameters += 0;
			break;
	}
	
	// peaking structures from FDC packages:
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
			for(int i=0; i<m_muFDC.size(); i++) {
				f_emptyWide->FixParameter(nParameters+0, 0.00);
				f_emptyWide->FixParameter(nParameters+1, m_muFDC[i]);
				f_emptyWide->FixParameter(nParameters+2, 0.02);
				nParameters += 3;
			}
			break;
		case 2:
			f_emptyWide->FixParameter(nParameters, 0.00);
			nParameters++;
			for(int i=0; i<m_muFDC.size(); i++) {
				f_emptyWide->FixParameter(nParameters, m_muFDC[i]-m_muFDC[0]);
				nParameters++;
			}
			break;
		case 3:
			for(int i=0; i<m_muFDC.size(); i++) {
				f_emptyWide->FixParameter(nParameters+0, 0.00);
				f_emptyWide->FixParameter(nParameters+1, m_muFDC[i]-m_muFDC[0]);
				nParameters += 2;
			}
			break;
		case 4:
			for(int i=0; i<m_muFDC_omega.size(); i++) {
				f_emptyWide->FixParameter(nParameters+0, 0.00);
				f_emptyWide->FixParameter(nParameters+1, m_muFDC_omega[i]);
				nParameters += 2;
			}
			for(int i=0; i<m_muFDC_eta.size(); i++) {
				f_emptyWide->FixParameter(nParameters+0, 0.00);
				f_emptyWide->FixParameter(nParameters+1, m_muFDC_eta[i]);
				nParameters += 2;
			}
			break;
	}
	
	f_emptyWide->FixParameter(nParameters, h_emptyWide->GetXaxis()->GetBinWidth(1));
	nParameters++;
	
	return nParameters;
}

double MggFitter::EmptyMggFitFunction(double *x, double *par)
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
	// Eta Mass Region:
	
	double fEta = 0.;
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
			
			fEta = N * NormGaus(locMgg, mu, sigma);
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			double N   = par[nParameters+0];
			double dmu = par[nParameters+1];
			nParameters += 2;
			
			fEta = N * f_etaLineshape->Eval(locMgg-dmu);
			break;
		}
	}
	
	//==================================================================================//
	// Omega Mass Region:
	
	double fOmega = 0.;
	switch(emptyFitOption_omega) {
		case 0:
		{
			// Don't try to fit the omega peak:
			break;
		}
		case 1:
		{
			double N     = par[nParameters + 0];
			double mu    = par[nParameters + 1];
			double sigma = par[nParameters + 2];
			double alpha = par[nParameters + 3];
			double n     = par[nParameters + 4];
			nParameters += 5;
			
			fOmega = N * NormCrystalBall(locMgg, mu, sigma, alpha, n);
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			double N     = par[nParameters + 0];
			double dmu   = par[nParameters + 1];
			nParameters += 2;
			
			fOmega = N * f_omegaLineshape->Eval(locMgg-dmu);
			break;
		}
	}
	
	//==================================================================================//
	// Smooth background:
	
	double fBkgd = 0.;
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// polynomial background:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				double locPar = par[nParameters+ipar];
				double   dPar = (double)(ipar);
				fBkgd += locPar*pow(locMgg,dPar);
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
			
			fBkgd = p0 * exp(p2*(locMgg - p1) + p3*pow(locMgg - p1,2.0));
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
			fBkgd = f_chebyshev->Eval(locMgg);
			break;
		}
		case 4:
		{
			break;
		}
	}
	//if(fBkgd<0.0) return -1.e6;
	
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
	
	//----------------------------------------------------------------------------------//
	
	double locBinWidth = par[nParameters];
	double fMgg = (fEta + fOmega + fFDC + fBkgd) * locBinWidth;
	return fMgg;
}
