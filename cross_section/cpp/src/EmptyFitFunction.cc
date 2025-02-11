#include "MggFitter.h"

int MggFitter::InitializeEmptyFitFunction(TF1 **f1, TString funcName)
{
	vector<TString> parNames; parNames.clear();
	
	int nEtaPars = 0;
	switch(emptyFitOption_eta) {
		case 0:
			break;
		case 1:
			nEtaPars = 3;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#mu_{#eta}");
			parNames.push_back("#sigma_{#eta}");
			break;
		case 2:
			nEtaPars = 2;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			break;
	}
	
	int nOmegaPars = 0;
	switch(emptyFitOption_omega) {
		case 0:
			break;
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
	}
	
	int nFDCPars = 0;
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
			for(int i=0; i<m_muFDC.size(); i++) {
				nFDCPars += 3;
				parNames.push_back(Form("N_{fdc,%d}",i+1));
				parNames.push_back(Form("#mu_{fdc,%d}",i+1));
				parNames.push_back(Form("#sigma_{fdc,%d}",i+1));
			}
			break;
		case 2:
			nFDCPars += 1;
			parNames.push_back("N_{fdc}");
			for(int i=0; i<m_muFDC.size(); i++) {
				nFDCPars++;
				parNames.push_back(Form("#Delta#mu_{fdc,%d}",i+1));
			}
			break;
		case 3:
			for(int i=0; i<m_muFDC.size(); i++) {
				nFDCPars += 2;
				parNames.push_back(Form("N_{fdc,%d}",i+1));
				parNames.push_back(Form("#Delta#mu_{fdc,%d}",i+1));
			}
			break;
	}
	
	int nBkgdPars = 0;
	switch(emptyFitOption_bkgd) {
		case 1:
			nBkgdPars = emptyFitOption_poly + 1;
			break;
		case 2:
			nBkgdPars = 5;
			break;
		case 3:
			nBkgdPars = emptyFitOption_poly + 1;
			break;
		case 4:
			nBkgdPars = 0;
			break;
	}
	for(int ipar=0; ipar<nBkgdPars; ipar++) parNames.push_back(Form("p%d",ipar));
	
	int nParameters = nEtaPars + nOmegaPars + nFDCPars + nBkgdPars;
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &MggFitter::EmptyMggFitFunction, minEmptyFitRange, maxEmptyFitRange, nParameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < parNames.size(); ipar++) (*f1)->SetParName(ipar, parNames[ipar]);
	return nParameters;
}

int MggFitter::InitializeEmptyFitParameters()
{
	// Initialize parameter values to something that will not produce nans:
	
	excludeRegions.clear();
	
	// eta fit parameters:
	int nParameters = 0;
	switch(emptyFitOption_eta) {
		case 0:
			break;
		case 1:
			nParameters = 3;
			f_empty->FixParameter(0, 0.00);
			f_empty->FixParameter(1, 0.54);
			f_empty->FixParameter(2, 0.02);
			break;
		case 2:
			nParameters = 2;
			f_empty->FixParameter(0, 0.00);
			f_empty->FixParameter(1,-0.01);
			break;
	}
	
	switch(emptyFitOption_omega) {
		case 0:
			break;
		case 1:
			f_empty->FixParameter(nParameters+0, 0.00);
			f_empty->FixParameter(nParameters+1, 0.75);
			f_empty->FixParameter(nParameters+2, 0.02);
			f_empty->FixParameter(nParameters+3, 1.00);
			f_empty->FixParameter(nParameters+4, 1.00);
			nParameters += 5;
			break;
		case 2:
			f_empty->FixParameter(nParameters+0, 0.00);
			f_empty->FixParameter(nParameters+1,-0.01);
			nParameters += 2;
			break;
	}
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
			for(int i=0; i<m_muFDC.size(); i++) {
				f_empty->FixParameter(nParameters+0, 0.00);
				f_empty->FixParameter(nParameters+1, m_muFDC[i]);
				f_empty->FixParameter(nParameters+2, 0.02);
				nParameters += 3;
			}
			break;
		case 2:
			f_empty->FixParameter(nParameters, 0.00);
			nParameters++;
			for(int i=0; i<m_muFDC.size(); i++) {
				f_empty->FixParameter(nParameters, m_muFDC[i]-m_muFDC[0]);
				nParameters++;
			}
			break;
		case 3:
			for(int i=0; i<m_muFDC.size(); i++) {
				f_empty->FixParameter(nParameters+0, 0.00);
				f_empty->FixParameter(nParameters+1, m_muFDC[i]-m_muFDC[0]);
				nParameters += 2;
			}
			break;
	}
	
	// background fit parameters:
	switch(emptyFitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) f_empty->FixParameter(nParameters+ipar, 0.0);
			nParameters += (emptyFitOption_poly+1);
			break;
		case 2:
			for(int ipar=0; ipar<5; ipar++) f_empty->FixParameter(nParameters+ipar, 0.0);
			nParameters += 5;
			break;
		case 3:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) f_empty->FixParameter(nParameters+ipar, 0.0001);
			nParameters += (emptyFitOption_poly+1);
			break;
		case 4:
			nParameters += 0;
			break;
	}
	
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
	//==================================================================================//
	// Eta Mass Region:
	
	double fEta = 0.;
	int nEtaParameters = 0;
	switch(emptyFitOption_eta) {
		case 0:
		{
			// Don't try to fit the eta peak:
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			nEtaParameters = 3;
			double     N_eta = par[0];
			double    mu_eta = par[1];
			double sigma_eta = par[2];
			double     A_eta = N_eta * binSize / sqrt(2.0*TMath::Pi()) / sigma_eta;
			
			double loc_x_eta = (locMgg - mu_eta)/sigma_eta;
			fEta = A_eta * exp(-0.5*pow(loc_x_eta,2.0));
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			nEtaParameters = 2;
			
			double   N_eta = par[0];
			double dmu_eta = par[1];
			
			fEta = N_eta * f_etaLineshape->Eval(locMgg-dmu_eta);
			break;
		}
	}
	
	//==================================================================================//
	// Omega Mass Region:
	
	double fOmega = 0.;
	int nOmegaParameters = 0;
	switch(emptyFitOption_omega) {
		case 0:
		{
			// Don't try to fit the omega peak:
			break;
		}
		case 1:
		{
			nOmegaParameters = 5;
			double     N_omega = par[nEtaParameters + 0];
			double    mu_omega = par[nEtaParameters + 1];
			double sigma_omega = par[nEtaParameters + 2];
			double     a_omega = par[nEtaParameters + 3];
			double     n_omega = par[nEtaParameters + 4];
			
			double Acb_omega = pow(n_omega/fabs(a_omega), n_omega) * exp(-0.5*pow(fabs(a_omega),2.0));
			double Bcb_omega = (n_omega/fabs(a_omega)) - fabs(a_omega);
			
			double loc_x_omega = (locMgg - mu_omega)/sigma_omega;
			
			if(loc_x_omega > -a_omega) {
				fOmega = N_omega * exp(-0.5*pow(loc_x_omega,2.0));
			} else {
				fOmega = N_omega * Acb_omega * pow(Bcb_omega - loc_x_omega, -n_omega);
			}
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			
			nOmegaParameters = 2;
			
			double   N_omega = par[nEtaParameters + 0];
			double dmu_omega = par[nEtaParameters + 1];
			
			fOmega = N_omega * f_omegaLineshape->Eval(locMgg-dmu_omega);
			break;
		}
	}
	
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
				double     N_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 0];
				double    mu_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 1];
				double sigma_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 2];
				double     A_fdc = N_fdc * binSize / sqrt(2.0*TMath::Pi()) / sigma_fdc;
				
				double loc_x_fdc = (locMgg - mu_fdc)/sigma_fdc;
				fFDC += A_fdc * exp(-0.5*pow(loc_x_fdc,2.0));
				
				nFDCParameters += 3;
			}
			break;
		}
		case 2:
		{
			// Use lineshape from simulation of omegas from 1st FDC package:
			// In this option, only one normalization parameter is used for all peaking structures
			
			double N_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 0];
			nFDCParameters += 1;
			
			for(int i=0; i<m_muFDC.size(); i++) {
				double dmu_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters];
				fFDC += N_fdc * f_fdcOmegaLineshape->Eval(locMgg-dmu_fdc);
				nFDCParameters += 1;
			}
			break;
		}
		case 3:
		{
			// Use lineshape from simulation of omegas from 1st FDC package:
			// In this option, separate normalization parameters are used for all peaking structures
			
			for(int i=0; i<m_muFDC.size(); i++) {
				double   N_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 0];
				double dmu_fdc = par[nEtaParameters + nOmegaParameters + nFDCParameters + 1];
				fFDC += N_fdc * f_fdcOmegaLineshape->Eval(locMgg-dmu_fdc);
				nFDCParameters += 2;
			}
			break;
		}
	}
	
	//-------------------------//
	// Smooth background:
	
	double fBkgd = 0.;
	int nBkgdPars = 0;
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// polynomial background:
			nBkgdPars = emptyFitOption_poly+1;
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				double locPar = par[nEtaParameters + nOmegaParameters + nFDCParameters + ipar];
				double   dPar = (double)(ipar);
				fBkgd += locPar*pow(locMgg,dPar);
			}
			break;
		}
		case 2:
		{
			// exponential background:
			nBkgdPars = 5;
			
			double p0 = par[nEtaParameters + nOmegaParameters + nFDCParameters + 0];
			double p1 = par[nEtaParameters + nOmegaParameters + nFDCParameters + 1];
			double p2 = par[nEtaParameters + nOmegaParameters + nFDCParameters + 2];
			double p3 = par[nEtaParameters + nOmegaParameters + nFDCParameters + 3];
			double p4 = par[nEtaParameters + nOmegaParameters + nFDCParameters + 4];
			
			fBkgd = p0 * exp(p2*(locMgg - p4) + p3*pow(locMgg - p4,2.0));
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			nBkgdPars = emptyFitOption_poly + 1;
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				double locPar = par[nEtaParameters + nOmegaParameters + nFDCParameters + ipar];
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
	
	//----------------------------------------------------------------------------------//
	
	double fMgg = fEta + fOmega + fFDC + fBkgd;
	return fMgg;
}
