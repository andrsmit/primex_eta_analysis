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
			parNames.push_back("N_{#eta,1}");
			parNames.push_back("N_{#eta,2}");
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
			parNames.push_back("N_{#eta,inc}");
			parNames.push_back("#mu_{#eta,inc}");
			parNames.push_back("#sigma_{#eta,inc}");
			break;
		case 7:
			nEtaPars = 4;
			parNames.push_back("N_{#eta}");
			parNames.push_back("#Delta#mu_{#eta}");
			parNames.push_back("N_{#eta#pi}");
			parNames.push_back("#Delta#mu_{#eta#pi}");
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
			nBkgdPars = 5;
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
			f_fit->FixParameter(4, 1.0);
			break;
		case 4:
			nParameters = 8;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, EtaAnalyzer::m_massEta);
			f_fit->FixParameter(2, 0.02);
			f_fit->FixParameter(3, 1.0);
			f_fit->FixParameter(4, 1.0);
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
			nParameters = 4;
			f_fit->FixParameter(0, 0.0);
			f_fit->FixParameter(1, 0.0);
			f_fit->FixParameter(2, 0.0);
			f_fit->FixParameter(3, 0.0);
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
			for(int ipar=0; ipar<5; ipar++) f_fit->FixParameter(nParameters+ipar, 0.0);
			nParameters += 5;
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
	if(fitOption_empty==1) f_fit->FixParameter(nParameters, 1.0);
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
			// double Gaussian:
			
			nSignalPars = 6;
			double     N1_eta = par[0];
			double     N2_eta = par[1];
			double    mu1_eta = par[2];
			double    mu2_eta = par[3] + mu1_eta;
			double sigma1_eta = par[4];
			double sigma2_eta = par[5];
			double     A1_eta = N1_eta * binSize / sqrt(2.0*TMath::Pi()) / sigma1_eta;
			double     A2_eta = N2_eta * binSize / sqrt(2.0*TMath::Pi()) / sigma2_eta;
			
			double loc_x1_eta = (locMgg - mu1_eta)/sigma1_eta;
			double loc_x2_eta = (locMgg - mu2_eta)/sigma2_eta;
			fEta = A1_eta*exp(-0.5*pow(loc_x1_eta,2.0)) + A2_eta*exp(-0.5*pow(loc_x2_eta,2.0));
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			
			nSignalPars = 5;
			
			double     N_eta = par[0];
			double    mu_eta = par[1];
			double sigma_eta = par[2];
			double     a_eta = par[3];
			double     n_eta = par[4];
			double     A_eta = N_eta * binSize / sqrt(2.0*TMath::Pi()) / sigma_eta;
			
			double   Acb_eta = pow(n_eta/fabs(a_eta), n_eta) * exp(-0.5*pow(fabs(a_eta),2.0));
			double   Bcb_eta = (n_eta/fabs(a_eta)) - fabs(a_eta);
			double loc_x_eta = (mu_eta - locMgg)/sigma_eta;
			if(loc_x_eta > -a_eta) {
				fEta = A_eta * exp(-0.5*pow(loc_x_eta,2.0));
			} else {
				fEta = A_eta * Acb_eta * pow(Bcb_eta - loc_x_eta, -n_eta);
			}
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
			
			double   N_eta = par[0];
			double dmu_eta = par[1];
			
			fEta = N_eta * h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu_eta));
			break;
		}
		case 6:
		{
			// Line shape from simulation:
			
			nSignalPars = 5;
			
			double   N_eta = par[0];
			double dmu_eta = par[1];
			
			double     N_eta_inc = par[2];
			double    mu_eta_inc = par[3];
			double sigma_eta_inc = par[4];
			
			fEta = N_eta * h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu_eta));
			
			double A_eta_inc = N_eta_inc * binSize / sqrt(2.0*TMath::Pi()) / sigma_eta_inc;
			double x_eta_inc = (locMgg - mu_eta_inc)/sigma_eta_inc;
			fEta += (A_eta_inc * exp(-0.5*pow(x_eta_inc,2.0)));
			break;
		}
		case 7:
		{
			// Line shape from simulation:
			
			nSignalPars = 4;
			
			double     N_eta = par[0];
			double   dmu_eta = par[1];
			double   N_etapi = par[2];
			double dmu_etapi = par[3];
			
			//fEta = (N_eta * h_etaLineshape->GetBinContent(h_etaLineshape->FindBin(locMgg-dmu_eta))
			//	+ N_etapi * f_etaPionLineshape->Eval(locMgg-dmu_eta));
			fEta = (N_eta*f_etaLineshape->Eval(locMgg-dmu_eta) + N_etapi*f_etaPionLineshape->Eval(locMgg-dmu_eta));
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
			double     N_omega = par[nSignalPars + 0];
			double    mu_omega = par[nSignalPars + 1];
			double sigma_omega = par[nSignalPars + 2];
			double     a_omega = par[nSignalPars + 3];
			double     n_omega = par[nSignalPars + 4];
			
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
			nOmegaPars = 2;
			double N_omega   = par[nSignalPars + 0];
			double dmu_omega = par[nSignalPars + 1];
			
			fOmega = N_omega * f_omegaLineshape->Eval(locMgg-dmu_omega);
			break;
		}
		case 3:
		{
			nOmegaPars = 2;
			double   N_omega = par[nSignalPars + 0];
			double dmu_omega = par[nSignalPars + 1];
			
			fOmega = N_omega * h_omegaLineshape->GetBinContent(h_omegaLineshape->FindBin(locMgg-dmu_omega));
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
			nBkgdPars = 5;
			
			double p0 = par[nSignalPars + nOmegaPars + 0];
			double p1 = par[nSignalPars + nOmegaPars + 1];
			double p2 = par[nSignalPars + nOmegaPars + 2];
			double p3 = par[nSignalPars + nOmegaPars + 3];
			double p4 = par[nSignalPars + nOmegaPars + 4];
			
			fBkgd = p0 * exp(p2*(locMgg - p4) + p3*pow(locMgg - p4,2.0));
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
		double     N_etap = par[nSignalPars + nOmegaPars + nBkgdPars + 0];
		double    mu_etap = par[nSignalPars + nOmegaPars + nBkgdPars + 1];
		double sigma_etap = par[nSignalPars + nOmegaPars + nBkgdPars + 2];
		double     A_etap = N_etap * binSize / sqrt(2.0*TMath::Pi()) / sigma_etap;
		
		fEtaPrime = A_etap*exp(-0.5*pow((locMgg-mu_etap)/sigma_etap, 2.0));
		nEtapPars = 3;
	}
	
	//==================================================================================//
	// Empty target:
	
	double fEmpty = 0.0;
	if(fitOption_empty==1) {
		double N_empty = par[nSignalPars + nOmegaPars + nBkgdPars + nEtapPars];
		fEmpty = N_empty * m_emptyRatio * f_empty->Eval(locMgg);
	}
	
	//==================================================================================//
	
	double fMgg = fEta + fOmega + fBkgd + fEtaPrime + fEmpty;
	return fMgg;
}
