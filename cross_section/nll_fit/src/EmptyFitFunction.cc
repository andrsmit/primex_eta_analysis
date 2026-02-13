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
	
	// eta fit parameters:
	int nParameters = 0;
	
	// Eta + Hadronic Bkgd from residual gas:
	switch(fitOption_signal) {
		case 1:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 1.000);
			f_emptyWide->FixParameter(nParameters+5, 0.000);
			nParameters += 6;
			break;
		}
		case 2:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 1.000);
			f_emptyWide->FixParameter(nParameters+5, 0.000);
			nParameters += 6;
			break;
		}
		case 3:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 1.000);
			nParameters += 5;
			break;
		}
		case 4:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, -5.00);
			f_emptyWide->FixParameter(nParameters+5, -5.00);
			nParameters += 6;
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, -5.00);
			f_emptyWide->FixParameter(nParameters+5, -5.00);
			nParameters += 6;
			break;
		}
	}
	
	// omega from gas:
	if(fitOption_rho==2) {
		f_emptyWide->FixParameter(nParameters+0, 0.000);
		f_emptyWide->FixParameter(nParameters+1, 0.000);
		nParameters += 2;
	}
	else {
		switch(fitOption_omega) {
			default:
				break;
			case 1:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.780);
				f_emptyWide->FixParameter(nParameters+2, 0.032);
				f_emptyWide->FixParameter(nParameters+3, 1.000);
				f_emptyWide->FixParameter(nParameters+4, 2.000);
				nParameters += 5;
				break;
			}
			case 2:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.000);
				nParameters += 2;
				break;
			}
			case 3:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.000);
				nParameters += 2;
				break;
			}
			case 4:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.780);
				f_emptyWide->FixParameter(nParameters+2, 0.020);
				f_emptyWide->FixParameter(nParameters+3, 1.000);
				f_emptyWide->FixParameter(nParameters+4, 2.000);
				f_emptyWide->FixParameter(nParameters+5, 0.000);
				f_emptyWide->FixParameter(nParameters+6, 0.040);
				f_emptyWide->FixParameter(nParameters+7, 1.000);
				f_emptyWide->FixParameter(nParameters+8, 2.000);
				f_emptyWide->FixParameter(nParameters+9, 0.000);
				nParameters += 10;
				break;
			}
			case 5:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.000);
				f_emptyWide->FixParameter(nParameters+2, 1.000);
				nParameters += 3;
				break;
			}
		}
		
		switch(fitOption_rho) {
			default:
				break;
			case 1:
			{
				f_emptyWide->FixParameter(nParameters+0, 0.000);
				f_emptyWide->FixParameter(nParameters+1, 0.000);
				nParameters += 2;
				break;
			}
		}
	}
	
	// (in-target) electromagnetic fit parameters:
	switch(fitOption_bkgd) {
		case 0:
			break;
		case 1:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0);
			nParameters += (fitOption_poly+1);
			break;
		case 2:
			for(int ipar=0; ipar<4; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0);
			nParameters += 4;
			break;
		case 3:
			for(int ipar=0; ipar<4; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0);
			nParameters += 4;
			break;
		case 4:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0001);
			nParameters += (fitOption_poly+1);
			break;
		case 5:
			f_emptyWide->FixParameter(nParameters, 0.0);
			nParameters += 1;
			break;
	}
	
	// eta-prime from gas:
	switch(fitOption_etap) {
		case 0:
			break;
		case 1:
		{
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.965);
			f_emptyWide->FixParameter(nParameters+2, 0.025);
			nParameters += 3;
			break;
		}
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
	
	// A_{empty}:
	f_emptyWide->FixParameter(nParameters, 0.0);
	nParameters++;
	
	// alpha_{flux}:
	f_emptyWide->FixParameter(nParameters, m_emptyFluxRatio);
	nParameters++;
	
	// alpha_{acc}:
	f_emptyWide->FixParameter(nParameters, 0.1);
	nParameters++;
	
	// binSize_empty:
	f_emptyWide->FixParameter(nParameters, h_emptyWide[0]->GetBinWidth(1));
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
	
	vector<TString> dummyVec;
	
	int nParsEta       = GetSignalParameters(dummyVec);
	int nParsOmega     = GetOmegaParameters(dummyVec);
	int nParsEM        = GetEMParameters(dummyVec);
	int nParsEtap      = GetEtaPrimeParameters(dummyVec);
	int nParsBeamline  = GetBeamlineParameters(dummyVec);
	
	int nParameters = nParsEta + nParsOmega + nParsEM + nParsEtap + nParsBeamline;
	
	double A_empty     = par[nParameters+0];
	double alpha_flux  = par[nParameters+1];
	double alpha_acc   = par[nParameters+2];
	
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
	
	vector<double> parsEM(nParsEM+1);
	for(int ipar=0; ipar<nParsEM; ipar++) {
		parsEM[ipar] = par[nParsEta + nParsOmega + ipar];
	}
	parsEM[nParsEM] = par[nParameters];
	
	double fEM = model_em(locMgg, parsEM.data());
	
	//==================================================================================//
	// Eta-prime:
	
	vector<double> parsEtap(nParsEtap+1);
	for(int ipar=0; ipar<nParsEtap; ipar++) {
		parsEtap[ipar] = par[nParsEta + nParsOmega + nParsEM + ipar];
	}
	parsEtap[nParsEtap] = par[nParameters];
	
	double fEtaPrime = model_etap(locMgg, parsEtap.data());
	
	//==================================================================================//
	// Empty target:
	
	vector<double> parsBeamline(nParsBeamline+1);
	for(int ipar=0; ipar<nParsBeamline; ipar++) {
		parsBeamline[ipar] = par[nParsEta + nParsOmega + nParsEM + nParsEtap + ipar];
	}
	parsBeamline[nParsBeamline] = par[nParameters];
	
	double fBeamline = model_beamline(locMgg, parsBeamline.data());
	
	//==================================================================================//
	// Accidentals:
	
	double fAcc = alpha_acc * h_emptyWide[1]->GetBinContent(h_emptyWide[1]->FindBin(locMgg));
	//double fAcc = alpha_acc * h_acc_fit_empty->GetBinContent(h_acc_fit_empty->FindBin(locMgg));
	
	//==================================================================================//
	
	double fMgg = (A_empty/alpha_flux) * (fEta + fOmega + fEM + fEtaPrime) + fBeamline + fAcc;
	return fMgg;
}
