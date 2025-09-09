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
			// single Gaussian:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.547);
			f_emptyWide->FixParameter(nParameters+2, 0.013);
			nParameters += 3;
			break;
		}
		case 2:
		{
			// double Gaussian:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.547);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 0.010);
			f_emptyWide->FixParameter(nParameters+5, 0.020);
			nParameters += 6;
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.547);
			f_emptyWide->FixParameter(nParameters+2, 0.013);
			f_emptyWide->FixParameter(nParameters+3, 1.000);
			f_emptyWide->FixParameter(nParameters+4, 2.000);
			nParameters += 5;
			break;
		}
		case 4:
		{
			break;
		}
		case 5:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			nParameters += 2;
			break;
		}
		case 6:
		{
			// Line shape from simulation + Crystal Ball for hadronic background:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 0.580);
			f_emptyWide->FixParameter(nParameters+4, 0.020);
			f_emptyWide->FixParameter(nParameters+5, 1.000);
			f_emptyWide->FixParameter(nParameters+6, 2.000);
			nParameters += 7;
			break;
		}
		case 7:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 1.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			nParameters += 4;
			break;
		}
		case 8:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 1.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			nParameters += 4;
			break;
		}
		case 9:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 1.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 0.000);
			nParameters += 5;
			break;
		}
		case 10:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 1.000);
			f_emptyWide->FixParameter(nParameters+3, 0.000);
			f_emptyWide->FixParameter(nParameters+4, 0.000);
			nParameters += 5;
			break;
		}
		case 11:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 1.000);
			f_emptyWide->FixParameter(nParameters+4, 0.000);
			nParameters += 5;
			break;
		}
		case 12:
		{
			// Line shape from simulation:
			f_emptyWide->FixParameter(nParameters+0, 0.000);
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			f_emptyWide->FixParameter(nParameters+3, 1.000);
			f_emptyWide->FixParameter(nParameters+4, 0.000);
			nParameters += 5;
			break;
		}
	}
	
	// omega from gas:
	switch(fitOption_omega) {
		case 0:
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
			f_emptyWide->FixParameter(nParameters+1, 0.000);
			f_emptyWide->FixParameter(nParameters+2, 0.000);
			nParameters += 3;
			break;
		}
		case 5:
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
			for(int ipar=0; ipar<=fitOption_poly; ipar++) f_emptyWide->FixParameter(nParameters+ipar, 0.0001);
			nParameters += (fitOption_poly+1);
			break;
		case 4:
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
	
	// binSize_empty:
	f_emptyWide->FixParameter(nParameters, h_emptyWide->GetBinWidth(1));
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
	
	int nParsEta    = GetSignalParameters(dummyVec);
	int nParsOmega  = GetOmegaParameters(dummyVec);
	int nParsBkgd   = GetBkgdParameters(dummyVec);
	int nParsEtap   = GetEtaPrimeParameters(dummyVec);
	int nParsEmpty  = GetEmptyParameters(dummyVec);
	
	int nParameters = nParsEta + nParsOmega + nParsBkgd + nParsEtap + nParsEmpty;
	
	double A_empty  = par[nParameters+0];
	nParameters++;
	
	//printf("\n\nEMPTY BIN SIZE: %f\n", par[nParameters]);
	
	//==================================================================================//
	// Eta:
	
	vector<double> parsEta(nParsEta+1);
	for(int ipar=0; ipar<nParsEta; ipar++) {
		parsEta[ipar] = par[ipar];
	}
	parsEta[nParsEta] = par[nParameters];
	
	double fEta = A_empty * model_eta(locMgg, parsEta.data());
	
	//==================================================================================//
	// omega->pi0+gamma:
	
	vector<double> parsOmega(nParsOmega+1);
	for(int ipar=0; ipar<nParsOmega; ipar++) {
		parsOmega[ipar] = par[nParsEta + ipar];
	}
	parsOmega[nParsOmega] = par[nParameters];
	
	double fOmega = A_empty * model_omega(locMgg, parsOmega.data());
	
	//==================================================================================//
	// Electromagnetic:
	
	vector<double> parsBkgd(nParsBkgd+1);
	for(int ipar=0; ipar<nParsBkgd; ipar++) {
		parsBkgd[ipar] = par[nParsEta + nParsOmega + ipar];
	}
	parsBkgd[nParsBkgd] = par[nParameters];
	
	double fBkgd = A_empty * model_em(locMgg, parsBkgd.data());
	
	//==================================================================================//
	// Eta-prime:
	
	vector<double> parsEtap(nParsEtap+1);
	for(int ipar=0; ipar<nParsEtap; ipar++) {
		parsEtap[ipar] = par[nParsEta + nParsOmega + nParsBkgd + ipar];
	}
	parsEtap[nParsEtap] = par[nParameters];
	
	double fEtaPrime = A_empty * model_etap(locMgg, parsEtap.data());
	
	//==================================================================================//
	// Empty target:
	
	vector<double> parsEmpty(nParsEmpty+1);
	for(int ipar=0; ipar<nParsEmpty; ipar++) {
		parsEmpty[ipar] = par[nParsEta + nParsOmega + nParsBkgd + nParsEtap + ipar];
	}
	parsEmpty[nParsEmpty] = par[nParameters];
	
	double fEmpty = model_empty(locMgg, parsEmpty.data());
	
	//==================================================================================//
	
	double fMgg = fEta + fOmega + fBkgd + fEtaPrime + fEmpty;
	return fMgg;
}
