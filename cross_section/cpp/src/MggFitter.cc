#include "MggFitter.h"
#include "EtaAnalyzer.h"
#include "CrossSection.h"

void MggFitter::DrawFitResults(TCanvas &canvas)
{
	canvas.cd();
	h_data->Draw("PE1");
	f_fit->Draw("same");
	canvas.Update();
	canvas.Modified();
	return;
}

void MggFitter::FitData()
{
	if(f_fit==NULL) InitializeFitFunction();
	InitializeFitParameters();
	
	//-----------------------------------------------------------------//
	// release omega parameters and fit the region around the peak:
	
	double minOmegaFit = 0.75;
	double maxOmegaFit = 0.85;
	
	h_data->GetXaxis()->SetRangeUser(minOmegaFit, maxOmegaFit);
	
	// get initial guesses for omega fit parameters:
	
	GuessOmegaParameters();
	
	// restrict range and fit omega parameters:
	
	f_fit->SetRange(minOmegaFit, maxOmegaFit);
	h_data->Fit(f_fit, fitOption);
	
	// fix omega fit parameters for now and focus on eta peak region:
	
	FixOmegaParameters();
	h_data->GetXaxis()->SetRangeUser(minFitRange, maxFitRange);
	
	excludeRegions.clear();
	f_fit->SetRange(0.52, 0.60);
	GuessEtaParameters();
	h_data->Fit(f_fit, fitOption);
	
	// widen fit range and release background parameters:
	
	f_fit->SetRange(minFitRange, maxFitRange);
	
	GuessBkgdParameters();
	h_data->Fit(f_fit, fitOption);
	
	// if desired fit range is within the omega peak, let those parameters float:
	
	if(maxFitRange > 0.70) {
		ReleaseOmegaParameters();
		h_data->Fit(f_fit, fitOption);
	}
	
	// try fitting eta':
	
	if(maxFitRange > EtaAnalyzer::m_massEtap) {
		GuessEtapParameters();
		h_data->Fit(f_fit, fitOption);
	}
	
	if(fitOption_signal==6) {
		
		// Use a second Gaussian to approximate shape of inclusive background:
		
		int     N_eta_par = f_fit->GetParNumber("N_{#eta,inc}");
		int    mu_eta_par = f_fit->GetParNumber("#mu_{#eta,inc}");
		int sigma_eta_par = f_fit->GetParNumber("#sigma_{#eta,inc}");
		
		f_fit->ReleaseParameter(    N_eta_par);
		f_fit->ReleaseParameter(   mu_eta_par);
		f_fit->ReleaseParameter(sigma_eta_par);
		
		f_fit->SetParameter(    N_eta_par, f_fit->GetParameter(N_eta_par)    );
		f_fit->SetParameter(   mu_eta_par, f_fit->GetParameter(mu_eta_par)   );
		f_fit->SetParameter(sigma_eta_par, f_fit->GetParameter(sigma_eta_par));
		
		f_fit->SetParLimits(    N_eta_par, 0.,   1.e5);
		f_fit->SetParLimits(   mu_eta_par, 0.55, 0.62);
		f_fit->SetParLimits(sigma_eta_par, 0.015, 0.05);
		
		h_data->Fit(f_fit, fitOption);
	}
	
	f_fit->SetRange(minFitRange, maxFitRange);
	
	return;
}

//==============================================================//
// Omega:

void MggFitter::GuessOmegaParameters()
{
	double minOmegaFit = 0.68;
	double maxOmegaFit = 0.85;
	
	double     nOmegaGuess = h_data->Integral(h_data->FindBin(minOmegaFit), h_data->FindBin(maxOmegaFit));
	double    muOmegaGuess = EtaAnalyzer::m_massOmega;
	double sigmaOmegaGuess = 0.025;
	double alphaOmegaGuess = 1.0;
	double    nnOmegaGuess = 2.0;
	
	double locOmegaMax = 0.0;
	for(int ibin=h_data->FindBin(minOmegaFit); ibin<=h_data->FindBin(maxOmegaFit); ibin++) {
		if(h_data->GetBinContent(ibin) > locOmegaMax) {
			locOmegaMax = h_data->GetBinContent(ibin);
			muOmegaGuess = h_data->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_omega) {
		case 1:
		{
			int     nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int    muOmegaPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaOmegaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaOmegaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int    nnOmegaPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->SetParameter(    nOmegaPar,     nOmegaGuess);
			f_fit->SetParameter(   muOmegaPar,    muOmegaGuess);
			f_fit->SetParameter(sigmaOmegaPar, sigmaOmegaGuess);
			f_fit->SetParameter(alphaOmegaPar, alphaOmegaGuess);
			f_fit->SetParameter(   nnOmegaPar,    nnOmegaGuess);
			
			f_fit->SetParLimits(    nOmegaPar, 0.000, 1.0e6);
			f_fit->SetParLimits(   muOmegaPar, 0.750, 0.800);
			f_fit->SetParLimits(sigmaOmegaPar, 0.015, 0.050);
			f_fit->SetParLimits(alphaOmegaPar, 0.500, 9.999);
			f_fit->SetParLimits(   nnOmegaPar, 0.500, 9.999);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(   muOmegaPar, m_OMEGA_MU);
				f_fit->FixParameter(sigmaOmegaPar, m_OMEGA_SIGMA);
				f_fit->FixParameter(alphaOmegaPar, m_OMEGA_ALPHA);
				f_fit->FixParameter(    nnOmegaPar, m_OMEGA_N);
			}
			*/
			break;
		}
		case 2:
		{
			// Using functional parameterization of MC lineshape:
			
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(  nOmegaPar, h_data->GetBinContent(h_data->FindBin(0.78))/f_omegaLineshape->Eval(0.78));
			f_fit->SetParameter(dmuOmegaPar, muOmegaGuess - h_omegaLineshape->GetBinCenter(h_omegaLineshape->GetMaximumBin()));
			
			f_fit->SetParLimits(  nOmegaPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuOmegaPar, -0.03, 0.03);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmu_omega_par, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 3:
		{
			// Lineshape fit:
			
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(  nOmegaPar, nOmegaGuess);
			f_fit->SetParameter(dmuOmegaPar, muOmegaGuess - h_omegaLineshape->GetBinCenter(h_omegaLineshape->GetMaximumBin()));
			
			f_fit->SetParLimits(  nOmegaPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuOmegaPar, -0.03, 0.03);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmu_omega_par, m_OMEGA_DMU);
			}
			*/
			break;
		}
	}
	return;
}

void MggFitter::FixOmegaParameters() 
{
	
	switch(fitOption_omega) {
		case 1:
		{
			int     nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int    muOmegaPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaOmegaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaOmegaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int    nnOmegaPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->FixParameter(    nOmegaPar, f_fit->GetParameter(    nOmegaPar));
			f_fit->FixParameter(   muOmegaPar, f_fit->GetParameter(   muOmegaPar));
			f_fit->FixParameter(sigmaOmegaPar, f_fit->GetParameter(sigmaOmegaPar));
			f_fit->FixParameter(alphaOmegaPar, f_fit->GetParameter(alphaOmegaPar));
			f_fit->FixParameter(   nnOmegaPar, f_fit->GetParameter(   nnOmegaPar));
			break;
		}
		case 2:
		{
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->FixParameter(  nOmegaPar, f_fit->GetParameter(  nOmegaPar));
			f_fit->FixParameter(dmuOmegaPar, f_fit->GetParameter(dmuOmegaPar));
			break;
		}
		case 3:
		{
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->FixParameter(  nOmegaPar, f_fit->GetParameter(  nOmegaPar));
			f_fit->FixParameter(dmuOmegaPar, f_fit->GetParameter(dmuOmegaPar));
			break;
		}
	}
	return;
}

void MggFitter::ReleaseOmegaParameters() {
	
	switch(fitOption_omega) {
		case 1:
		{
			int     nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int    muOmegaPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaOmegaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaOmegaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int    nnOmegaPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->ReleaseParameter(nOmegaPar);
			f_fit->SetParLimits(nOmegaPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(muOmegaPar);
			f_fit->SetParLimits(muOmegaPar, 0.750, 0.800);
			
			f_fit->ReleaseParameter(sigmaOmegaPar);
			f_fit->SetParLimits(sigmaOmegaPar, 0.015, 0.050);
			
			f_fit->ReleaseParameter(alphaOmegaPar);
			f_fit->SetParLimits(alphaOmegaPar, 0.500, 9.999);
			
			f_fit->ReleaseParameter(nnOmegaPar);
			f_fit->SetParLimits(nnOmegaPar, 0.500, 9.999);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(   muOmegaPar, m_OMEGA_MU);
				f_fit->FixParameter(sigmaOmegaPar, m_OMEGA_SIGMA);
				f_fit->FixParameter(alphaOmegaPar, m_OMEGA_ALPHA);
				f_fit->FixParameter(   nnOmegaPar, m_OMEGA_N);
			}
			*/
		}
		case 2:
		{
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->ReleaseParameter(nOmegaPar);
			f_fit->SetParLimits(nOmegaPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(dmuOmegaPar);
			f_fit->SetParLimits(dmuOmegaPar, -0.01, 0.01);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmu_omega_par, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 3:
		{
			int   nOmegaPar = f_fit->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->ReleaseParameter(nOmegaPar);
			f_fit->SetParLimits(nOmegaPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(dmuOmegaPar);
			f_fit->SetParLimits(dmuOmegaPar, -0.01, 0.01);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmu_omega_par, m_OMEGA_DMU);
			}
			*/
			break;
		}
	}
	return;
}

//==============================================================//
// Background:

void MggFitter::GuessBkgdParameters() {
	
	if(fitOption_bkgd==4) return;
	
	int p0Par = f_fit->GetParNumber("p0");
	int p1Par = f_fit->GetParNumber("p1");
	int p2Par = f_fit->GetParNumber("p2");
	int p3Par = f_fit->GetParNumber("p3");
	
	double p0Guess = 0.;
	double p1Guess = 0., p2Guess = 0., p3Guess = 0., p4Guess = 0.;
	
	if(fitOption_bkgd==2) {
		p0Guess = h_data->GetBinContent(h_data->FindBin(minFitRange)) - f_fit->Eval(minFitRange);
		p1Guess =  0.0;
		p2Guess = -1.0;
		p3Guess =  0.0;
		p4Guess = minFitRange;
	}
	
	switch(fitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->SetParameter(locParIndex, 0.0);
				f_fit->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p4Par = f_fit->GetParNumber("p4");
			
			f_fit->SetParameter(p0Par, p0Guess);
			f_fit->SetParameter(p1Par, p1Guess);
			f_fit->SetParameter(p2Par, p2Guess);
			//f_fit->SetParameter(p3Par, p3Guess);
			f_fit->SetParameter(p4Par, p4Guess);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e4);
			f_fit->SetParLimits(p1Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p2Par, -1.e3, 1.e3);
			//f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p4Par,  0.00, 0.50);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->SetParameter(locParIndex, 0.0001);
				f_fit->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
	}
	
	return;
}

void MggFitter::FixBkgdParameters() {
	
	if(fitOption_bkgd==4) return;
	
	switch(fitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->FixParameter(locParIndex, f_fit->GetParameter(locParIndex));
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p4_par = f_fit->GetParNumber("p4");
			for(int ipar=0; ipar<5; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->FixParameter(locParIndex, f_fit->GetParameter(locParIndex));
			}
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->FixParameter(locParIndex, f_fit->GetParameter(locParIndex));
			}
			break;
		}
	}
	
	return;
}

void MggFitter::ReleaseBkgdParameters() {
	
	if(fitOption_bkgd==4) return;
	
	switch(fitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->ReleaseParameter(locParIndex);
				f_fit->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p0Par = f_fit->GetParNumber("p0");
			int p1Par = f_fit->GetParNumber("p1");
			int p2Par = f_fit->GetParNumber("p2");
			int p3Par = f_fit->GetParNumber("p3");
			int p4Par = f_fit->GetParNumber("p4");
			
			f_fit->ReleaseParameter(p0Par);
			f_fit->ReleaseParameter(p1Par);
			f_fit->ReleaseParameter(p2Par);
			f_fit->ReleaseParameter(p3Par);
			f_fit->ReleaseParameter(p4Par);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e4);
			f_fit->SetParLimits(p1Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p2Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p4Par,  0.00, 0.50);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->ReleaseParameter(locParIndex);
				f_fit->SetParLimits(locParIndex, 0.0, 1.e4);
			}
			break;
		}
	}
	
	return;
}

//==============================================================//
// Eta-Prime:

void MggFitter::GuessEtapParameters()
{
	if(fitOption_etap==0) return;
	
	int     nEtapPar = f_fit->GetParNumber("N_{#eta'}");
	int    muEtapPar = f_fit->GetParNumber("#mu_{#eta'}");
	int sigmaEtapPar = f_fit->GetParNumber("#sigma_{#eta'}");
	
	// guess number of eta' by integrating histogram and subtracting background:
	
	double minEtapFit = 0.91;
	double maxEtapFit = 1.02;
	
	double     nEtapGuess = h_data->Integral(h_data->FindBin(minEtapFit), h_data->FindBin(maxEtapFit));
	double    muEtapGuess = EtaAnalyzer::m_massEtap;
	double sigmaEtapGuess = 0.025;
	
	f_fit->SetParameter(    nEtapPar,     nEtapGuess);
	f_fit->SetParameter(   muEtapPar,    muEtapGuess);
	f_fit->SetParameter(sigmaEtapPar, sigmaEtapGuess);
	
	f_fit->SetParLimits(    nEtapPar, 0.000, 1.0e4);
	f_fit->SetParLimits(   muEtapPar, 0.920, 0.990);
	f_fit->SetParLimits(sigmaEtapPar, 0.015, 0.050);
	return;
}

//==============================================================//
// Signal:

void MggFitter::GuessEtaParameters() {
	
	// guess number of eta' by integrating histogram and subtracting background:
	
	double minEtaFit = 0.52;
	double maxEtaFit = 0.59;
	
	double     nEtaGuess = h_data->Integral(h_data->FindBin(minEtaFit), h_data->FindBin(maxEtaFit));
	double    muEtaGuess = EtaAnalyzer::m_massEta;
	double sigmaEtaGuess = 0.015;
	
	double locEtaMax = 0.0;
	for(int ibin=h_data->FindBin(minEtaFit); ibin<=h_data->FindBin(maxEtaFit); ibin++) {
		if(h_data->GetBinContent(ibin) > locEtaMax) {
			locEtaMax  = h_data->GetBinContent(ibin);
			muEtaGuess = h_data->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_signal) {
		case 1:
		{
			// single Gaussian:
			
			int     nEtaPar = f_fit->GetParNumber("N_{#eta}");
			int    muEtaPar = f_fit->GetParNumber("#mu_{#eta}");
			int sigmaEtaPar = f_fit->GetParNumber("#sigma_{#eta}");
			
			f_fit->SetParameter(    nEtaPar,     nEtaGuess);
			f_fit->SetParameter(   muEtaPar,    muEtaGuess);
			f_fit->SetParameter(sigmaEtaPar, sigmaEtaGuess);
			
			f_fit->SetParLimits(    nEtaPar, 0.00, 1.e5);
			f_fit->SetParLimits(   muEtaPar, 0.54, 0.62);
			f_fit->SetParLimits(sigmaEtaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			int     n1EtaPar = f_fit->GetParNumber("N_{#eta,1}");
			int     n2EtaPar = f_fit->GetParNumber("N_{#eta,2}");
			int    mu1EtaPar = f_fit->GetParNumber("#mu_{#eta,1}");
			int    dmuEtaPar = f_fit->GetParNumber("#mu_{#eta,2}-#mu_{#eta,1}");
			int sigma1EtaPar = f_fit->GetParNumber("#sigma_{#eta,1}");
			int sigma2EtaPar = f_fit->GetParNumber("#sigma_{#eta,2}");
			
			f_fit->SetParameter(    n1EtaPar, 0.9*nEtaGuess);
			f_fit->SetParameter(    n2EtaPar, 0.1*nEtaGuess);
			f_fit->SetParameter(   mu1EtaPar, muEtaGuess);
			f_fit->SetParameter(   dmuEtaPar, 0.0);
			f_fit->SetParameter(sigma1EtaPar, sigmaEtaGuess);
			f_fit->SetParameter(sigma2EtaPar, 2.0*sigmaEtaGuess);
			
			f_fit->SetParLimits(    n1EtaPar,  0.00, 1.e5);
			f_fit->SetParLimits(    n2EtaPar,  0.00, 0.5*nEtaGuess);
			f_fit->SetParLimits(   mu1EtaPar,  0.54, 0.62);
			f_fit->SetParLimits(   dmuEtaPar,  0.00, 0.03);
			f_fit->SetParLimits(sigma1EtaPar,  0.01, 0.03);
			f_fit->SetParLimits(sigma2EtaPar,  0.01, 0.05);
			
			/*
			if(angle<1.0) {
				f_fit->FixParameter(    n2EtaPar, 0.0);
				f_fit->FixParameter(   dmuEtaPar, 0.0);
				f_fit->FixParameter(sigma2EtaPar, 2.0*sigmaEtaGuess);
			}
			*/
			break;
		}
		case 3:
		{
			// Crsytal ball:
			
			int     nEtaPar = f_fit->GetParNumber("N_{#eta}");
			int    muEtaPar = f_fit->GetParNumber("#mu_{#eta}");
			int sigmaEtaPar = f_fit->GetParNumber("#sigma_{#eta}");
			int alphaEtaPar = f_fit->GetParNumber("#alpha_{#eta}");
			int    nnEtaPar = f_fit->GetParNumber("n_{#eta}");
			
			f_fit->SetParameter(    nEtaPar,     nEtaGuess);
			f_fit->SetParameter(   muEtaPar,    muEtaGuess);
			f_fit->SetParameter(sigmaEtaPar, sigmaEtaGuess);
			f_fit->SetParameter(alphaEtaPar,           1.0);
			f_fit->SetParameter(   nnEtaPar,           2.0);
			
			f_fit->SetParLimits(    nEtaPar, 0.00,  1.e5);
			f_fit->SetParLimits(   muEtaPar, 0.54,  0.62);
			f_fit->SetParLimits(sigmaEtaPar, 0.01,  0.03);
			f_fit->SetParLimits(alphaEtaPar, 0.01,  9.99);
			f_fit->SetParLimits(   nnEtaPar, 0.01, 99.99);
			break;
		}
		case 4:
		{
			
			break;
		}
		case 5:
		{
			// Lineshape fit:
			
			int   nEtaPar = f_fit->GetParNumber("N_{#eta}");
			int dmuEtaPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  nEtaPar, nEtaGuess);
			f_fit->SetParameter(dmuEtaPar, 0.005);
			
			f_fit->SetParLimits(  nEtaPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuEtaPar, -0.01, 0.02);
			break;
		}
		case 6:
		{
			// Lineshape fit:
			
			int   nEtaPar = f_fit->GetParNumber("N_{#eta}");
			int dmuEtaPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  nEtaPar, nEtaGuess);
			f_fit->SetParameter(dmuEtaPar, 0.005);
			
			f_fit->SetParLimits(  nEtaPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuEtaPar, -0.01, 0.02);
			break;
		}
	}
	
	return;
}

void MggFitter::FitOmegaLineshape(int drawOption)
{
	if(fitOption_omega!=2) return;
	if(h_omegaLineshape==NULL) return;
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBall, 0.2, 0.9, 5);
	
	double     nOmegaGuess = h_omegaLineshape->Integral(h_omegaLineshape->FindBin(0.75), 
		h_omegaLineshape->FindBin(0.85));
	double    muOmegaGuess = 0.78;
	double sigmaOmegaGuess = 0.025;
	double alphaOmegaGuess = 1.0;
	double    nnOmegaGuess = 2.0;
	
	double locOmegaMax = 0.0;
	for(int ibin=h_omegaLineshape->FindBin(0.75); ibin<=h_omegaLineshape->FindBin(0.85); ibin++) {
		if(h_omegaLineshape->GetBinContent(ibin) > locOmegaMax) {
			locOmegaMax  = h_omegaLineshape->GetBinContent(ibin);
			muOmegaGuess = h_omegaLineshape->GetBinCenter(ibin);
		}
	}
	
	fOmega1->SetParameter(0,     nOmegaGuess);
	fOmega1->SetParameter(1,    muOmegaGuess);
	fOmega1->SetParameter(2, sigmaOmegaGuess);
	fOmega1->SetParameter(3, alphaOmegaGuess);
	fOmega1->SetParameter(4,    nnOmegaGuess);
	
	fOmega1->SetParLimits(0, 0.000, 1.e6);
	fOmega1->SetParLimits(1, 0.750, 0.800);
	fOmega1->SetParLimits(2, 0.015, 0.050);
	fOmega1->SetParLimits(3, 0.500, 9.999);
	fOmega1->SetParLimits(4, 0.100, 9.999);
	
	h_omegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", CrystalBall2, 0.2, 0.9, 10);
	for(int i=0; i<10; i++) {
		fOmega2->SetParameter(i, fOmega1->GetParameter(i%5));
	}
	
	fOmega2->SetParameter(5, 0.0);
	fOmega2->SetParameter(7, fOmega1->GetParameter(2)*2.0);
	
	fOmega2->SetParLimits(0, 0.000, 1.e6);
	fOmega2->SetParLimits(1, 0.750, 0.800);
	fOmega2->SetParLimits(2, 0.015, 0.050);
	fOmega2->SetParLimits(3, 0.500, 9.999);
	fOmega2->SetParLimits(4, 0.100, 9.999);
	
	fOmega2->SetParLimits(5, 0.000, 1.e6);
	fOmega2->SetParLimits(6, 0.720, 0.800);
	fOmega2->SetParLimits(7, 0.015, 0.100);
	fOmega2->SetParLimits(8, 0.500, 9.999);
	fOmega2->SetParLimits(9, 0.100, 9.999);
	
	h_omegaLineshape->Fit(fOmega2, "R0QL");
	
	if(drawOption) {
		TCanvas *cOmegaLS = new TCanvas("cOmegaLS", "cOmegaLS", 950, 700);
		cOmegaLS->SetLogy();
		h_omegaLineshape->Draw();
		fOmega2->Draw("same");
		
		TF1 *fOmegaA = new TF1("fOmegaA", CrystalBall, 0.2, 1.0, 5);
		for(int i=0; i<5; i++) fOmegaA->SetParameter(i, fOmega2->GetParameter(i));
		fOmegaA->SetLineColor(kBlue);
		fOmegaA->SetLineStyle(2);
		
		TF1 *fOmegaB = new TF1("fOmegaB", CrystalBall, 0.2, 1.0, 5);
		for(int i=0; i<5; i++) fOmegaB->SetParameter(i, fOmega2->GetParameter(i+5));
		fOmegaB->SetLineColor(kMagenta);
		fOmegaB->SetLineStyle(2);
		
		fOmegaA->Draw("same");
		fOmegaB->Draw("same");
		
		cOmegaLS->Update();
		getchar();
		delete cOmegaLS;
	}
	
	f_omegaLineshape = new TF1("f_omegaLineshape", CrystalBall2, minFitRange, maxFitRange, 10);
	f_omegaLineshape->SetParameters(fOmega2->GetParameters());
	
	fOmega1->Delete();
	fOmega2->Delete();
	
	return;
}

void MggFitter::GetYield(double &yield, double &yieldErr, int useFitPars) {
	
	// if useFitPars==0 (default), extract yield by integrating counts and subtracting bkgd parameters
	// if useFitPars==1, extract yield by integrating signal lineshape
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut);
	int maxMggBin = h_data->FindBin(maxMggCut)-1;
	
	double locMinMggCut = h_data->GetBinCenter(minMggBin) - 0.5*binSize;
	double locMaxMggCut = h_data->GetBinCenter(maxMggBin) + 0.5*binSize;
	
	// check that locMinMggCut and locMaxMggCut align with desired cut range:
	if((fabs(locMinMggCut-minMggCut)>1.e-6) || (fabs(locMaxMggCut-maxMggCut)>1.e-6)) {
		printf("\n\nWarning: Set mgg cut range does not overlap with histogram bin edges.\n");
		printf("  Cut range: %f GeV - %f GeV\n", minMggCut, maxMggCut);
		printf("  Bin edges used in signal integration: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
	}
	//-----------------------------------------------//
	
	TF1 *locfBkgd;
	InitializeFitFunction(&locfBkgd, "locBkgdClone");
	locfBkgd->SetParameters(f_fit->GetParameters());
	ZeroSignalPars(locfBkgd);
	
	if(useFitPars==0) {
		for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
			double locBkgd = locfBkgd->Eval(h_data->GetBinCenter(ibin));
			yield    += h_data->GetBinContent(ibin) - locBkgd;
			yieldErr += h_data->GetBinContent(ibin) + locBkgd;
		}
		yieldErr = sqrt(yieldErr);
	}
	else {
		switch(fitOption_signal) {
			case 1:
			{
				// Single Gaussian:
				double locN     = f_fit->GetParameter("N_{#eta}");
				double locMu    = f_fit->GetParameter("#mu_{#eta}");
				double locSigma = f_fit->GetParameter("#sigma_{#eta}");
				double locA     = locN * binSize / sqrt(2.0*TMath::Pi()) / locSigma;
				
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				yield    = IntegrateGaussian(locA, locMu, locSigma, locMinMggCut, locMaxMggCut) / binSize;
				yieldErr = f_fit->GetParError(yieldPar) * (yield / locN);
				break;
			}
			case 2:
			{
				// Double Gaussian:
				double locN1     = f_fit->GetParameter("N_{#eta,1}");
				double locMu1    = f_fit->GetParameter("#mu_{#eta,1}");
				double locSigma1 = f_fit->GetParameter("#sigma_{#eta,1}");
				double locA1     = locN1 * binSize / sqrt(2.0*TMath::Pi()) / locSigma1;
				
				double locN2     = f_fit->GetParameter("N_{#eta,2}");
				double locMu2    = f_fit->GetParameter("#mu_{#eta,2}-#mu_{#eta,1}") + locMu1;
				double locSigma2 = f_fit->GetParameter("#sigma_{#eta,2}");
				double locA2     = locN2 * binSize / sqrt(2.0*TMath::Pi()) / locSigma2;
				
				int yieldPar1 = f_fit->GetParNumber("N_{#eta,1}");
				int yieldPar2 = f_fit->GetParNumber("N_{#eta,2}");
				
				yield    = (IntegrateGaussian(locA1, locMu1, locSigma1, locMinMggCut, locMaxMggCut)
					+ IntegrateGaussian(locA2, locMu2, locSigma2, locMinMggCut, locMaxMggCut)) / binSize;
				yieldErr = sqrt(pow(f_fit->GetParError(yieldPar1),2.0) + pow(f_fit->GetParError(yieldPar1),2.0)) 
					* (yield / (locN1+locN2));
				break;
			}
			case 3:
			{
				// Crystal Ball:
				TF1 *locfSignal;
				InitializeFitFunction(&locfSignal, "locSignalClone");
				locfSignal->SetParameters(f_fit->GetParameters());
				ZeroBkgdPars(locfSignal);
				
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				yield    = locfSignal->Integral(locMinMggCut, locMaxMggCut) / binSize;
				yieldErr = (f_fit->GetParError(yieldPar)/f_fit->GetParameter(yieldPar)) * yield;
				delete locfSignal;
				break;
			}
			case 4:
			{
				// Crystal Ball + Gaussian:
				TF1 *locfSignal;
				InitializeFitFunction(&locfSignal, "locSignalClone");
				locfSignal->SetParameters(f_fit->GetParameters());
				ZeroBkgdPars(locfSignal);
				
				int yieldPar1 = f_fit->GetParNumber("N_{#eta,1}");
				int yieldPar2 = f_fit->GetParNumber("N_{#eta,2}");
				
				yield    = locfSignal->Integral(locMinMggCut, locMaxMggCut) / binSize;
				yieldErr = sqrt(pow(f_fit->GetParError(yieldPar1),2.0) + pow(f_fit->GetParError(yieldPar2),2.0))
					* (yield / (f_fit->GetParameter(yieldPar1) + f_fit->GetParameter(yieldPar2)));
				delete locfSignal;
				break;
			}
			case 5:
			{
				// Lineshape fit:
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				int locMinMggBin = h_etaLineshape->FindBin(locMinMggCut+0.5*h_etaLineshape->GetBinWidth(1));
				int locMaxMggBin = h_etaLineshape->FindBin(locMaxMggCut-0.5*h_etaLineshape->GetBinWidth(1));
				
				double lsMinMggCut = h_etaLineshape->GetBinCenter(locMinMggBin) - 0.5*h_etaLineshape->GetBinWidth(1);
				double lsMaxMggCut = h_etaLineshape->GetBinCenter(locMaxMggBin) + 0.5*h_etaLineshape->GetBinWidth(1);
				
				// check that lsMinMggCut and lsMaxMggCut align with desired cut range:
				if((fabs(lsMinMggCut-locMinMggCut)>1.e-6) || (fabs(lsMaxMggCut-locMaxMggCut)>1.e-6)) {
					printf("\n\nWarning: Mgg binning of eta lineshape does not overlap with data.\n");
					printf("  Lineshape cut range: %f GeV - %f GeV\n", lsMinMggCut, lsMaxMggCut);
					printf("  Data cut range: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
				}
				
				double binSizeRatio = h_etaLineshape->GetBinWidth(1) / binSize;
				
				yield    = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				yieldErr = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				break;
			}
			case 6:
			{
				// Lineshape + Extra Gaussian:
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				int locMinMggBin = h_etaLineshape->FindBin(locMinMggCut);
				int locMaxMggBin = h_etaLineshape->FindBin(locMaxMggCut)-1;
				
				double lsMinMggCut = h_etaLineshape->GetBinCenter(locMinMggBin) - 0.5*h_etaLineshape->GetBinWidth(1);
				double lsMaxMggCut = h_etaLineshape->GetBinCenter(locMaxMggBin) + 0.5*h_etaLineshape->GetBinWidth(1);
				
				// check that lsMinMggCut and lsMaxMggCut align with desired cut range:
				if((fabs(lsMinMggCut-locMinMggCut)>1.e-6) || (fabs(lsMaxMggCut-locMaxMggCut)>1.e-6)) {
					printf("\n\nWarning: Mgg binning of eta lineshape does not overlap with data.\n");
					printf("  Lineshape cut range: %f GeV - %f GeV\n", lsMinMggCut, lsMaxMggCut);
					printf("  Data cut range: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
				}
				
				double binSizeRatio = h_etaLineshape->GetBinWidth(1) / binSize;
				
				yield    = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				yieldErr = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				break;
			}
		}
	}
	
	delete locfBkgd;
	return;
}

void MggFitter::ZeroSignalPars(TF1 *f1)
{
	switch(fitOption_signal) {
		case 1:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 2:
			f1->SetParameter("N_{#eta,1}", 0.0);
			f1->SetParameter("N_{#eta,2}", 0.0);
			break;
		case 3:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 4:
			f1->SetParameter("N_{#eta,1}", 0.0);
			f1->SetParameter("N_{#eta,2}", 0.0);
			break;
		case 5:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 6:
			f1->SetParameter("N_{#eta}",     0.0);
			//f1->SetParameter("N_{#eta,inc}", 0.0);
			break;
	}
	return;
}

void MggFitter::ZeroBkgdPars(TF1 *f1)
{
	switch(fitOption_bkgd) {
		case 1:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
		case 2:
			for(int ipar=0; ipar<5; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
		case 3:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				f1->SetParameter(Form("p%d",ipar), 0.0);
			}
			break;
	}
	f1->SetParameter("N_{#omega}", 0.0);
	f1->SetParameter("N_{#eta'}", 0.0);
	return;
}

void MggFitter::FillPull(TH1F *h_pull)
{
	if(h_pull==NULL) {
		return;
	}
	
	// Just check that the binning of supplied histogram and h_data are consistent:
	if(h_pull->GetXaxis()->GetNbins() != h_data->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of pull histogram.\n" << endl;
	}
	
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_data->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_data->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_data->GetBinContent(ibin) - f_fit->Eval(loc_mgg))/loc_unc);
		h_pull->SetBinError(ibin, 1.0);
	}
	h_pull->GetYaxis()->SetRangeUser(-6.5, 6.5);
	
	h_pull->GetXaxis()->SetTitleSize(0.15);
	h_pull->GetXaxis()->SetLabelSize(0.10);
	h_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma}");
	h_pull->GetYaxis()->SetTitleSize(0.125);
	h_pull->GetYaxis()->SetTitleOffset(0.3);
	h_pull->GetYaxis()->SetLabelSize(0.10);
	
	return;
}

double MggFitter::IntegrateGaussian(double A, double mu, double sigma, double x1, double x2)
{
	double erf1 = TMath::Erf((x1-mu)/(sqrt(2.0)*sigma));
	double erf2 = TMath::Erf((x2-mu)/(sqrt(2.0)*sigma));
	double integral = A * sigma * sqrt(TMath::Pi()/2.0) * (erf2 - erf1);
	return integral;
}
