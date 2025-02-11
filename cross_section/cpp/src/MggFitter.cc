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
	InitializeFitFunction(&f_fit);
	m_nParameters = InitializeFitParameters();
	
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
	
	// if we're fitting with the empty target pdf, release that normalization here:
	if(fitOption_empty) {
		//if(emptyFitOption_eta>1)   f_empty->SetParameter("N_{#eta}",   0.0);
		//if(emptyFitOption_omega>1) f_empty->SetParameter("N_{#omega}", 0.0);
		
		GuessEmptyParameters();
		h_data->Fit(f_fit, fitOption);
	}
	
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
	
	if(fitOption_signal==6 && angle>0.5) {
		
		// Use a second Gaussian to approximate shape of inclusive background:
		
		int     N_eta_par = f_fit->GetParNumber("N_{#eta,inc}");
		int    mu_eta_par = f_fit->GetParNumber("#mu_{#eta,inc}");
		int sigma_eta_par = f_fit->GetParNumber("#sigma_{#eta,inc}");
		
		f_fit->ReleaseParameter(    N_eta_par);
		f_fit->ReleaseParameter(   mu_eta_par);
		f_fit->ReleaseParameter(sigma_eta_par);
		
		f_fit->SetParLimits(    N_eta_par, 0.,   1.e5);
		f_fit->SetParLimits(   mu_eta_par, 0.55, 0.62);
		f_fit->SetParLimits(sigma_eta_par, 0.015, 0.05);
		
		h_data->Fit(f_fit, fitOption);
	}
	
	if(fitOption_signal==7) {
		
		FixBkgdParameters();
		
		// Use a eta+pion lineshape from bggen:
		f_fit->SetRange(0.5, 0.62);
		
		int N_etapi_par = f_fit->GetParNumber("N_{#eta#pi}");
		f_fit->ReleaseParameter(N_etapi_par);
		f_fit->SetParLimits(N_etapi_par, 0., 1.e5);
		
		h_data->Fit(f_fit, fitOption);
		
		f_fit->SetRange(minFitRange, maxFitRange);
		ReleaseBkgdParameters();
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
			f_fit->SetParameter(p3Par, p3Guess);
			f_fit->SetParameter(p4Par, p4Guess);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e4);
			f_fit->SetParLimits(p1Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p2Par, -1.e3, 0.);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
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
			f_fit->SetParLimits(p2Par, -1.e3, 0.);
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
				f_fit->SetParLimits(locParIndex, -1.e4, 1.e4);
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
	
	double minEtaFit = 0.50;
	double maxEtaFit = 0.60;
	
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
			
			if(angle<0.5) {
				f_fit->FixParameter(    n2EtaPar, 0.0);
				f_fit->FixParameter(   dmuEtaPar, 0.0);
				f_fit->FixParameter(sigma2EtaPar, 2.0*sigmaEtaGuess);
			}
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
			// Start with Lineshape fit:
			
			int   nEtaPar = f_fit->GetParNumber("N_{#eta}");
			int dmuEtaPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  nEtaPar, nEtaGuess);
			f_fit->SetParameter(dmuEtaPar, 0.005);
			
			f_fit->SetParLimits(  nEtaPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuEtaPar, -0.01, 0.02);
			break;
		}
		case 7:
		{
			// Start with Lineshape fit:
			
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


void MggFitter::GuessEmptyParameters()
{
	if(fitOption_empty!=1) return;
	
	int nEmptyPar = f_fit->GetParNumber("N_{empty}");
	f_fit->SetParameter(nEmptyPar, m_emptyRatio);
	f_fit->SetParLimits(nEmptyPar, m_emptyRatio-m_emptyRatioErr, m_emptyRatio+m_emptyRatioErr);
	
	return;
}

void MggFitter::FixEmptyParameters()
{
	if(fitOption_empty!=1) return;
	
	int nEmptyPar = f_fit->GetParNumber("N_{empty}");
	f_fit->FixParameter(nEmptyPar, f_fit->GetParameter(nEmptyPar));
	return;
}

void MggFitter::FitEtaLineshape(int drawOption)
{
	if(h_etaLineshape==NULL) return;
	
	
	TF1 *fEta1 = new TF1("fEta1", CrystalBall, 0.4, 0.65, 5);
	
	double     nEtaGuess = h_etaLineshape->Integral(h_etaLineshape->FindBin(0.5), 
		h_etaLineshape->FindBin(0.6));
	double    muEtaGuess = 0.547;
	double sigmaEtaGuess = 0.015;
	double alphaEtaGuess = 1.0;
	double    nnEtaGuess = 2.0;
	
	fEta1->SetParameter(0,     nEtaGuess);
	fEta1->SetParameter(1,    muEtaGuess);
	fEta1->SetParameter(2, sigmaEtaGuess);
	fEta1->SetParameter(3, alphaEtaGuess);
	fEta1->SetParameter(4,    nnEtaGuess);
	
	fEta1->SetParLimits(0, 0.000, 1.e6);
	fEta1->SetParLimits(1, 0.530, 0.560);
	fEta1->SetParLimits(2, 0.005, 0.050);
	fEta1->SetParLimits(3, 0.500, 9.999);
	fEta1->SetParLimits(4, 0.100, 9.999);
	
	h_etaLineshape->Fit(fEta1, "R0QL");
	
	TF1 *fEta2 = new TF1("fEta2", CrystalBall2, 0.4, 0.65, 10);
	for(int i=0; i<10; i++) {
		fEta2->SetParameter(i, fEta1->GetParameter(i%5));
	}
	
	fEta2->SetParameter(5, 0.0);
	fEta2->SetParameter(7, fEta1->GetParameter(2)*2.0);
	
	fEta2->SetParLimits(0, 0.000, 1.e6);
	fEta2->SetParLimits(1, 0.530, 0.560);
	fEta2->SetParLimits(2, 0.005, 0.050);
	fEta2->SetParLimits(3, 0.500, 9.999);
	fEta2->SetParLimits(4, 0.100, 9.999);
	
	fEta2->SetParLimits(5, 0.000, 1.e6);
	fEta2->SetParLimits(6, 0.500, 0.600);
	fEta2->SetParLimits(7, 0.010, 0.100);
	fEta2->SetParLimits(8, 0.500, 9.999);
	fEta2->SetParLimits(9, 0.100, 9.999);
	
	h_etaLineshape->Fit(fEta2, "R0QL");
	
	if(drawOption) {
		TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
		//cEtaLS->SetLogy();
		h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
		h_etaLineshape->Draw();
		fEta2->Draw("same");
		
		TF1 *fEtaA = new TF1("fEtaA", CrystalBall, 0.4, 0.7, 5);
		for(int i=0; i<5; i++) fEtaA->SetParameter(i, fEta2->GetParameter(i));
		fEtaA->SetLineColor(kBlue);
		fEtaA->SetLineStyle(2);
		
		TF1 *fEtaB = new TF1("fEtaB", CrystalBall, 0.4, 0.7, 5);
		for(int i=0; i<5; i++) fEtaB->SetParameter(i, fEta2->GetParameter(i+5));
		fEtaB->SetLineColor(kMagenta);
		fEtaB->SetLineStyle(2);
		
		fEtaA->Draw("same");
		fEtaB->Draw("same");
		
		cEtaLS->Update();
		getchar();
		delete cEtaLS;
	}
	
	f_etaLineshape = new TF1("f_etaLineshape", CrystalBall2, minFitRange, maxFitRange, 10);
	f_etaLineshape->SetParameters(fEta2->GetParameters());
	
	fEta1->Delete();
	fEta2->Delete();
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	TF1 *fEta = new TF1("fEta", DoubleGaus, 0.45, 0.65, 6);
	
	fEta->FixParameter(3, 0.0);
	fEta->FixParameter(4, 0.547);
	fEta->FixParameter(5, 0.02);
	
	double aEtaGuess     = h_etaLineshape->GetMaximum();
	double muEtaGuess    = h_etaLineshape->GetBinCenter(h_etaLineshape->GetMaximumBin());
	double sigmaEtaGuess = 0.01;
	
	// fit only region around peak first:
	
	fEta->SetRange(muEtaGuess-1.5*sigmaEtaGuess, muEtaGuess+1.5*sigmaEtaGuess);
	
	fEta->SetParameter(0,     aEtaGuess);
	fEta->SetParameter(1,    muEtaGuess);
	fEta->SetParameter(2, sigmaEtaGuess);
	
	fEta->SetParLimits(0, 0.000, 1.000);
	fEta->SetParLimits(1, 0.530, 0.580);
	fEta->SetParLimits(2, 0.005, 0.030);
	
	h_etaLineshape->Fit(fEta, "R0Q");
	
	// Now let second Gaussian float:
	
	fEta->ReleaseParameter(3);
	fEta->SetParameter(3, 0.1*fEta->GetParameter(0));
	fEta->SetParLimits(3, 0.0, 1.0);
	
	fEta->ReleaseParameter(4);
	fEta->SetParameter(4, fEta->GetParameter(1));
	fEta->SetParLimits(4, 0.530, 0.580);
	
	fEta->ReleaseParameter(5);
	fEta->SetParameter(5, 2.0*fEta->GetParameter(2));
	fEta->SetParLimits(5, 0.010, 0.050);
	
	fEta->SetRange(0.5, 0.6);
	
	h_etaLineshape->Fit(fEta, "R0Q");
	
	if(drawOption) {
		TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
		cEtaLS->SetLogy();
		h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
		h_etaLineshape->Draw();
		fEta->SetRange(0.4, 0.7);
		fEta->Draw("same");
		
		TF1 *fEtaA = new TF1("fEtaA", "gaus", 0.4, 0.7);
		for(int i=0; i<3; i++) fEtaA->SetParameter(i, fEta->GetParameter(i));
		fEtaA->SetLineColor(kBlue);
		fEtaA->SetLineStyle(2);
		
		TF1 *fEtaB = new TF1("fEtaB", "gaus", 0.4, 0.7);
		for(int i=0; i<3; i++) fEtaB->SetParameter(i, fEta->GetParameter(i+3));
		fEtaB->SetLineColor(kMagenta);
		fEtaB->SetLineStyle(2);
		
		//fEtaA->Draw("same");
		//fEtaB->Draw("same");
		
		cEtaLS->Update();
		getchar();
		delete cEtaLS;
		delete fEtaA;
		delete fEtaB;
	}
	
	f_etaLineshape = new TF1("f_etaLineshape", DoubleGaus, minFitRange, maxFitRange, 6);
	f_etaLineshape->SetParameters(fEta->GetParameters());
	
	fEta->Delete();
	*/
	return;
}

void MggFitter::FitEtaPionLineshape(int drawOption)
{
	if(fitOption_signal<7) return;
	if(h_etaPionLineshape==NULL) return;
	
	TF1 *fEtaPion1 = new TF1("fEtaPion1", CrystalBall_flip, 0.5, 0.7, 5);
	
	double     nEtaPionGuess = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(0.5), 
		h_etaPionLineshape->FindBin(0.625));
	double    muEtaPionGuess =  0.56;
	double sigmaEtaPionGuess =  0.015;
	double alphaEtaPionGuess =  1.0;
	double    nnEtaPionGuess =  2.0;
	
	double locEtaPionMax = 0.0;
	for(int ibin=h_etaPionLineshape->FindBin(0.5); ibin<=h_etaPionLineshape->FindBin(0.6); ibin++) {
		if(h_etaPionLineshape->GetBinContent(ibin) > locEtaPionMax) {
			locEtaPionMax  = h_etaPionLineshape->GetBinContent(ibin);
			muEtaPionGuess = h_etaPionLineshape->GetBinCenter(ibin);
		}
	}
	
	fEtaPion1->SetParameter(0,     nEtaPionGuess);
	fEtaPion1->SetParameter(1,    muEtaPionGuess);
	fEtaPion1->SetParameter(2, sigmaEtaPionGuess);
	fEtaPion1->SetParameter(3, alphaEtaPionGuess);
	fEtaPion1->SetParameter(4,    nnEtaPionGuess);
	
	fEtaPion1->SetParLimits(0,  0.000,  1.e6);
	fEtaPion1->SetParLimits(1,  0.540,  0.600);
	fEtaPion1->SetParLimits(2,  0.010,  0.050);
	fEtaPion1->SetParLimits(3,  0.500,  9.999);
	fEtaPion1->SetParLimits(4,  0.100,  9.999);
	
	h_etaPionLineshape->Fit(fEtaPion1, "R0QL");
	
	TF1 *fEtaPion2 = new TF1("fEtaPion2", CrystalBall2_flip, 0.5, 0.7, 10);
	for(int i=0; i<10; i++) {
		fEtaPion2->SetParameter(i, fEtaPion1->GetParameter(i%5));
	}
	
	fEtaPion2->SetParameter(5, 0.0);
	fEtaPion2->SetParameter(7, fEtaPion1->GetParameter(2)*2.0);
	
	fEtaPion2->SetParLimits(0,  0.000,  1.e6);
	fEtaPion2->SetParLimits(1,  0.540,  0.600);
	fEtaPion2->SetParLimits(2,  0.010,  0.050);
	fEtaPion2->SetParLimits(3,  0.500,  9.999);
	fEtaPion2->SetParLimits(4,  0.100,  9.999);
	
	fEtaPion2->SetParLimits(5,  0.000,  1.e6);
	fEtaPion2->SetParLimits(6,  0.540,  0.620);
	fEtaPion2->SetParLimits(7,  0.010,  0.100);
	fEtaPion2->SetParLimits(8,  0.500,  9.999);
	fEtaPion2->SetParLimits(9,  0.100,  9.999);
	
	h_etaPionLineshape->Fit(fEtaPion2, "R0QL");
	
	if(drawOption) {
		TCanvas *cEtaPionLS = new TCanvas("cEtaPionLS", "cEtaPionLS", 950, 700);
		//cEtaPionLS->SetLogy();
		h_etaPionLineshape->Draw();
		fEtaPion2->Draw("same");
		
		TF1 *fEtaPionA = new TF1("fEtaPionA", CrystalBall_flip, 0.2, 1.0, 5);
		for(int i=0; i<5; i++) fEtaPionA->SetParameter(i, fEtaPion2->GetParameter(i));
		fEtaPionA->SetLineColor(kBlue);
		fEtaPionA->SetLineStyle(2);
		
		TF1 *fEtaPionB = new TF1("fEtaPionB", CrystalBall_flip, 0.2, 1.0, 5);
		for(int i=0; i<5; i++) fEtaPionB->SetParameter(i, fEtaPion2->GetParameter(i+5));
		fEtaPionB->SetLineColor(kMagenta);
		fEtaPionB->SetLineStyle(2);
		
		fEtaPionA->Draw("same");
		fEtaPionB->Draw("same");
		
		cEtaPionLS->Update();
		getchar();
		delete cEtaPionLS;
	}
	
	f_etaPionLineshape = new TF1("f_etaPionLineshape", CrystalBall2_flip, minFitRange, maxFitRange, 10);
	f_etaPionLineshape->SetParameters(fEtaPion2->GetParameters());
	
	fEtaPion1->Delete();
	fEtaPion2->Delete();
	
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

void MggFitter::FitFDCOmegaLineshape(int drawOption)
{
	if(h_fdcOmegaLineshape==NULL) return;
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBall, 0.2, 0.7, 5);
	
	double    muOmegaGuess = h_fdcOmegaLineshape->GetBinCenter(h_fdcOmegaLineshape->GetMaximumBin());
	double     nOmegaGuess = h_fdcOmegaLineshape->Integral(h_fdcOmegaLineshape->FindBin(muOmegaGuess-0.05), 
		h_fdcOmegaLineshape->FindBin(muOmegaGuess+0.05));
	double sigmaOmegaGuess = 0.025;
	double alphaOmegaGuess = 1.0;
	double    nnOmegaGuess = 2.0;
	
	fOmega1->SetParameter(0,     nOmegaGuess);
	fOmega1->SetParameter(1,    muOmegaGuess);
	fOmega1->SetParameter(2, sigmaOmegaGuess);
	fOmega1->SetParameter(3, alphaOmegaGuess);
	fOmega1->SetParameter(4,    nnOmegaGuess);
	
	fOmega1->SetParLimits(0, 0.000, 1.e6);
	fOmega1->SetParLimits(1, muOmegaGuess-0.05, muOmegaGuess+0.05);
	fOmega1->SetParLimits(2, 0.005, 0.050);
	
	h_fdcOmegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", CrystalBall2, 0.2, 0.7, 10);
	for(int i=0; i<10; i++) {
		fOmega2->SetParameter(i, fOmega1->GetParameter(i%5));
	}
	
	fOmega2->SetParameter(5, 0.0);
	fOmega2->SetParameter(7, fOmega1->GetParameter(2)*2.0);
	
	fOmega2->SetParLimits(0, 0.000, 1.e6);
	fOmega2->SetParLimits(1, muOmegaGuess-0.05, muOmegaGuess+0.05);
	fOmega2->SetParLimits(2, 0.005, 0.050);
	
	fOmega2->SetParLimits(5, 0.000, 1.e6);
	fOmega2->SetParLimits(6, muOmegaGuess-0.10, muOmegaGuess+0.05);
	fOmega2->SetParLimits(7, 0.005, 0.200);
	
	h_fdcOmegaLineshape->Fit(fOmega2, "R0QL");
	/*
	TF1 *fOmega3 = new TF1("fOmega3", CrystalBall3, 0.2, 0.9, 10);
	for(int i=0; i<15; i++) {
		fOmega3->SetParameter(i, fOmega2->GetParameter(i%5));
	}
	
	fOmega3->SetParameter(10, 0.0);
	fOmega3->SetParameter(12, fOmega2->GetParameter(2)*2.0);
	
	fOmega2->SetParLimits(0, 0.000, 1.e6);
	fOmega1->SetParLimits(1, muOmegaGuess-0.05, muOmegaGuess+0.05);
	fOmega2->SetParLimits(2, 0.005, 0.050);
	
	fOmega2->SetParLimits(5, 0.000, 1.e6);
	fOmega1->SetParLimits(6, muOmegaGuess-0.1, muOmegaGuess+0.1);
	fOmega2->SetParLimits(7, 0.005, 0.100);
	
	h_fdcOmegaLineshape->Fit(fOmega2, "R0QL");
	*/
	if(drawOption) {
		TCanvas *cFDCOmegaLS = new TCanvas("cFDCOmegaLS", "cFDCOmegaLS", 950, 700);
		//cFDCOmegaLS->SetLogy();
		h_fdcOmegaLineshape->Draw();
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
		
		cFDCOmegaLS->Update();
		getchar();
		delete cFDCOmegaLS;
	}
	
	f_fdcOmegaLineshape = new TF1("f_fdcOmegaLineshape", CrystalBall2, minFitRange, maxFitRange, 10);
	f_fdcOmegaLineshape->SetParameters(fOmega2->GetParameters());
	
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
	int nParameters = InitializeFitFunction(&locfBkgd, "locBkgdClone");
	locfBkgd->SetParameters(f_fit->GetParameters());
	ZeroSignalPars(locfBkgd);
	
	if(useFitPars==0) {
		for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
			double locData = h_data->GetBinContent(ibin);
			double locBkgd = locfBkgd->Eval(h_data->GetBinCenter(ibin));
			yield    += locData - locBkgd;
			yieldErr += pow(h_data->GetBinError(ibin),2.0) + locBkgd;
			//printf("  %f: %f  %f\n", h_data->GetBinCenter(ibin), locData, locBkgd);
		}
		yieldErr = sqrt(yieldErr);
		/*
		// For debugging:
		printf("\n\n");
		for(int ipar=0; ipar<nParameters; ipar++) {
			printf(" %s: %f +/- %f\n", locfBkgd->GetParName(ipar), locfBkgd->GetParameter(ipar), locfBkgd->GetParError(ipar));
		}
		*/
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
				
				int locMinMggBin = h_etaLineshape->FindBin(locMinMggCut + 0.5*h_etaLineshape->GetBinWidth(1));
				int locMaxMggBin = h_etaLineshape->FindBin(locMaxMggCut - 0.5*h_etaLineshape->GetBinWidth(1));
				
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
				
				// add in the yield measured by the additional Gaussian:
				
				double locN     = f_fit->GetParameter("N_{#eta,inc}");
				double locMu    = f_fit->GetParameter("#mu_{#eta,inc}");
				double locSigma = f_fit->GetParameter("#sigma_{#eta,inc}");
				double locA     = locN * binSize / sqrt(2.0*TMath::Pi()) / locSigma;
				
				yield += (IntegrateGaussian(locA, locMu, locSigma, locMinMggCut, locMaxMggCut) / binSize);
				
				yieldErr = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				break;
			}
			case 7:
			{
				// Eta Lineshape + Eta+Pion Lineshape:
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				int locMinMggBin = h_etaLineshape->FindBin(locMinMggCut + 0.5*h_etaLineshape->GetBinWidth(1));
				int locMaxMggBin = h_etaLineshape->FindBin(locMaxMggCut - 0.5*h_etaLineshape->GetBinWidth(1));
				
				double lsMinMggCut = h_etaLineshape->GetBinCenter(locMinMggBin) - 0.5*h_etaLineshape->GetBinWidth(1);
				double lsMaxMggCut = h_etaLineshape->GetBinCenter(locMaxMggBin) + 0.5*h_etaLineshape->GetBinWidth(1);
				
				// check that lsMinMggCut and lsMaxMggCut align with desired cut range:
				if((fabs(lsMinMggCut-locMinMggCut)>1.e-6) || (fabs(lsMaxMggCut-locMaxMggCut)>1.e-6)) {
					printf("\n\nWarning: Mgg binning of eta lineshape does not overlap with data.\n");
					printf("  Lineshape cut range: %f GeV - %f GeV\n", lsMinMggCut, lsMaxMggCut);
					printf("  Data cut range: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
				}
				
				double binSizeRatio = h_etaLineshape->GetBinWidth(1) / binSize;
				
				//yield = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar) * binSizeRatio;
				yield = (f_etaLineshape->Integral(0.5,0.6)/binSize) * f_fit->GetParameter(yieldPar);
				
				// add in the yield measured by the additional Gaussian:
				
				int yieldPar2 = f_fit->GetParNumber("N_{#eta#pi}");
				
				int locMinMggBin2 = h_etaPionLineshape->FindBin(locMinMggCut + 0.5*h_etaPionLineshape->GetBinWidth(1));
				int locMaxMggBin2 = h_etaPionLineshape->FindBin(locMaxMggCut - 0.5*h_etaPionLineshape->GetBinWidth(1));
				
				double lsMinMggCut2 = h_etaPionLineshape->GetBinCenter(locMinMggBin2) - 0.5*h_etaPionLineshape->GetBinWidth(1);
				double lsMaxMggCut2 = h_etaPionLineshape->GetBinCenter(locMaxMggBin2) + 0.5*h_etaPionLineshape->GetBinWidth(1);
				
				// check that lsMinMggCut and lsMaxMggCut align with desired cut range:
				if((fabs(lsMinMggCut2-locMinMggCut)>1.e-6) || (fabs(lsMaxMggCut2-locMaxMggCut)>1.e-6)) {
					printf("\n\nWarning: Mgg binning of eta+pion lineshape does not overlap with data.\n");
					printf("  Lineshape cut range: %f GeV - %f GeV\n", lsMinMggCut2, lsMaxMggCut2);
					printf("  Data cut range: %f GeV - %f GeV\n\n", locMinMggCut, locMaxMggCut);
				}
				
				double binSizeRatio2 = h_etaPionLineshape->GetBinWidth(1) / binSize;
				
				//yield += (h_etaPionLineshape->Integral(locMinMggBin2, locMaxMggBin2) * f_fit->GetParameter(yieldPar2) * binSizeRatio2);
				//yield += ((f_etaPionLineshape->Integral(0.5,0.6)/binSize) * f_fit->GetParameter(yieldPar2));
				
				
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
			f1->SetParameter("N_{#eta,inc}", 0.0);
			break;
		case 7:
			f1->SetParameter("N_{#eta}",    0.0);
			f1->SetParameter("N_{#eta#pi}", 0.0);
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
	f1->SetParameter("N_{#omega}",  0.0);
	f1->SetParameter("N_{#eta'}",   0.0);
	f1->SetParameter("N_{#eta#pi}", 0.0);
	f1->SetParameter("N_{empty}",   0.0);
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

void MggFitter::FillEmptyPull(TH1F *h_pull)
{
	if(h_pull==NULL) {
		return;
	}
	
	// Just check that the binning of supplied histogram and h_data are consistent:
	if(h_empty->GetXaxis()->GetNbins() != h_empty->GetXaxis()->GetNbins()) {
		cout << "\nWarning: Issue with binning of empty pull histogram.\n" << endl;
	}
	
	for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_mgg = h_empty->GetXaxis()->GetBinCenter(ibin);
		double loc_unc = h_empty->GetBinError(ibin);
		if(loc_unc <= 1.0) loc_unc = 1.0;
		h_pull->SetBinContent(ibin, (h_empty->GetBinContent(ibin) - f_empty->Eval(loc_mgg))/loc_unc);
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

void MggFitter::DumpFitParameters()
{
	printf("\n\nFit parameters:\n");
	for(int ipar=0; ipar<m_nParameters; ipar++) {
		printf("  p%d (%s): %f +/- %f\n", ipar, f_fit->GetParName(ipar), f_fit->GetParameter(ipar), f_fit->GetParError(ipar));
	}
	printf("\n\n");
	return;
}
