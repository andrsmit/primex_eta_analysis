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
	/* get initial guesses for omega fit parameters */
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
		
		int frac_etapi_par = f_fit->GetParNumber("frac_{#eta#pi}");
		f_fit->ReleaseParameter(frac_etapi_par);
		f_fit->SetParLimits(frac_etapi_par, 0.0, 0.5);
		/*
		f_fit->SetParLimits(frac_etapi_par, 
			m_etaPionFraction - m_etaPionFraction, 
			m_etaPionFraction + m_etaPionFraction);
		*/
		h_data->Fit(f_fit, fitOption);
		
		//---------------------------------//
		// allow 'shift' of eta+pion background to change compared with exclusive signal:
		/*
		int dmu_etapi_par = f_fit->GetParNumber("#Delta#mu_{#eta#pi}");
		f_fit->ReleaseParameter(dmu_etapi_par);
		f_fit->SetParLimits(dmu_etapi_par, -0.01, 0.01);
		
		h_data->Fit(f_fit, fitOption);
		*/
		//---------------------------------//
		
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
	
	double     NGuess = h_data->Integral(h_data->FindBin(minOmegaFit), h_data->FindBin(maxOmegaFit));
	double    muGuess = EtaAnalyzer::m_massOmega;
	double sigmaGuess = 0.025;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	double locMax = 0.0;
	for(int ibin=h_data->FindBin(0.75); ibin<=h_data->FindBin(maxOmegaFit); ibin++) {
		if(h_data->GetBinContent(ibin) > locMax) {
			locMax  = h_data->GetBinContent(ibin);
			muGuess = h_data->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_omega) {
		case 1:
		{
			int     NPar = f_fit->GetParNumber("N_{#omega}");
			int    muPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int     nPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->SetParameter(    NPar,     NGuess);
			f_fit->SetParameter(   muPar,    muGuess);
			f_fit->SetParameter(sigmaPar, sigmaGuess);
			f_fit->SetParameter(alphaPar, alphaGuess);
			f_fit->SetParameter(    nPar,    nGuess);
			
			f_fit->SetParLimits(    NPar, 0.000,  1.0e6);
			f_fit->SetParLimits(   muPar, 0.750,  0.800);
			f_fit->SetParLimits(sigmaPar, 0.015,  0.050);
			f_fit->SetParLimits(alphaPar, 0.200,  9.999);
			f_fit->SetParLimits(    nPar, 1.100, 49.999);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(   muPar, m_OMEGA_MU);
				f_fit->FixParameter(sigmaPar, m_OMEGA_SIGMA);
				f_fit->FixParameter(alphaPar, m_OMEGA_ALPHA);
				f_fit->FixParameter(    nPar, m_OMEGA_N);
			}
			*/
			break;
		}
		case 2:
		{
			// Using functional parameterization of MC lineshape:
			
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, muGuess - h_omegaLineshape->GetBinCenter(h_omegaLineshape->GetMaximumBin()));
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e6);
			f_fit->SetParLimits(dmuPar, -0.03, 0.03);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 3:
		{
			// Lineshape fit:
			
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, muGuess - h_omegaLineshape->GetBinCenter(h_omegaLineshape->GetMaximumBin()));
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuPar, -0.03, 0.03);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
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
			int     NPar = f_fit->GetParNumber("N_{#omega}");
			int    muPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int     nPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->FixParameter(    NPar, f_fit->GetParameter(    NPar));
			f_fit->FixParameter(   muPar, f_fit->GetParameter(   muPar));
			f_fit->FixParameter(sigmaPar, f_fit->GetParameter(sigmaPar));
			f_fit->FixParameter(alphaPar, f_fit->GetParameter(alphaPar));
			f_fit->FixParameter(    nPar, f_fit->GetParameter(    nPar));
			break;
		}
		case 2:
		{
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->FixParameter(  NPar, f_fit->GetParameter(  NPar));
			f_fit->FixParameter(dmuPar, f_fit->GetParameter(dmuPar));
			break;
		}
		case 3:
		{
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->FixParameter(  NPar, f_fit->GetParameter(  NPar));
			f_fit->FixParameter(dmuPar, f_fit->GetParameter(dmuPar));
			break;
		}
	}
	return;
}

void MggFitter::ReleaseOmegaParameters()
{
	switch(fitOption_omega) {
		case 1:
		{
			int     NPar = f_fit->GetParNumber("N_{#omega}");
			int    muPar = f_fit->GetParNumber("#mu_{#omega}");
			int sigmaPar = f_fit->GetParNumber("#sigma_{#omega}");
			int alphaPar = f_fit->GetParNumber("#alpha_{#omega}");
			int     nPar = f_fit->GetParNumber("n_{#omega}");
			
			f_fit->ReleaseParameter(NPar);
			f_fit->SetParLimits(NPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(muPar);
			f_fit->SetParLimits(muPar, 0.750, 0.800);
			
			f_fit->ReleaseParameter(sigmaPar);
			f_fit->SetParLimits(sigmaPar, 0.015, 0.050);
			
			f_fit->ReleaseParameter(alphaPar);
			f_fit->SetParLimits(alphaPar, 0.200, 9.999);
			
			f_fit->ReleaseParameter(nPar);
			f_fit->SetParLimits(nPar, 1.100, 49.999);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(   muPar, m_OMEGA_MU);
				f_fit->FixParameter(sigmaPar, m_OMEGA_SIGMA);
				f_fit->FixParameter(alphaPar, m_OMEGA_ALPHA);
				f_fit->FixParameter(    nPar, m_OMEGA_N);
			}
			*/
		}
		case 2:
		{
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->ReleaseParameter(NPar);
			f_fit->SetParLimits(NPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(dmuPar);
			f_fit->SetParLimits(dmuPar, -0.01, 0.01);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 3:
		{
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->ReleaseParameter(NPar);
			f_fit->SetParLimits(NPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(dmuPar);
			f_fit->SetParLimits(dmuPar, -0.01, 0.01);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
	}
	return;
}

//==============================================================//
// Background:

void MggFitter::GuessBkgdParameters()
{
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
			
			f_fit->SetParameter(p0Par, p0Guess);
			f_fit->SetParameter(p1Par, p1Guess);
			f_fit->SetParameter(p2Par, p2Guess);
			f_fit->SetParameter(p3Par, p3Guess);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e5);
			f_fit->SetParLimits(p1Par,  0.00, 1.00);
			f_fit->SetParLimits(p2Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
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

void MggFitter::FixBkgdParameters()
{
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
			for(int ipar=0; ipar<4; ipar++) {
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

void MggFitter::ReleaseBkgdParameters()
{
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
			
			f_fit->ReleaseParameter(p0Par);
			f_fit->ReleaseParameter(p1Par);
			f_fit->ReleaseParameter(p2Par);
			f_fit->ReleaseParameter(p3Par);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e5);
			f_fit->SetParLimits(p1Par,  0.00, 1.00);
			f_fit->SetParLimits(p2Par, -1.e3, 1.e3);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
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

void MggFitter::GuessEtaParameters()
{
	// guess number of eta' by integrating histogram and subtracting background:
	
	double minEtaFit = 0.50;
	double maxEtaFit = 0.60;
	
	double     NGuess = h_data->Integral(h_data->FindBin(minEtaFit), h_data->FindBin(maxEtaFit));
	double    muGuess = EtaAnalyzer::m_massEta;
	double sigmaGuess = 0.015;
	
	double locMax = 0.0;
	for(int ibin=h_data->FindBin(minEtaFit); ibin<=h_data->FindBin(maxEtaFit); ibin++) {
		if(h_data->GetBinContent(ibin) > locMax) {
			locMax  = h_data->GetBinContent(ibin);
			muGuess = h_data->GetBinCenter(ibin);
		}
	}
	
	switch(fitOption_signal) {
		case 1:
		{
			// single Gaussian:
			
			int     NPar = f_fit->GetParNumber("N_{#eta}");
			int    muPar = f_fit->GetParNumber("#mu_{#eta}");
			int sigmaPar = f_fit->GetParNumber("#sigma_{#eta}");
			
			f_fit->SetParameter(    NPar,     NGuess);
			f_fit->SetParameter(   muPar,    muGuess);
			f_fit->SetParameter(sigmaPar, sigmaGuess);
			
			f_fit->SetParLimits(    NPar, 0.00, 1.e5);
			f_fit->SetParLimits(   muPar, 0.54, 0.62);
			f_fit->SetParLimits(sigmaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			int     N1Par = f_fit->GetParNumber("N_{#eta,1}");
			int     N2Par = f_fit->GetParNumber("N_{#eta,2}");
			int    mu1Par = f_fit->GetParNumber("#mu_{#eta,1}");
			int    dmuPar = f_fit->GetParNumber("#mu_{#eta,2}-#mu_{#eta,1}");
			int sigma1Par = f_fit->GetParNumber("#sigma_{#eta,1}");
			int sigma2Par = f_fit->GetParNumber("#sigma_{#eta,2}");
			
			f_fit->SetParameter(    N1Par, 0.9*NGuess);
			f_fit->SetParameter(    N2Par, 0.1*NGuess);
			f_fit->SetParameter(   mu1Par, muGuess);
			f_fit->SetParameter(   dmuPar, 0.0);
			f_fit->SetParameter(sigma1Par, sigmaGuess);
			f_fit->SetParameter(sigma2Par, sigmaGuess*2.0);
			
			f_fit->SetParLimits(    N1Par,  0.00, 1.e5);
			f_fit->SetParLimits(    N2Par,  0.00, 0.5*NGuess);
			f_fit->SetParLimits(   mu1Par,  0.54, 0.62);
			f_fit->SetParLimits(   dmuPar,  0.00, 0.03);
			f_fit->SetParLimits(sigma1Par,  0.01, 0.03);
			f_fit->SetParLimits(sigma2Par,  0.01, 0.05);
			
			if(angle<0.5) {
				f_fit->FixParameter(    N2Par, 0.0);
				f_fit->FixParameter(   dmuPar, 0.0);
				f_fit->FixParameter(sigma2Par, 2.0*sigmaGuess);
			}
			break;
		}
		case 3:
		{
			// Crsytal ball:
			
			int     NPar = f_fit->GetParNumber("N_{#eta}");
			int    muPar = f_fit->GetParNumber("#mu_{#eta}");
			int sigmaPar = f_fit->GetParNumber("#sigma_{#eta}");
			int alphaPar = f_fit->GetParNumber("#alpha_{#eta}");
			int     nPar = f_fit->GetParNumber("n_{#eta}");
			
			f_fit->SetParameter(    NPar,     NGuess);
			f_fit->SetParameter(   muPar,    muGuess);
			f_fit->SetParameter(sigmaPar, sigmaGuess);
			f_fit->SetParameter(alphaPar,        1.0);
			f_fit->SetParameter(   nPar,        2.0);
			
			f_fit->SetParLimits(    NPar, 0.00,  1.e5);
			f_fit->SetParLimits(   muPar, 0.54,  0.62);
			f_fit->SetParLimits(sigmaPar, 0.01,  0.03);
			f_fit->SetParLimits(alphaPar, 0.20,  9.99);
			f_fit->SetParLimits(    nPar, 1.10, 99.99);
			break;
		}
		case 4:
		{
			
			break;
		}
		case 5:
		{
			// Lineshape fit:
			
			int   NPar = f_fit->GetParNumber("N_{#eta}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, 0.005);
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuPar, -0.01, 0.02);
			break;
		}
		case 6:
		{
			// Start with Lineshape fit:
			
			int   NPar = f_fit->GetParNumber("N_{#eta}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, 0.005);
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuPar, -0.01, 0.02);
			break;
		}
		case 7:
		{
			// Start with Lineshape fit:
			
			int   NPar = f_fit->GetParNumber("N_{#eta}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, 0.005);
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuPar, -0.01, 0.02);
			break;
		}
	}
	return;
}

//==============================================================//
// Empty:

void MggFitter::GuessEmptyParameters()
{
	if(fitOption_empty!=1) return;
	
	int NPar = f_fit->GetParNumber("N_{empty}");
	f_fit->SetParameter(NPar, m_emptyRatio);
	f_fit->SetParLimits(NPar, m_emptyRatio-m_emptyRatioErr, m_emptyRatio+m_emptyRatioErr);
	return;
}

void MggFitter::FixEmptyParameters()
{
	if(fitOption_empty!=1) return;
	
	int NPar = f_fit->GetParNumber("N_{empty}");
	f_fit->FixParameter(NPar, f_fit->GetParameter(NPar));
	return;
}

//==============================================================//
// Lineshape Fits:

void MggFitter::FitEtaLineshape(int drawOption)
{
	if(h_etaLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	CheckBinSize(h_etaLineshape, "Eta Lineshape");
	
	TF1 *locfEta = new TF1("locfEta", DoubleCrystalBallPDF_oneflip, 0.45, 0.65, 10);
	locfEta->SetParameters(
		0.538, //    mu1
		0.007, // sigma1
		0.925, // alpha1
		6.122, //     n1
		0.000, //    mu2-mu1
		0.010, // sigma2
		1.455, // alpha2
		9.266, //     n2
		0.300  // fraction
	);
	locfEta->FixParameter(9, h_etaLineshape->GetXaxis()->GetBinWidth(1));
	
	locfEta->SetParLimits(0, 0.530,  0.560);
	locfEta->SetParLimits(1, 0.004,  0.050);
	locfEta->SetParLimits(2, 0.200,  9.999);
	locfEta->SetParLimits(3, 1.100, 49.999);
	
	locfEta->SetParLimits(4,-0.050,  0.050);
	locfEta->SetParLimits(5, 0.004,  0.050);
	locfEta->SetParLimits(6, 0.200,  9.999);
	locfEta->SetParLimits(7, 1.100, 49.999);
	
	locfEta->SetParLimits(8, 0.0, 1.0);
	
	// Force both crystal ball shapes to have the same mean:
	locfEta->FixParameter(4, 0.0);
	
	h_etaLineshape->Fit(locfEta, "R0QL");
	
	if(drawOption) {
		TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
		h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
		h_etaLineshape->Draw();
		locfEta->Draw("same");
		h_etaLineshape->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
		
		/*
		TF1 *fEtaA = new TF1("fEtaA", CrystalBall, 0.4, 0.7, 5);
		for(int i=0; i<5; i++) fEtaA->SetParameter(i, locfEta->GetParameter(i));
		fEtaA->SetLineColor(kBlue);
		fEtaA->SetLineStyle(2);
		
		TF1 *fEtaB = new TF1("fEtaB", CrystalBall_flip, 0.4, 0.7, 5);
		for(int i=0; i<5; i++) fEtaB->SetParameter(i, locfEta->GetParameter(i+5));
		fEtaB->SetLineColor(kMagenta);
		fEtaB->SetLineStyle(2);
		
		printf("Eta lineshape fit parameters:\n");
		for(int i=0; i<10; i++) {
			printf(" p%d = %f\n", i, locfEta->GetParameter(i));
		}
		fEtaA->Draw("same");
		fEtaB->Draw("same");
		*/
		
		cEtaLS->Update();
		getchar();
		cEtaLS->SetLogy();
		cEtaLS->Update();
		getchar();
		delete cEtaLS;
	}
	
	f_etaLineshape = new TF1("f_etaLineshape", DoubleCrystalBallPDF_oneflip, minFitRange, maxFitRange, 10);
	f_etaLineshape->SetParameters(locfEta->GetParameters());
	f_etaLineshape->FixParameter(9, 1.0);
	
	locfEta->Delete();
	
	return;
}

void MggFitter::FitEtaPionLineshape(int drawOption)
{
	if(fitOption_signal<7) return;
	if(h_etaPionLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	CheckBinSize(h_etaPionLineshape, "Eta+Pion Lineshape");
	
	TF1 *fEtaPion1 = new TF1("fEtaPion1", CrystalBallPDF_flip, 0.5, 0.7, 5);
	
	double    muGuess = h_etaPionLineshape->GetBinCenter(h_etaPionLineshape->GetMaximumBin());
	double sigmaGuess = 0.015;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fEtaPion1->SetParameter(0,    muGuess);
	fEtaPion1->SetParameter(1, sigmaGuess);
	fEtaPion1->SetParameter(2, alphaGuess);
	fEtaPion1->SetParameter(3,     nGuess);
	
	fEtaPion1->SetParLimits(0, 0.540,  0.610);
	fEtaPion1->SetParLimits(1, 0.010,  0.050);
	fEtaPion1->SetParLimits(2, 0.200,  9.999);
	fEtaPion1->SetParLimits(3, 1.100, 49.999);
	
	fEtaPion1->FixParameter(4, h_etaPionLineshape->GetBinWidth(1));
	
	h_etaPionLineshape->Fit(fEtaPion1, "R0QL");
	
	TF1 *fEtaPion2 = new TF1("fEtaPion2", DoubleCrystalBallPDF_flip, 0.5, 0.7, 10);
	fEtaPion2->SetParameter(0, fEtaPion1->GetParameter(0));
	fEtaPion2->SetParameter(1, fEtaPion1->GetParameter(1));
	fEtaPion2->SetParameter(2, fEtaPion1->GetParameter(2));
	fEtaPion2->SetParameter(3, fEtaPion1->GetParameter(3));
	fEtaPion2->SetParameter(4, 0.0);
	fEtaPion2->SetParameter(5, fEtaPion1->GetParameter(1)*2.0);
	fEtaPion2->SetParameter(6, fEtaPion1->GetParameter(2));
	fEtaPion2->SetParameter(7, fEtaPion1->GetParameter(3));
	fEtaPion2->SetParameter(8, 0.0);
	
	fEtaPion2->SetParLimits(0,  0.540,  0.610);
	fEtaPion2->SetParLimits(1,  0.010,  0.050);
	fEtaPion2->SetParLimits(2,  0.200,  9.999);
	fEtaPion2->SetParLimits(3,  1.100, 49.999);
	fEtaPion2->SetParLimits(4, -0.050,  0.050);
	fEtaPion2->SetParLimits(5,  0.010,  0.050);
	fEtaPion2->SetParLimits(6,  0.200,  9.999);
	fEtaPion2->SetParLimits(7,  1.100, 49.999);
	fEtaPion2->SetParLimits(8,  0.000,  1.000);
	
	fEtaPion2->FixParameter(9, h_etaPionLineshape->GetBinWidth(1));
	
	h_etaPionLineshape->Fit(fEtaPion2, "R0QL");
	
	if(drawOption) {
		TCanvas *cEtaPionLS = new TCanvas("cEtaPionLS", "cEtaPionLS", 950, 700);
		//cEtaPionLS->SetLogy();
		h_etaPionLineshape->Draw();
		fEtaPion2->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = fEtaPion2->Integral(0.5,0.6) / h_etaPionLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of eta+pion lineshape within mgg cut: %f\n", fracAccepted);
		
		cEtaPionLS->Update();
		getchar();
		delete cEtaPionLS;
	}
	
	f_etaPionLineshape = new TF1("f_etaPionLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
	f_etaPionLineshape->SetParameters(fEtaPion2->GetParameters());
	f_etaPionLineshape->FixParameter(9, 1.0);
	
	fEtaPion1->Delete();
	fEtaPion2->Delete();
	
	return;
}

void MggFitter::FitOmegaLineshape(int drawOption)
{
	if(fitOption_omega!=2) return;
	if(h_omegaLineshape==nullptr) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_omegaLineshape, "Omega Lineshape");
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBallPDF, 0.2, 0.9, 5);
	
	double    muGuess = 0.78;
	double sigmaGuess = 0.025;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	double locMax = 0.0;
	for(int ibin=h_omegaLineshape->FindBin(0.75); ibin<=h_omegaLineshape->FindBin(0.85); ibin++) {
		if(h_omegaLineshape->GetBinContent(ibin) > locMax) {
			locMax  = h_omegaLineshape->GetBinContent(ibin);
			muGuess = h_omegaLineshape->GetBinCenter(ibin);
		}
	}
	
	fOmega1->SetParameter(0,    muGuess);
	fOmega1->SetParameter(1, sigmaGuess);
	fOmega1->SetParameter(2, alphaGuess);
	fOmega1->SetParameter(3,     nGuess);
	
	fOmega1->SetParLimits(0, 0.750,  0.800);
	fOmega1->SetParLimits(1, 0.015,  0.050);
	fOmega1->SetParLimits(2, 0.200,  9.999);
	fOmega1->SetParLimits(3, 1.100, 49.999);
	
	fOmega1->FixParameter(4, h_omegaLineshape->GetBinWidth(1));
	
	h_omegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", DoubleCrystalBallPDF, 0.2, 0.9, 10);
	fOmega2->SetParameter(0, fOmega1->GetParameter(0));
	fOmega2->SetParameter(1, fOmega1->GetParameter(1));
	fOmega2->SetParameter(2, fOmega1->GetParameter(2));
	fOmega2->SetParameter(3, fOmega1->GetParameter(3));
	fOmega2->SetParameter(4, 0.0);
	fOmega2->SetParameter(5, fOmega1->GetParameter(1)*2.0);
	fOmega2->SetParameter(6, fOmega1->GetParameter(2));
	fOmega2->SetParameter(7, fOmega1->GetParameter(3));
	fOmega2->SetParameter(8, 0.0);
	
	fOmega2->SetParLimits(0,  0.750,  0.800);
	fOmega2->SetParLimits(1,  0.015,  0.050);
	fOmega2->SetParLimits(2,  0.200,  9.999);
	fOmega2->SetParLimits(3,  1.100, 49.999);
	fOmega2->SetParLimits(4, -0.050,  0.050);
	fOmega2->SetParLimits(5,  0.010,  0.100);
	fOmega2->SetParLimits(6,  0.200,  9.999);
	fOmega2->SetParLimits(7,  1.100, 49.999);
	fOmega2->SetParLimits(8,  0.000,  1.000);
	
	fOmega2->FixParameter(9, h_omegaLineshape->GetBinWidth(1));
	
	h_omegaLineshape->Fit(fOmega2, "R0QL");
	
	if(drawOption) {
		TCanvas *cOmegaLS = new TCanvas("cOmegaLS", "cOmegaLS", 950, 700);
		cOmegaLS->SetLogy();
		h_omegaLineshape->Draw();
		fOmega2->Draw("same");
		
		cOmegaLS->Update();
		getchar();
		delete cOmegaLS;
	}
	
	f_omegaLineshape = new TF1("f_omegaLineshape", DoubleCrystalBallPDF, minFitRange, maxFitRange, 10);
	f_omegaLineshape->SetParameters(fOmega2->GetParameters());
	f_omegaLineshape->FixParameter(9, 1.0);
	
	fOmega1->Delete();
	fOmega2->Delete();
	
	return;
}

void MggFitter::FitFDCOmegaLineshape(int drawOption)
{
	if(h_fdcOmegaLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_fdcOmegaLineshape, "FDC Omega Lineshape");
	
	TF1 *fOmega1 = new TF1("fOmega1", CrystalBallPDF, 0.2, 0.7, 5);
	
	double    muGuess = h_fdcOmegaLineshape->GetBinCenter(h_fdcOmegaLineshape->GetMaximumBin());
	double sigmaGuess = 0.025;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fOmega1->SetParameter(0,    muGuess);
	fOmega1->SetParameter(1, sigmaGuess);
	fOmega1->SetParameter(2, alphaGuess);
	fOmega1->SetParameter(3,     nGuess);
	
	fOmega1->SetParLimits(0, muGuess-0.05, muGuess+0.05);
	fOmega1->SetParLimits(1, 0.005,  0.050);
	fOmega1->SetParLimits(2, 0.200,  9.999);
	fOmega1->SetParLimits(3, 1.100, 49.999);
	
	fOmega1->FixParameter(4, h_fdcOmegaLineshape->GetBinWidth(1));
	
	h_fdcOmegaLineshape->Fit(fOmega1, "R0QL");
	
	TF1 *fOmega2 = new TF1("fOmega2", DoubleCrystalBallPDF, 0.2, 0.7, 10);
	fOmega2->SetParameter(0, fOmega1->GetParameter(0));
	fOmega2->SetParameter(1, fOmega1->GetParameter(1));
	fOmega2->SetParameter(2, fOmega1->GetParameter(2));
	fOmega2->SetParameter(3, fOmega1->GetParameter(3));
	fOmega2->SetParameter(4, 0.0);
	fOmega2->SetParameter(5, fOmega1->GetParameter(1)*2.0);
	fOmega2->SetParameter(6, fOmega1->GetParameter(2));
	fOmega2->SetParameter(7, fOmega1->GetParameter(3));
	fOmega2->SetParameter(8, 0.0);
	
	fOmega2->SetParLimits(0,  muGuess-0.05, muGuess+0.05);
	fOmega2->SetParLimits(1,  0.005,  0.050);
	fOmega2->SetParLimits(2,  0.200,  9.999);
	fOmega2->SetParLimits(3,  1.100, 49.999);
	fOmega2->SetParLimits(4, -0.050,  0.050);
	fOmega2->SetParLimits(5,  0.005,  0.150);
	fOmega2->SetParLimits(6,  0.200,  9.999);
	fOmega2->SetParLimits(7,  1.100, 49.999);
	fOmega2->SetParLimits(8,  0.000,  1.000);
	
	fOmega2->FixParameter(9, h_fdcOmegaLineshape->GetBinWidth(1));
	
	h_fdcOmegaLineshape->Fit(fOmega2, "R0QL");
	
	if(drawOption) {
		TCanvas *cFDCOmegaLS = new TCanvas("cFDCOmegaLS", "cFDCOmegaLS", 950, 700);
		//cFDCOmegaLS->SetLogy();
		h_fdcOmegaLineshape->Draw();
		fOmega2->Draw("same");
		
		cFDCOmegaLS->Update();
		getchar();
		delete cFDCOmegaLS;
	}
	
	f_fdcOmegaLineshape = new TF1("f_fdcOmegaLineshape", DoubleCrystalBallPDF, minFitRange, maxFitRange, 10);
	f_fdcOmegaLineshape->SetParameters(fOmega2->GetParameters());
	f_fdcOmegaLineshape->FixParameter(9, 1.0);
	
	fOmega1->Delete();
	fOmega2->Delete();
	
	return;
}

void MggFitter::GetYield(double &yield, double &yieldErr, int useFitPars, int subtractEtaPion) {
	
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
	ZeroSignalPars(locfBkgd, subtractEtaPion);
	
	if(useFitPars==0) {
		for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
			double locData = h_data->GetBinContent(ibin);
			double locBkgd = locfBkgd->Eval(h_data->GetBinCenter(ibin));
			yield    += locData - locBkgd;
			yieldErr += pow(h_data->GetBinError(ibin),2.0) + locBkgd;
		}
		yieldErr = sqrt(yieldErr);
		
		// For debugging:
		/*
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
				int yieldPar    = f_fit->GetParNumber("N_{#eta}");
				
				double locMu    = f_fit->GetParameter("#mu_{#eta}");
				double locSigma = f_fit->GetParameter("#sigma_{#eta}");
				double locInt   = IntegrateGaussian(locMu, locSigma, locMinMggCut, locMaxMggCut);
				
				double locN     = f_fit->GetParameter(yieldPar);
				double locNErr  = f_fit->GetParError(yieldPar);
				
				yield    = locN    * locInt;
				yieldErr = locNErr * locInt;
				break;
			}
			case 2:
			{
				// Double Gaussian:
				int yieldPar     = f_fit->GetParNumber("N_{#eta}");
				int fractionPar  = f_fit->GetParNumber("fraction_{#eta}");
				
				double locYield       = f_fit->GetParameter(yieldPar);
				double locYieldErr    = f_fit->GetParError(yieldPar);
				double locFraction    = f_fit->GetParameter(fractionPar);
				double locFractionErr = f_fit->GetParError(fractionPar);
				
				double locMu1    = f_fit->GetParameter("#mu_{#eta,1}");
				double locSigma1 = f_fit->GetParameter("#sigma_{#eta,1}");
				double locInt1   = IntegrateGaussian(locMu1, locSigma1, locMinMggCut, locMaxMggCut);
				
				double locMu2    = f_fit->GetParameter("#mu_{#eta,2}-#mu_{#eta,1}") + locMu1;
				double locSigma2 = f_fit->GetParameter("#sigma_{#eta,2}");
				double locInt2   = IntegrateGaussian(locMu2, locSigma2, locMinMggCut, locMaxMggCut);
				
				double locInt    = (1.0 - locFraction)*locInt1 + locFraction*locInt2;
				double locIntErr = locFractionErr * fabs(locInt2 - locInt1);
				
				yield    = locYield * locInt;
				yieldErr = sqrt(pow(locYieldErr*locInt,2.0) + pow(locYield*locIntErr,2.0));
				break;
			}
			case 3:
			{
				// Crystal Ball:
				TF1 *locPDF = new TF1("locCrystalBallPDF", CrystalBallPDF_flip, 0.5, 0.6, 5);
				locPDF->SetParameters(
					f_fit->GetParameter("#mu_{#eta}"),
					f_fit->GetParameter("#sigma_{#eta}"),
					f_fit->GetParameter("#alpha_{#eta}"),
					f_fit->GetParameter("n_{#eta}"),
					1.0
				);
				double locInt = locPDF->Integral(locMinMggCut, locMaxMggCut);
				
				int yieldPar = f_fit->GetParNumber("N_{#eta}");
				
				yield    = locInt * f_fit->GetParameter(yieldPar);
				yieldErr = locInt * f_fit->GetParError(yieldPar);
				delete locPDF;
				break;
			}
			case 4:
			{
				// Crystal Ball + Gaussian (not yet fully implemented):
				yield    = 0.0;
				yieldErr = 0.0;
				break;
			}
			case 5:
			{
				// Lineshape fit (with histogram):
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
				
				yield    = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParameter(yieldPar);
				yieldErr = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin) * f_fit->GetParError(yieldPar);
				break;
			}
			case 6:
			{
				// Lineshape (with histogram) + Extra Gaussian:
				int yieldPar    = f_fit->GetParNumber("N_{#eta}");
				int fractionPar = f_fit->GetParNumber("fraction_{#eta#pi}");
				
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
				
				double locIntExc      = h_etaLineshape->Integral(locMinMggBin, locMaxMggBin);
				double locYieldExc    = locIntExc * f_fit->GetParameter(yieldPar);
				double locYieldExcErr = locIntExc * f_fit->GetParError(yieldPar);
				
				double locIntEtapi = IntegrateGaussian(
					f_fit->GetParameter("#mu_{#eta#pi}"), f_fit->GetParameter("#sigma_{#eta#pi}"), 
					locMinMggCut, locMaxMggCut);
				double locYieldEtapi    = locIntEtapi * f_fit->GetParameter(yieldPar) * f_fit->GetParameter(fractionPar);
				double locYieldEtapiErr = locIntEtapi * sqrt(
					pow(f_fit->GetParError( yieldPar)*f_fit->GetParameter(fractionPar),2.0) +
					pow(f_fit->GetParameter(yieldPar)*f_fit->GetParError( fractionPar), 2.0)
				);
				
				if(subtractEtaPion) {
					yield    = locYieldExc;
					yieldErr = locYieldExcErr;
				} else {
					yield    = locYieldExc + locYieldEtapi;
					yieldErr = sqrt(pow(locYieldExcErr,2.0) + pow(locYieldEtapiErr,2.0));
				}
				break;
			}
			case 7:
			{
				// Eta Lineshape (PDF) + Eta+Pion Lineshape (PDF):
				
				int yieldPar    = f_fit->GetParNumber("N_{#eta}");
				int fractionPar = f_fit->GetParNumber("frac_{#eta#pi}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double frac_etapi    = f_fit->GetParameter(fractionPar);
				double frac_etapiErr = f_fit->GetParError(fractionPar);
				
				double N_etapi    = N_eta * frac_etapi;
				double N_etapiErr = sqrt(pow(N_etaErr*frac_etapi,2.0) + pow(N_eta*frac_etapiErr,2.0));
				
				/*
				'N_eta' and 'N_etapi' above represent the yield of exclusive eta's and eta+pions integrated over all mgg.
				But to be consistent with our efficiency correction that will be applied later, we need to correct
				this for the small fraction of events that fall outside our mgg cut range.
				
				Recall that at this point f_etaLineshape and f_etaPionLineshape are normalized to have unit-integrals.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut, maxMggCut);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrection_etapi = f_etaPionLineshape->Integral(minMggCut, maxMggCut);
				double locYieldEtaPi    = locCorrection_etapi * N_etapi;
				double locYieldEtaPiErr = locCorrection_etapi * N_etapiErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPi;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPiErr,2.0));
				
				if(subtractEtaPion) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
		}
	}
	
	delete locfBkgd;
	return;
}

void MggFitter::GetEmptyYield(double &yield, double &yieldErr) {
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeEmptyFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_empty->GetParameters());
	locfEta->SetParameter(locfEta->GetNpar()-1, emptyBinSize);
	
	// zero the parameters not associated with peaking structure in eta mass region:
	
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		//printf(" p%d: %s = %f\n", ipar, locParName.Data(), locfEta->GetParameter(ipar));
		
		if(locParName.Contains("eta") || locParName.Contains("sigma")) continue;
		
		if(locParName.Contains("fdc,2")) {
			//
			// Comment out the following line to exclude the 
			// portion of the empty fit function which is meant to 
			// describe the omegas coming from the second FDC package
			// from the yield of etas estimated from the empty target
			// runs.
			// 
			//continue;
		}
		
		locfEta->SetParameter(ipar,0.0);
		locfEta->SetParError(ipar,0.0);
	}
	/*
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		printf(" p%d: %s = %f\n", ipar, locParName.Data(), locfEta->GetParameter(ipar));
	}
	*/
	
	//-----------------------------------------------//
	
	int minMggBin = h_empty->FindBin(minMggCut);
	int maxMggBin = h_empty->FindBin(maxMggCut)-1;
	
	// integrate empty fit function and divide by bin size (make sure h_empty has same bin size as h_data):
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locEtas = locfEta->Eval(h_empty->GetBinCenter(ibin));
		yield += locEtas;
	}
	
	// correct this yield by normalization parameter found by fit to full data:
	
	double A_empty = f_fit->GetParameter("N_{empty}");
	
	yield   *= (A_empty * m_emptyRatio);
	yieldErr = sqrt(yield);
	
	delete locfEta;
	return;
}

void MggFitter::ZeroSignalPars(TF1 *f1, int subtractEtaPion)
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
			if(subtractEtaPion) {
				// If we are subtracting the eta+pion background, we
				// don't consider it as 'signal' parameter. So when
				// we zero the signal parameters, we only zero the portion
				// due to exclusive etas.
				f1->SetParameter("frac_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}", 0.0);
			}
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
	f1->SetParameter("frac_{#eta#pi}", 0.0);
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

double MggFitter::IntegrateGaussian(double mu, double sigma, double x1, double x2)
{
	// Integrates normalized Gaussian PDF of mean, mu, and width, sigma, between x1 and x2
	// Normalized means that if x1=-inf and x2=+inf, Integral = 1.
	
	double erf1 = TMath::Erf((x1-mu)/(sqrt(2.0)*sigma));
	double erf2 = TMath::Erf((x2-mu)/(sqrt(2.0)*sigma));
	double integral = 0.5 * (erf2 - erf1);
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

void MggFitter::GetEmptyEtaFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	// get yield of etas from fit to empty target:
	
	double locYieldEmpty, locYieldEmptyErr;
	GetEmptyYield(locYieldEmpty, locYieldEmptyErr);
	
	// get yield of etas from fit to full target:
	
	double locYield, locYieldErr;
	GetYield(locYield, locYieldErr);
	
	fraction    = locYieldEmpty / locYield;
	fractionErr = sqrt(pow(locYieldEmptyErr/locYield,2.0) 
		+ pow(locYieldErr*locYieldEmpty/(locYield*locYield),2.0));
	return;
}

void MggFitter::GetEtaPionFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	// redundant check:
	if(fitOption_signal!=7) return;
	
	int fractionPar = f_fit->GetParNumber("frac_{#eta#pi}");
	
	fraction    = f_fit->GetParameter(fractionPar);
	fractionErr = f_fit->GetParError(fractionPar);
	return;
}

void MggFitter::CheckBinSize(TH1F *h1, TString histTitle)
{
	double locBinWidth = h1->GetXaxis()->GetBinWidth(1);
	if(fabs(locBinWidth-binSize)>1.e-6) {
		printf("\n\nWARNING IN FIT OF %s:\n", histTitle.Data());
		printf("  Bin width does not match MggFitter object (%.3f)\n\n", locBinWidth);
	}
	
	return;
}
