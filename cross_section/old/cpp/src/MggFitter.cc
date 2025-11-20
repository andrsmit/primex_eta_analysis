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
	fitOption = "R0Q";
	
	InitializeFitFunction(&f_fit);
	m_nParameters = InitializeFitParameters();
	
	int hadronicBkgdFracPar = f_fit->GetParNumber("frac_{bkgd}");
	int  etaPionBkgdFracPar = f_fit->GetParNumber("frac_{#eta#pi}");
	
	if(useRawMass) {
		// If fitting the raw invariant mass, don't try to separate hadronic background:
		if((fitOption_signal==7) || (fitOption_signal==8)) {
			f_fit->FixParameter(hadronicBkgdFracPar, 0.0);
		}
		else if(fitOption_signal>8) {
			f_fit->FixParameter( etaPionBkgdFracPar, 0.0);
			f_fit->FixParameter(hadronicBkgdFracPar, 0.0);
		}
	}
	
	//-----------------------------------------------------------------//
	// release omega parameters and fit the region around the peak:
	
	// TEMPORARILY COMMENT OUT THE FIT TO THE OMEGA PEAK:
	
	if(fitOption_omega>0) {
		double minOmegaFit = 0.60;
		double maxOmegaFit = 0.90;
		
		h_data->GetXaxis()->SetRangeUser(minOmegaFit, maxOmegaFit);
		
		GuessOmegaParameters();
		f_fit->SetRange(minOmegaFit, maxOmegaFit);
		h_data->Fit(f_fit, fitOption);
		FixOmegaParameters();
	}
	
	h_data->GetXaxis()->SetRangeUser(minFitRange, maxFitRange);
	
	excludeRegions.clear();
	f_fit->SetRange(0.50, 0.60);
	GuessEtaParameters();
	
	int deltaMuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
	f_fit->SetParameter(deltaMuPar, 0.002);
	f_fit->SetParLimits(deltaMuPar,-0.01,0.02);
	
	// shift fitted from veto option 7 with hadronic bkgd fractions fixed from bggen mc (and using coherent lineshape for all angles):
	//double deltaMuGuess = 0.00176841 + (-0.000136172*angle) + (0.000218961*pow(angle,2.0));
	
	double deltaMuGuess = 0.553356
		+ ((-0.00492198)*pow(angle,1.0))
		+ (( 0.00391119)*pow(angle,2.0))
		+ ((-0.00103892)*pow(angle,3.0))
		+ ((9.21721e-05)*pow(angle,4.0));
	
	double muMC = 0.549863
		+ ((-0.000510435)*pow(angle,1.0))
		+ (( 0.000512479)*pow(angle,2.0))
		+ ((-0.000198997)*pow(angle,3.0))
		+ (( 2.42066e-05)*pow(angle,4.0));
	
	deltaMuGuess -= muMC;
	
	/*
	f_fit->SetParameter(deltaMuPar, 0.002);
	h_data->Fit(f_fit, fitOption);
	
	f_fit->FixParameter(deltaMuPar, f_fit->GetParameter(deltaMuPar));
	*/
	deltaMuGuess = 0.0025;
	f_fit->FixParameter(deltaMuPar, deltaMuGuess);
	
	
	// widen fit range and release background parameters:
	
	f_fit->SetRange(minFitRange, 0.5);
	
	// if we're fitting with the empty target pdf, release that normalization here:
	/*
	if(fitOption_empty) {
		GuessEmptyParameters();
		h_data->Fit(f_fit, fitOption);
	}
	*/
	
	GuessBkgdParameters();
	h_data->Fit(f_fit, fitOption);
	
	f_fit->SetRange(minFitRange, maxFitRange);
	
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
		
		FixBkgdParameters();
		f_fit->SetRange(0.5, maxFitRange);
		
		int     NBkgdPar = f_fit->GetParNumber(     "N_{#eta,bkgd}");
		int    muBkgdPar = f_fit->GetParNumber(   "#mu_{#eta,bkgd}");
		int sigmaBkgdPar = f_fit->GetParNumber("#sigma_{#eta,bkgd}");
		int alphaBkgdPar = f_fit->GetParNumber("#alpha_{#eta,bkgd}");
		int     nBkgdPar = f_fit->GetParNumber(     "n_{#eta,bkgd}");
		
		double NBkgdGuess = f_fit->GetParameter("N_{#eta}") * m_etaPionBkgdFrac;
		
		f_fit->ReleaseParameter(    NBkgdPar);
		f_fit->ReleaseParameter(   muBkgdPar);
		f_fit->ReleaseParameter(sigmaBkgdPar);
		f_fit->ReleaseParameter(alphaBkgdPar);
		f_fit->ReleaseParameter(    nBkgdPar);
		
		f_fit->SetParameter(    NBkgdPar, NBkgdGuess);
		f_fit->SetParameter(   muBkgdPar, 0.58);
		f_fit->SetParameter(sigmaBkgdPar, 0.02);
		f_fit->SetParameter(alphaBkgdPar, 1.0);
		f_fit->SetParameter(    nBkgdPar, 2.0);
		
		f_fit->SetParLimits(    NBkgdPar, 0.000,  1.0e5);
		f_fit->SetParLimits(   muBkgdPar, 0.560,  0.585);
		f_fit->SetParLimits(sigmaBkgdPar, 0.010,  0.025);
		f_fit->SetParLimits(alphaBkgdPar, 0.500,  9.999);
		f_fit->SetParLimits(    nBkgdPar, 1.100, 49.999);
		
		double muBkgdGuess = f_hadronicBkgdLineshape->GetMaximumX() + 0.00175;
		f_fit->FixParameter(muBkgdPar, muBkgdGuess);
		
		//f_fit->FixParameter(alphaBkgdPar, 1.e6);
		//f_fit->FixParameter(    nBkgdPar, 2.0);
		
		h_data->Fit(f_fit, fitOption);
		
		f_fit->SetRange(minFitRange, maxFitRange);
		ReleaseBkgdParameters();
		h_data->Fit(f_fit, fitOption);
		
		f_fit->SetRange(minMggCut, maxFitRange);
		FixBkgdParameters();
		FixOmegaParameters();
		h_data->Fit(f_fit, fitOption);
	}
	else if((fitOption_signal==7) || (fitOption_signal==8)) {
		
		FixBkgdParameters();
		
		f_fit->SetRange(0.5, maxFitRange);
		
		if(!useRawMass) {
			f_fit->ReleaseParameter(hadronicBkgdFracPar);
			f_fit->SetParLimits(hadronicBkgdFracPar, 0.0, 3.0);
			h_data->Fit(f_fit, fitOption);
		}
		
		//---------------------------------//
		
		f_fit->SetRange(minFitRange, maxFitRange);
		ReleaseBkgdParameters();
		h_data->Fit(f_fit, fitOption);
		
		// If at this point the eta+pion fraction is near 0, fix it as such to get a better estimation on the signal yield uncertainty:
		if(f_fit->GetParameter(hadronicBkgdFracPar)<0.001) {
			f_fit->FixParameter(hadronicBkgdFracPar, 0.0);
			FixBkgdParameters();
			h_data->Fit(f_fit, fitOption);
		}
		h_data->Fit(f_fit, fitOption);
		
		f_fit->SetRange(minMggCut, maxFitRange);
		FixBkgdParameters();
		FixOmegaParameters();
		h_data->Fit(f_fit, fitOption);
	}
	else if((fitOption_signal==9) || (fitOption_signal==10)) {
		
		FixBkgdParameters();
		
		f_fit->SetRange(0.5, maxFitRange);
		
		if(!useRawMass) {
			f_fit->ReleaseParameter(etaPionBkgdFracPar);
			f_fit->SetParLimits(etaPionBkgdFracPar, 0.0, 3.0);
			
			if(m_hadronicBkgdYieldBGGEN>10.0) {
				f_fit->ReleaseParameter(hadronicBkgdFracPar);
				f_fit->SetParLimits(hadronicBkgdFracPar, 0.0, 3.0);
			}
			h_data->Fit(f_fit, fitOption);
		}
		
		//---------------------------------//
		
		f_fit->SetRange(minFitRange, maxFitRange);
		ReleaseBkgdParameters();
		h_data->Fit(f_fit, fitOption);
		
		// If at this point the eta+pion fraction is near 0, fix it as such to get a better estimation on the signal yield uncertainty:
		if(f_fit->GetParameter(etaPionBkgdFracPar)<0.01) {
			f_fit->FixParameter(etaPionBkgdFracPar, 0.0);
			FixBkgdParameters();
			h_data->Fit(f_fit, fitOption);
		}
		if(f_fit->GetParameter(hadronicBkgdFracPar)<0.01) {
			f_fit->FixParameter(hadronicBkgdFracPar, 0.0);
			FixBkgdParameters();
			h_data->Fit(f_fit, fitOption);
		}
		h_data->Fit(f_fit, fitOption);
		
		f_fit->SetRange(minMggCut, maxFitRange);
		FixBkgdParameters();
		FixOmegaParameters();
		h_data->Fit(f_fit, fitOption);
	}
	else if(fitOption_signal>10) {
		
		FixBkgdParameters();
		FixOmegaParameters();
		
		//f_fit->ReleaseParameter(deltaMuPar);
		//f_fit->SetParLimits(deltaMuPar, 0.00, 0.01);
		
		f_fit->SetRange(0.5, maxFitRange);
		
		if(!useRawMass) {
			etaPionBkgdFracPar = f_fit->GetParNumber("A_{#eta#pi}");
			hadronicBkgdFracPar = f_fit->GetParNumber("A_{#eta#pi#pi}");
			
			printf("\n\nETAPIONYIELD: %f\n", m_etaPionYieldBGGEN);
			printf("HADRONICBKGDYIELD: %f\n", m_hadronicBkgdYieldBGGEN);
			printf("ETAYIELD: %f\n\n", f_fit->GetParameter("N_{#eta}"));
			
			f_fit->ReleaseParameter(etaPionBkgdFracPar);
			f_fit->SetParameter(etaPionBkgdFracPar, 1.0);
			f_fit->SetParLimits(etaPionBkgdFracPar, 0.0, 10.0);
			
			if(m_hadronicBkgdYieldBGGEN>(0.02*f_fit->GetParameter("N_{#eta}"))) {
			//if(m_hadronicBkgdYieldBGGEN>10.0) {
				f_fit->ReleaseParameter(hadronicBkgdFracPar);
				f_fit->SetParameter(hadronicBkgdFracPar, 1.0);
				f_fit->SetParLimits(hadronicBkgdFracPar, 0.0, 10.0);
			}
			
			if(vetoOption>=6) {
				f_fit->FixParameter(hadronicBkgdFracPar, 1.0);
			}
			
			h_data->Fit(f_fit, fitOption);
			
			//f_fit->FixParameter(deltaMuPar, f_fit->GetParameter(deltaMuPar));
		}
		
		//---------------------------------//
		
		f_fit->SetRange(minFitRange, maxFitRange);
		
		ReleaseBkgdParameters();
		if(maxFitRange > 0.70) {
			ReleaseOmegaParameters();
		}
		h_data->Fit(f_fit, fitOption);
		
		// If at this point the eta+pion fraction is near 0, fix it as such to get a better estimation on the signal yield uncertainty:
		/*
		if(f_fit->GetParameter(etaPionBkgdFracPar)<0.01) {
			f_fit->FixParameter(etaPionBkgdFracPar, 0.0);
			FixBkgdParameters();
			h_data->Fit(f_fit, fitOption);
		}
		if(f_fit->GetParameter(hadronicBkgdFracPar)<0.01) {
			f_fit->FixParameter(hadronicBkgdFracPar, 0.0);
			FixBkgdParameters();
			h_data->Fit(f_fit, fitOption);
		}
		
		h_data->Fit(f_fit, fitOption);
		*/
		
		double locMax = maxFitRange > 0.7 ? 0.7 : maxFitRange;
		f_fit->SetRange(minFitRange, locMax);
		
		FixBkgdParameters();
		FixOmegaParameters();
		h_data->Fit(f_fit, fitOption);
	}
	
	f_fit->SetRange(0.30, 1.20);
	
	return;
}

//==============================================================//
// Omega:

void MggFitter::GuessOmegaParameters()
{
	if(fitOption_omega<=0) return;
	
	double minOmegaFit = 0.68;
	double maxOmegaFit = 0.85;
	
	double     NGuess = h_data->Integral(h_data->FindBin(minOmegaFit), h_data->FindBin(maxOmegaFit));
	double    muGuess = EtaAnalyzer::m_massOmega;
	double sigmaGuess = 0.032;
	double alphaGuess = 1.0;
	double     nGuess = 6.0;
	
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
			f_fit->FixParameter(sigmaPar, 0.0325);
			f_fit->FixParameter(alphaPar, 1.0);
			f_fit->FixParameter(nPar, exp(1.87475 - 0.47071*angle));
			break;
		}
		case 2:
		{
			// Using functional parameterization of MC lineshape:
			
			int   NPar = f_fit->GetParNumber("N_{#omega}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(  NPar, NGuess);
			f_fit->SetParameter(dmuPar, 0.01);
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e6);
			f_fit->SetParLimits(dmuPar, -0.02, 0.02);
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
			f_fit->SetParameter(dmuPar, 0.005);
			
			f_fit->SetParLimits(  NPar,  0.00, 1.e5);
			f_fit->SetParLimits(dmuPar, -0.02, 0.02);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 4:
		{
			// Lineshape fit:
			
			int NomegaPar = f_fit->GetParNumber("N_{#omega}");
			int   NrhoPar = f_fit->GetParNumber("N_{#rho}");
			int    dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->SetParameter(NomegaPar, NGuess);
			f_fit->SetParameter(  NrhoPar, 0.0);
			f_fit->SetParameter(   dmuPar, 0.01);
			
			f_fit->SetParLimits(NomegaPar,  0.00, 1.e5);
			f_fit->SetParLimits(  NrhoPar,  0.00, 1.e5);
			f_fit->SetParLimits(   dmuPar, -0.02, 0.02);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 5:
		{
			// Using lineshape fit result as initial guess
			
			// Let number of omega's and center of peak float freely:
			int  NPar = f_fit->GetParNumber(  "N_{#omega}");
			int muPar = f_fit->GetParNumber("#mu_{#omega}");
			
			f_fit->SetParameter( NPar, NGuess);
			f_fit->SetParLimits( NPar,  0.00, 1.e5);
			f_fit->SetParameter(muPar, f_omegaLineshape->GetParameter(0));
			f_fit->SetParLimits(muPar,  0.76, 0.80);
			
			// for all other parameters, fix according to lineshape:
			for(int ipar=2; ipar<10; ipar++) {
				f_fit->FixParameter(NPar+ipar, f_omegaLineshape->GetParameter(ipar-1));
			}
			break;
		}
	}
	return;
}

void MggFitter::FixOmegaParameters()
{
	switch(fitOption_omega) {
		case 0:
			break;
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
		case 4:
		{
			int NomegaPar = f_fit->GetParNumber("N_{#omega}");
			int   NrhoPar = f_fit->GetParNumber("N_{#rho}");
			int    dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->FixParameter(NomegaPar, f_fit->GetParameter(NomegaPar));
			f_fit->FixParameter(  NrhoPar, f_fit->GetParameter(  NrhoPar));
			f_fit->FixParameter(   dmuPar, f_fit->GetParameter(   dmuPar));
			break;
		}
		case 5:
		{
			int NPar = f_fit->GetParNumber("N_{#omega}");
			for(int ipar=0; ipar<10; ipar++) {
				f_fit->FixParameter(NPar+ipar, f_fit->GetParameter(NPar+ipar));
			}
			break;
		}
	}
	return;
}

void MggFitter::ReleaseOmegaParameters()
{
	switch(fitOption_omega) {
		case 0:
			break;
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
			
			//f_fit->ReleaseParameter(sigmaPar);
			//f_fit->SetParLimits(sigmaPar, 0.015, 0.050);
			
			//f_fit->ReleaseParameter(alphaPar);
			//f_fit->SetParLimits(alphaPar, 0.200, 9.999);
			
			//f_fit->ReleaseParameter(nPar);
			//f_fit->SetParLimits(nPar, 1.100, 49.999);
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
			f_fit->SetParLimits(dmuPar, -0.02, 0.02);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 4:
		{
			// Lineshape fit:
			
			int NomegaPar = f_fit->GetParNumber("N_{#omega}");
			int   NrhoPar = f_fit->GetParNumber("N_{#rho}");
			int    dmuPar = f_fit->GetParNumber("#Delta#mu_{#omega}");
			
			f_fit->ReleaseParameter(NomegaPar);
			f_fit->ReleaseParameter(  NrhoPar);
			f_fit->ReleaseParameter(   dmuPar);
			
			f_fit->SetParameter(NomegaPar, f_fit->GetParameter(NomegaPar));
			f_fit->SetParameter(  NrhoPar, f_fit->GetParameter(  NrhoPar));
			f_fit->SetParameter(   dmuPar, f_fit->GetParameter(   dmuPar));
			
			f_fit->SetParLimits(NomegaPar,  0.00, 1.e5);
			f_fit->SetParLimits(  NrhoPar,  0.00, 1.e5);
			f_fit->SetParLimits(   dmuPar, -0.02, 0.02);
			/*
			if(m_FIX_OMEGA_PARS) {
				f_fit->FixParameter(dmuPar, m_OMEGA_DMU);
			}
			*/
			break;
		}
		case 5:
		{
			int      NPar = f_fit->GetParNumber("N_{#omega}");
			int sigma1Par = f_fit->GetParNumber("#sigma_{#omega,1}");
			int sigma2Par = f_fit->GetParNumber("#sigma_{#omega,2}");
			
			f_fit->ReleaseParameter(NPar);
			f_fit->SetParLimits(NPar, 0., 1.e6);
			
			f_fit->ReleaseParameter(sigma1Par);
			f_fit->SetParLimits(sigma1Par, 0.02, 0.05);
			f_fit->ReleaseParameter(sigma2Par);
			f_fit->SetParLimits(sigma1Par, 0.02, 0.05);
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
		p0Guess =  h_data->GetBinContent(h_data->FindBin(minFitRange)) - f_fit->Eval(minFitRange);
		p1Guess =  minFitRange;
		p2Guess = -0.63;
		p3Guess =  0.00;
	}
	
	switch(fitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->SetParameter(locParIndex, 0.0);
				f_fit->SetParLimits(locParIndex, -1.e6, 1.e6);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			
			f_fit->SetParameter(p0Par, p0Guess);
			f_fit->FixParameter(p1Par, p1Guess);
			f_fit->SetParameter(p2Par, p2Guess);
			f_fit->SetParameter(p3Par, p3Guess);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e6);
			//f_fit->SetParLimits(p1Par,  0.00, 1.00);
			f_fit->SetParLimits(p2Par, -1.e3, 0.0);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->SetParameter(locParIndex, 0.0001);
				f_fit->SetParLimits(locParIndex, -1.e5, 1.e5);
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
				f_fit->SetParLimits(locParIndex, -1.e6, 1.e6);
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
			//f_fit->ReleaseParameter(p1Par);
			f_fit->ReleaseParameter(p2Par);
			f_fit->ReleaseParameter(p3Par);
			
			f_fit->SetParLimits(p0Par,  0.00, 1.e6);
			//f_fit->SetParLimits(p1Par,  0.00, 1.00);
			f_fit->SetParLimits(p2Par, -1.e3, -0.0);
			f_fit->SetParLimits(p3Par, -1.e3, 1.e3);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=fitOption_poly; ipar++) {
				int locParIndex = f_fit->GetParNumber(Form("p%d",ipar));
				f_fit->ReleaseParameter(locParIndex);
				f_fit->SetParLimits(locParIndex, -1.e5, 1.e5);
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
			
			int      NPar = f_fit->GetParNumber("N_{#eta}");
			int   fracPar = f_fit->GetParNumber("fraction_{#eta}");
			int    mu1Par = f_fit->GetParNumber("#mu_{#eta,1}");
			int    dmuPar = f_fit->GetParNumber("#mu_{#eta,2}-#mu_{#eta,1}");
			int sigma1Par = f_fit->GetParNumber("#sigma_{#eta,1}");
			int sigma2Par = f_fit->GetParNumber("#sigma_{#eta,2}");
			
			f_fit->SetParameter(     NPar, NGuess);
			f_fit->SetParameter(  fracPar, 0.2);
			f_fit->SetParameter(   mu1Par, muGuess);
			f_fit->FixParameter(   dmuPar, 0.0);
			f_fit->SetParameter(sigma1Par, sigmaGuess);
			f_fit->SetParameter(sigma2Par, sigmaGuess*2.0);
			
			f_fit->SetParLimits(     NPar,  0.00, 1.e5);
			f_fit->SetParLimits(  fracPar,  0.00, 1.00);
			f_fit->SetParLimits(   mu1Par,  0.54, 0.62);
			//f_fit->SetParLimits(   dmuPar,  0.00, 0.03);
			f_fit->SetParLimits(sigma1Par,  0.01, 0.03);
			f_fit->SetParLimits(sigma2Par,  0.01, 0.05);
			/*
			if(angle<0.5) {
				f_fit->FixParameter(  fracPar, 0.0);
				f_fit->FixParameter(   dmuPar, 0.0);
				f_fit->FixParameter(sigma2Par, 2.0*sigmaGuess);
			}
			*/
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
			
			int   NPar = f_fit->GetParNumber(        "N_{#eta}");
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			f_fit->SetParameter( NPar, NGuess);
			f_fit->SetParLimits( NPar, 0.00, 1.e5);
			//f_fit->SetParameter(dmuPar, 0.0);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.02);
			
			/*
			int     NBkgdPar = f_fit->GetParNumber(     "N_{#eta,bkgd}");
			int    muBkgdPar = f_fit->GetParNumber(   "#mu_{#eta,bkgd}");
			int sigmaBkgdPar = f_fit->GetParNumber("#sigma_{#eta,bkgd}");
			int alphaBkgdPar = f_fit->GetParNumber("#alpha_{#eta,bkgd}");
			int     nBkgdPar = f_fit->GetParNumber(     "n_{#eta,bkgd}");
			
			double NBkgdGuess = NGuess * m_etaPionBkgdFrac;
			
			f_fit->SetParameter(    NBkgdPar, NBkgdGuess);
			f_fit->SetParameter(   muBkgdPar, 0.58);
			f_fit->SetParameter(sigmaBkgdPar, 0.02);
			f_fit->SetParameter(alphaBkgdPar, 1.0);
			f_fit->SetParameter(    nBkgdPar, 2.0);
			
			f_fit->SetParLimits(    NBkgdPar, 0.000,  1.0e5);
			f_fit->SetParLimits(   muBkgdPar, 0.560,  0.585);
			f_fit->SetParLimits(sigmaBkgdPar, 0.010,  0.050);
			f_fit->SetParLimits(alphaBkgdPar, 0.200,  9.999);
			f_fit->SetParLimits(    nBkgdPar, 1.100, 49.999);
			
			f_fit->FixParameter(alphaBkgdPar, 1.e6);
			*/
			break;
		}
		case 7:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
			break;
		}
		case 8:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
			break;
		}
		case 9:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
			break;
		}
		case 10:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
			break;
		}
		case 11:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
			break;
		}
		case 12:
		{
			// Start with Lineshape fit:
			
			int NPar = f_fit->GetParNumber("N_{#eta}");
			f_fit->SetParameter(NPar, NGuess);
			f_fit->SetParLimits(NPar,  0.00, 1.e5);
			
			int dmuPar = f_fit->GetParNumber("#Delta#mu_{#eta}");
			
			// Fixed by a sloped line:
			//f_fit->FixParameter(dmuPar, 0.00127989 + 0.000870044*angle);
			
			// Floating:
			//f_fit->SetParameter(dmuPar, 0.0020);
			//f_fit->SetParLimits(dmuPar, -0.01, 0.03);
			
			// Fixed at a flat value:
			f_fit->FixParameter(dmuPar, 0.00);
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
	
	int shiftPar = f_fit->GetParNumber("shift_{empty}");
	f_fit->SetParameter(shiftPar, 0.0045);
	f_fit->SetParLimits(shiftPar, 0.0, 0.01);
	return;
}

void MggFitter::FixEmptyParameters()
{
	if(fitOption_empty!=1) return;
	
	int NPar = f_fit->GetParNumber("N_{empty}");
	f_fit->FixParameter(NPar, f_fit->GetParameter(NPar));
	
	int shiftPar = f_fit->GetParNumber("shift_{empty}");
	f_fit->FixParameter(shiftPar, f_fit->GetParameter(shiftPar));
	
	return;
}

//==============================================================//
// Lineshape Fits:

void MggFitter::FitEtaLineshape(int drawOption)
{
	if(h_etaLineshape==NULL) return;
	
	// first, check that the bin width of h_etaLineshape matches 'binSize':
	//CheckBinSize(h_etaLineshape, "Eta Lineshape");
	/*
	for(int ibin=1; ibin<=h_etaLineshape->GetXaxis()->GetNbins(); ibin++) {
		double locX = h_etaLineshape->GetBinCenter(ibin);
		if((locX<minMggCut) || (locX>maxMggCut)) {
			h_etaLineshape->SetBinContent(ibin, 0.0);
			h_etaLineshape->SetBinError(ibin, 0.0);
		}
	}
	h_etaLineshape->Scale(1.0/h_etaLineshape->Integral());
	*/
	
	h_etaLineshape->SetMarkerStyle(8);
	h_etaLineshape->SetMarkerSize(0.7);
	h_etaLineshape->SetMarkerColor(kBlue);
	h_etaLineshape->SetLineColor(kBlue);
	h_etaLineshape->GetYaxis()->SetTitle(Form("Normalized Counts / %d MeV/c^{2}", (int)(1.e3 * h_etaLineshape->GetBinWidth(1))));
	h_etaLineshape->GetYaxis()->SetTitleSize(0.05);
	h_etaLineshape->GetYaxis()->SetTitleOffset(1.2);
	h_etaLineshape->GetYaxis()->CenterTitle(true);
	h_etaLineshape->GetXaxis()->SetTitleSize(0.05);
	h_etaLineshape->GetXaxis()->SetTitleOffset(1.0);
	h_etaLineshape->GetXaxis()->CenterTitle(true);
	h_etaLineshape->SetTitle("");
	
	switch(useRawMass) {
		case 0:
		{
			// Two Crystal Ball Functions + One Gaussian:
			
			TF1 *locfEta = new TF1("locfEta", DoubleCrystalBallPlusGausPDF, 0.30, 0.80, 13);
			locfEta->SetParameter(0, 0.540);  // mu1
			locfEta->SetParameter(1, 0.0075); // sigma1
			locfEta->SetParameter(2, 1.400);  // alpha1
			locfEta->SetParameter(3, 9.000);  // n1
			locfEta->SetParameter(4, 0.015);  // mu2-mu1
			locfEta->SetParameter(5, 0.010);  // sigma2
			locfEta->SetParameter(6, 1.400);  // alpha2
			locfEta->SetParameter(7, 13.50);  // n2
			locfEta->SetParameter(8, 0.535);  // mu3
			locfEta->SetParameter(9, 0.007);  // sigma3
			locfEta->SetParameter(10, 0.63);   // fraction1
			locfEta->SetParameter(11, 0.20);   // fraction2
			locfEta->FixParameter(12, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParName(0, "#mu_{1}");
			locfEta->SetParName(1, "#sigma_{1}");
			locfEta->SetParName(2, "#alpha_{1}");
			locfEta->SetParName(3, "n_{1}");
			locfEta->SetParName(4, "#mu_{2}-#mu_{1}");
			locfEta->SetParName(5, "#sigma_{2}");
			locfEta->SetParName(6, "#alpha_{2}");
			locfEta->SetParName(7, "n_{2}");
			locfEta->SetParName(8, "#mu_{3}-#mu_{1}");
			locfEta->SetParName(9, "#sigma_{3}");
			locfEta->SetParName(10, "frac1");
			locfEta->SetParName(11, "frac2");
			
			locfEta->SetParLimits(0, 0.530,  0.560);
			locfEta->SetParLimits(1, 0.004,  0.050);
			locfEta->SetParLimits(2, 0.200,  9.999);
			locfEta->SetParLimits(3, 1.100, 49.999);
			
			locfEta->SetParLimits(4,-0.050,  0.050);
			locfEta->SetParLimits(5, 0.004,  0.050);
			locfEta->SetParLimits(6, 0.200,  9.999);
			locfEta->SetParLimits(7, 1.100, 49.999);
			
			locfEta->SetParLimits(10, 0.0, 1.0);
			locfEta->SetParLimits(11, 0.0, 1.0);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(4, 0.0);
			//locfEta->FixParameter(8, 0.0);
			
			locfEta->FixParameter(8, 0.0);
			locfEta->FixParameter(9, 0.007);
			locfEta->FixParameter(11, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0Q");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				cEtaLS->SetLeftMargin(0.13); cEtaLS->SetRightMargin(0.07);
				cEtaLS->SetBottomMargin(0.13); cEtaLS->SetTopMargin(0.07);
				//gStyle->SetOptFit(0);
				//h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.5,0.6);
				h_etaLineshape->Draw();
				locfEta->SetRange(0.4,0.7);
				locfEta->SetNpx(1000);
				locfEta->Draw("same");
				//h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
				printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
				
				TLatex locLat;
				locLat.SetTextFont(42);
				locLat.SetTextSize(0.045);
				
				double locMinAngle = angle - angleWidth;
				double locMaxAngle = angle + angleWidth;
				
				locLat.SetTextColor(kBlue);
				locLat.DrawLatexNDC(0.155, 0.875, 
					Form("#scale[1.0]{%.2f#circ < #theta_{#gamma#gamma} < %.2f#circ}", locMinAngle, locMaxAngle));
				
				locLat.SetTextColor(kBlack);
				locLat.DrawLatexNDC(0.165, 0.725, Form("#int_{%.2f}^{%.2f}#color[632]{f_{#eta}}dm_{#gamma#gamma} = %.4f", 
					minMggCut, maxMggCut, fracAccepted));
				
				cEtaLS->Update();
				getchar();
				
				/*
				gStyle->SetOptFit(0);
				cEtaLS->SetLogy();
				cEtaLS->Update();
				//cEtaLS->SaveAs(Form("mgg_lineshape_%.2fdeg_%.2fdeg_log.pdf", locMinAngle, locMaxAngle));
				getchar();
				*/
				delete cEtaLS;
			}
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleCrystalBallPlusGausPDF, minFitRange, maxFitRange, 13);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(12, 1.0);
			
			// Widen the lineshape by 10% (test):
			//f_etaLineshape->SetParameter(1, 1.1*f_etaLineshape->GetParameter(1));
			//f_etaLineshape->SetParameter(5, 1.1*f_etaLineshape->GetParameter(5));
			
			locfEta->Delete();
			break;
		}
		case 2:
		{
			// Two Gaussian Functions:
			
			TF1 *locfEta = new TF1("locfEta", DoubleGausPDF, 0.52, 0.59, 6);
			locfEta->SetParameters(
				0.547, //      mu1
				0.010, //   sigma1
				0.000, //  mu2-mu1
				0.015, //   sigma2
				0.300  // fraction
			);
			locfEta->FixParameter(5, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParLimits(0, 0.530,  0.570);
			locfEta->SetParLimits(1, 0.005,  0.050);
			locfEta->SetParLimits(2,-0.050,  0.050);
			locfEta->SetParLimits(3, 0.010,  0.050);
			locfEta->SetParLimits(4, 0.000,  1.000);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(2, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0QL");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				locfEta->SetRange(0.4,0.7);
				h_etaLineshape->Draw();
				locfEta->Draw("same");
				h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
				printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
				
				cEtaLS->Update();
				getchar();
				cEtaLS->SetLogy();
				cEtaLS->Update();
				getchar();
				delete cEtaLS;
			}
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleGausPDF, minFitRange, maxFitRange, 6);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(5, 1.0);
			
			locfEta->Delete();
			break;
		}
		case 1:
		{
			// Two Gaussian Functions:
			
			TF1 *locfEta = new TF1("locfEta", DoubleGausPDF, 0.40, 0.70, 6);
			locfEta->SetParameters(
				0.547, //      mu1
				0.015, //   sigma1
				0.000, //  mu2-mu1
				0.025, //   sigma2
				0.300  // fraction
			);
			locfEta->FixParameter(5, h_etaLineshape->GetXaxis()->GetBinWidth(1));
			
			locfEta->SetParLimits(0, 0.530,  0.570);
			locfEta->SetParLimits(1, 0.005,  0.050);
			locfEta->SetParLimits(2,-0.050,  0.050);
			locfEta->SetParLimits(3, 0.010,  0.050);
			locfEta->SetParLimits(4, 0.000,  1.000);
			
			// Force both crystal ball shapes to have the same mean:
			//locfEta->FixParameter(2, 0.0);
			
			h_etaLineshape->Fit(locfEta, "R0QL");
			
			if(drawOption) {
				TCanvas *cEtaLS = new TCanvas("cEtaLS", "cEtaLS", 950, 700);
				h_etaLineshape->GetXaxis()->SetRangeUser(0.4,0.7);
				locfEta->SetRange(0.4,0.7);
				h_etaLineshape->Draw();
				locfEta->Draw("same");
				h_etaLineshape->Draw("same");
				
				// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
				double fracAccepted = locfEta->Integral(0.5,0.6) / h_etaLineshape->GetXaxis()->GetBinWidth(1);
				printf("  fraction of signal lineshape within mgg cut: %f\n", fracAccepted);
				
				cEtaLS->Update();
				getchar();
				cEtaLS->SetLogy();
				cEtaLS->Update();
				getchar();
				delete cEtaLS;
			}
			
			f_etaLineshape = new TF1("f_etaLineshape", DoubleGausPDF, minFitRange, maxFitRange, 6);
			f_etaLineshape->SetParameters(locfEta->GetParameters());
			f_etaLineshape->FixParameter(5, 1.0);
			
			locfEta->Delete();
			break;
		}
	}
	return;
}

void MggFitter::FitHadronicBkgdLineshape(int drawOption)
{
	if(h_hadronicBkgdLineshape==NULL || m_hadronicBkgdYieldBGGEN<1.0) {
		f_hadronicBkgdLineshape = new TF1("f_hadronicBkgdLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
		f_hadronicBkgdLineshape->SetParameters(0.57, 0.015, 1.0, 5.0, 0.02, 0.020, 1.0, 5.0, 0.0, 1.0);
		return;
	}
	
	TF1 *fHadronicBkgd1 = new TF1("fHadronicBkgd1", CrystalBallPDF_flip, 0.5, 0.7, 5);
	
	double    muGuess = h_hadronicBkgdLineshape->GetBinCenter(h_hadronicBkgdLineshape->GetMaximumBin());
	double sigmaGuess = 0.015;
	double alphaGuess = 1.0;
	double     nGuess = 2.0;
	
	fHadronicBkgd1->SetParameter(0,    muGuess);
	fHadronicBkgd1->SetParameter(1, sigmaGuess);
	fHadronicBkgd1->SetParameter(2, alphaGuess);
	fHadronicBkgd1->SetParameter(3,     nGuess);
	
	fHadronicBkgd1->SetParLimits(0, 0.540,  0.610);
	fHadronicBkgd1->SetParLimits(1, 0.010,  0.050);
	fHadronicBkgd1->SetParLimits(2, 0.200,  9.999);
	fHadronicBkgd1->SetParLimits(3, 1.100, 49.999);
	
	fHadronicBkgd1->FixParameter(4, h_hadronicBkgdLineshape->GetBinWidth(1));
	
	h_hadronicBkgdLineshape->Fit(fHadronicBkgd1, "R0QL");
	
	TF1 *fHadronicBkgd2 = new TF1("fHadronicBkgd2", DoubleCrystalBallPDF_flip, 0.5, 0.75, 10);
	fHadronicBkgd2->SetParameter(0, fHadronicBkgd1->GetParameter(0));
	fHadronicBkgd2->SetParameter(1, fHadronicBkgd1->GetParameter(1));
	fHadronicBkgd2->SetParameter(2, fHadronicBkgd1->GetParameter(2));
	fHadronicBkgd2->SetParameter(3, fHadronicBkgd1->GetParameter(3));
	fHadronicBkgd2->SetParameter(4, 0.02);
	fHadronicBkgd2->SetParameter(5, fHadronicBkgd1->GetParameter(1)*2.0);
	fHadronicBkgd2->SetParameter(6, fHadronicBkgd1->GetParameter(2));
	fHadronicBkgd2->SetParameter(7, fHadronicBkgd1->GetParameter(3));
	fHadronicBkgd2->SetParameter(8, 0.0);
	
	fHadronicBkgd2->SetParLimits(0,  0.540,  0.610);
	fHadronicBkgd2->SetParLimits(1,  0.010,  0.050);
	fHadronicBkgd2->SetParLimits(2,  0.200,  9.999);
	fHadronicBkgd2->SetParLimits(3,  1.100, 49.999);
	fHadronicBkgd2->SetParLimits(4, -0.050,  0.050);
	fHadronicBkgd2->SetParLimits(5,  0.010,  0.050);
	fHadronicBkgd2->SetParLimits(6,  0.200,  9.999);
	fHadronicBkgd2->SetParLimits(7,  1.100, 49.999);
	fHadronicBkgd2->SetParLimits(8,  0.000,  1.000);
	
	fHadronicBkgd2->FixParameter(9, h_hadronicBkgdLineshape->GetBinWidth(1));
	
	h_hadronicBkgdLineshape->Fit(fHadronicBkgd2, "R0Q");
	
	if(drawOption) {
		TCanvas *cHadronicBkgdLS = new TCanvas("cHadronicBkgdLS", "cHadronicBkgdLS", 950, 700);
		//cHadronicBkgdLS->SetLogy();
		h_hadronicBkgdLineshape->Draw();
		fHadronicBkgd2->Draw("same");
		
		// Calculate fraction of PDF between 0.5 and 0.6 GeV/c2:
		double fracAccepted = fHadronicBkgd2->Integral(0.5,0.6) / h_hadronicBkgdLineshape->GetXaxis()->GetBinWidth(1);
		printf("  fraction of hadronic bkgd lineshape within mgg cut: %f\n", fracAccepted);
		
		cHadronicBkgdLS->Update();
		getchar();
		delete cHadronicBkgdLS;
	}
	
	f_hadronicBkgdLineshape = new TF1("f_hadronicBkgdLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
	f_hadronicBkgdLineshape->SetParameters(fHadronicBkgd2->GetParameters());
	f_hadronicBkgdLineshape->FixParameter(9, 1.0);
	
	fHadronicBkgd1->Delete();
	fHadronicBkgd2->Delete();
	
	return;
}

void MggFitter::FitEtaPionLineshape(int drawOption)
{
	if(h_etaPionLineshape==NULL || m_etaPionYieldBGGEN<10.0) {
		f_etaPionLineshape = new TF1("f_etaPionLineshape", DoubleCrystalBallPDF_flip, minFitRange, maxFitRange, 10);
		f_etaPionLineshape->SetParameters(0.57, 0.015, 1.0, 5.0, 0.02, 0.020, 1.0, 5.0, 0.0, 1.0);
		return;
	}
	
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
	
	TF1 *fEtaPion2 = new TF1("fEtaPion2", DoubleCrystalBallPDF_flip, 0.5, 0.75, 10);
	fEtaPion2->SetParameter(0, fEtaPion1->GetParameter(0));
	fEtaPion2->SetParameter(1, fEtaPion1->GetParameter(1));
	fEtaPion2->SetParameter(2, fEtaPion1->GetParameter(2));
	fEtaPion2->SetParameter(3, fEtaPion1->GetParameter(3));
	fEtaPion2->SetParameter(4, 0.02);
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
	
	h_etaPionLineshape->Fit(fEtaPion2, "R0Q");
	
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
	
	fOmega2->SetParName(0, "#mu_{1}");
	fOmega2->SetParName(1, "#sigma_{1}");
	fOmega2->SetParName(2, "#alpha_{1}");
	fOmega2->SetParName(3, "n_{1}");
	fOmega2->SetParName(4, "#mu_{2}-#mu_{1}");
	fOmega2->SetParName(5, "#sigma_{2}");
	fOmega2->SetParName(6, "#alpha_{2}");
	fOmega2->SetParName(7, "n_{2}");
	fOmega2->SetParName(8, "frac");
	
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
	h_omegaLineshape->SetTitle("#gamma+p(n)#rightarrow#omega+p(n)");
	
	if(drawOption) {
		TCanvas *cOmegaLS = new TCanvas("cOmegaLS", "cOmegaLS", 950, 700);
		//cOmegaLS->SetLogy();
		h_omegaLineshape->Draw();
		fOmega2->Draw("same");
		
		cOmegaLS->Update();
		/*
		gStyle->SetOptStat(0);
		cOmegaLS->Modified();
		cOmegaLS->SaveAs("omega_linshape_fit.pdf");
		*/
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

void MggFitter::GetYield(double &yield, double &yieldErr, int useFitPars, int subtractHadronicBkgd) {
	
	// if useFitPars==0 (default), extract yield by integrating counts and subtracting bkgd parameters
	// if useFitPars==1, extract yield by integrating signal lineshape
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut+0.001*binSize);
	int maxMggBin = h_data->FindBin(maxMggCut-0.001*binSize);
	
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
	ZeroSignalPars(locfBkgd, subtractHadronicBkgd);
	
	if(useFitPars==0) {
		double dataCounts = 0.0;
		double  fitCounts = 0.0;
		for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
			double locData = h_data->GetBinContent(ibin);
			double locBkgd = locfBkgd->Eval(h_data->GetBinCenter(ibin));
			double locFit  = f_fit->Eval(h_data->GetBinCenter(ibin));
			fitCounts  += locFit;
			dataCounts += locData;
			yield    += locData - locBkgd;
			yieldErr += pow(h_data->GetBinError(ibin),2.0) + locBkgd;
			//printf("  mgg, data, fit: %f %f %f\n", h_data->GetBinCenter(ibin), locData, locFit);
		}
		yieldErr = sqrt(yieldErr);
		
		printf("Difference between data and fit: %f sigma\n", (dataCounts-fitCounts)/yieldErr);
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
				
				printf("Int1, Int2, locInt = %f, %f, %f\n", locInt1, locInt2, locInt);
				
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
				// Eta Lineshape (PDF) + Eta+Pion Lineshape (PDF):
				
				int yieldPar          = f_fit->GetParNumber("N_{#eta}");
				double N_eta          = f_fit->GetParameter(yieldPar);
				double N_etaErr       = f_fit->GetParError(yieldPar);
				double locCorrection  = f_etaLineshape->Integral(minMggCut, maxMggCut);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				int bkgdPar           = f_fit->GetParNumber("N_{#eta,bkgd}");
				double N_etaBkgd      = f_fit->GetParameter(bkgdPar);
				double N_etaBkgdErr   = f_fit->GetParError(bkgdPar);
				
				TF1 *flocClone;
				InitializeFitFunction(&flocClone,"locf1");
				flocClone->SetParameters(f_fit->GetParameters());
				ZeroSignalPars(flocClone);
				ZeroBkgdPars(flocClone);
				flocClone->SetParameter(bkgdPar, f_fit->GetParameter(bkgdPar));
				
				double locYieldBkgd    = flocClone->Integral(minMggCut,maxMggCut);
				double locYieldBkgdErr = locYieldBkgd * (N_etaBkgdErr/N_etaBkgd);
				
				double locYieldInc    = locYieldEta + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 7:
			{
				// Eta Lineshape (PDF) + Hadronic Lineshape (PDF):
				
				int yieldPar    = f_fit->GetParNumber("N_{#eta}");
				int fractionPar = f_fit->GetParNumber("frac_{bkgd}");
				
				double lsShift  = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double frac_bkgd    = f_fit->GetParameter(fractionPar);
				double frac_bkgdErr = f_fit->GetParError(fractionPar);
				
				double N_bkgd    = N_eta * frac_bkgd;
				double N_bkgdErr = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta' and 'N_etapi' above represent the yield of exclusive eta's and eta+pions integrated over all mgg.
				But to be consistent with our efficiency correction that will be applied later, we need to correct
				this for the small fraction of events that fall outside our mgg cut range.
				
				Recall that at this point f_etaLineshape and f_etaPionLineshape are normalized to have unit-integrals.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrection_bkgd = f_hadronicBkgdLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldHadronicBkgd    = locCorrection_bkgd * N_bkgd;
				double locYieldHadronicBkgdErr = locCorrection_bkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldHadronicBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldHadronicBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 8:
			{
				// Eta Lineshape (PDF) + Hadronic Lineshape (PDF):
				
				int yieldPar    = f_fit->GetParNumber("N_{#eta}");
				int fractionPar = f_fit->GetParNumber("frac_{bkgd}");
				
				double lsShift  = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double frac_bkgd    = f_fit->GetParameter(fractionPar);
				double frac_bkgdErr = f_fit->GetParError(fractionPar);
				
				double N_bkgd    = N_eta * frac_bkgd;
				double N_bkgdErr = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta' and 'N_etapi' above represent the yield of exclusive eta's and eta+pions integrated over all mgg.
				But to be consistent with our efficiency correction that will be applied later, we need to correct
				this for the small fraction of events that fall outside our mgg cut range.
				
				Recall that at this point f_etaLineshape and f_etaPionLineshape are normalized to have unit-integrals.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrection_bkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldHadronicBkgd    = locCorrection_bkgd * N_bkgd;
				double locYieldHadronicBkgdErr = locCorrection_bkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldHadronicBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldHadronicBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 9:
			{
				// Eta (Histogram) + Eta+Pion Background (Histogram) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar         = f_fit->GetParNumber("N_{#eta}");
				int fractionParEtaPi = f_fit->GetParNumber("frac_{#eta#pi}");
				int fractionParBkgd  = f_fit->GetParNumber("frac_{bkgd}");
				
				double lsShift    = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double frac_etapi    = f_fit->GetParameter(fractionParEtaPi);
				double frac_etapiErr = f_fit->GetParError(fractionParEtaPi);
				
				double frac_bkgd     = f_fit->GetParameter(fractionParBkgd);
				double frac_bkgdErr  = f_fit->GetParError(fractionParBkgd);
				
				double N_etapi    = N_eta * frac_etapi;
				double N_etapiErr = sqrt(pow(N_etaErr*frac_etapi,2.0) + pow(N_eta*frac_etapiErr,2.0));
				
				double N_bkgd     = N_eta * frac_bkgd;
				double N_bkgdErr  = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = h_etaLineshape->Integral(h_etaLineshape->FindBin(minMggCut-lsShift), h_etaLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = f_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift), 
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion    = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 10:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar         = f_fit->GetParNumber("N_{#eta}");
				int fractionParEtaPi = f_fit->GetParNumber("frac_{#eta#pi}");
				int fractionParBkgd  = f_fit->GetParNumber("frac_{bkgd}");
				
				double lsShift    = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double frac_etapi    = f_fit->GetParameter(fractionParEtaPi);
				double frac_etapiErr = f_fit->GetParError(fractionParEtaPi);
				
				double frac_bkgd     = f_fit->GetParameter(fractionParBkgd);
				double frac_bkgdErr  = f_fit->GetParError(fractionParBkgd);
				
				double N_etapi    = N_eta * frac_etapi;
				double N_etapiErr = sqrt(pow(N_etaErr*frac_etapi,2.0) + pow(N_eta*frac_etapiErr,2.0));
				
				double N_bkgd     = N_eta * frac_bkgd;
				double N_bkgdErr  = sqrt(pow(N_etaErr*frac_bkgd,2.0) + pow(N_eta*frac_bkgdErr,2.0));
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				//double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locCorrection  = h_etaLineshape->Integral(h_etaLineshape->FindBin(minMggCut-lsShift), 
					h_etaLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift), 
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 11:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar      = f_fit->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_fit->GetParNumber("A_{#eta#pi}");
				int yieldParBkgd  = f_fit->GetParNumber("A_{#eta#pi#pi}");
				
				double lsShift    = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double N_etapi    = f_fit->GetParameter(yieldParEtaPi) * m_etaPionYieldBGGEN;
				double N_etapiErr = f_fit->GetParError(yieldParEtaPi) * m_etaPionYieldBGGEN;
				
				double N_bkgd     = f_fit->GetParameter(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				double N_bkgdErr  = f_fit->GetParError(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = f_etaPionLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
					yield    = locYieldEta;
					yieldErr = locYieldEtaErr;
				} else {
					yield    = locYieldInc;
					yieldErr = locYieldIncErr;
				}
				break;
			}
			case 12:
			{
				// Eta (Lineshape) + Eta+Pion Background (Lineshape) + Other Hadronic Bkgd (Histogram):
				
				int yieldPar      = f_fit->GetParNumber("N_{#eta}");
				int yieldParEtaPi = f_fit->GetParNumber("A_{#eta#pi}");
				int yieldParBkgd  = f_fit->GetParNumber("A_{#eta#pi#pi}");
				
				double lsShift    = f_fit->GetParameter("#Delta#mu_{#eta}");
				
				double N_eta    = f_fit->GetParameter(yieldPar);
				double N_etaErr = f_fit->GetParError(yieldPar);
				
				double N_etapi    = f_fit->GetParameter(yieldParEtaPi) * m_etaPionYieldBGGEN;
				double N_etapiErr = f_fit->GetParError(yieldParEtaPi) * m_etaPionYieldBGGEN;
				
				double N_bkgd     = f_fit->GetParameter(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				double N_bkgdErr  = f_fit->GetParError(yieldParBkgd) * m_hadronicBkgdYieldBGGEN;
				
				/*
				'N_eta', 'N_etapi', and 'N_bkgd' above represent the yield of exclusive eta's, eta+pions, and other 
				hadronic backgrounds integrated over all mgg. But to be consistent with our efficiency correction 
				that will be applied later, we need to correct this for the fraction of events that fall outside our mgg cut range.
				*/
				
				double locCorrection  = f_etaLineshape->Integral(minMggCut-lsShift, maxMggCut-lsShift);
				double locYieldEta    = locCorrection * N_eta;
				double locYieldEtaErr = locCorrection * N_etaErr;
				
				double locCorrectionEtaPion = h_etaPionLineshape->Integral(h_etaPionLineshape->FindBin(minMggCut-lsShift),
					h_etaPionLineshape->FindBin(maxMggCut-lsShift));
				double locYieldEtaPion      = locCorrectionEtaPion * N_etapi;
				double locYieldEtaPionErr   = locCorrectionEtaPion * N_etapiErr;
				
				double locCorrectionBkgd = h_hadronicBkgdLineshape->Integral(h_hadronicBkgdLineshape->FindBin(minMggCut-lsShift), 
					h_hadronicBkgdLineshape->FindBin(maxMggCut-lsShift));
				double locYieldBkgd    = locCorrectionBkgd * N_bkgd;
				double locYieldBkgdErr = locCorrectionBkgd * N_bkgdErr;
				
				double locYieldInc    = locYieldEta + locYieldEtaPion + locYieldBkgd;
				double locYieldIncErr = sqrt(pow(locYieldEtaErr,2.0) + pow(locYieldEtaPion,2.0) + pow(locYieldBkgdErr,2.0));
				
				if(subtractHadronicBkgd) {
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

void MggFitter::GetEmptyYield(double &yield, double &yieldErr, int excludeNonPeaking) {
	
	yield    = 0.0;
	yieldErr = 0.0;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeEmptyFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_empty->GetParameters());
	
	// zero the parameters not associated with peaking structure in eta mass region:
	
	if(excludeNonPeaking){
		for(int ipar=0; ipar<nParameters; ipar++) {
			TString locParName = Form("%s",locfEta->GetParName(ipar));
			//printf(" p%d: %s = %f\n", ipar, locParName.Data(), locfEta->GetParameter(ipar));
			
			if(locParName.Contains("eta") || locParName.Contains("sigma") || locParName.Contains("alpha")) continue;
			if(locParName.Contains("binWidth") || locParName.Contains("angle")) continue;
			if(locParName.Contains("fdc,2")) {
				//
				// Comment out the following line to exclude the 
				// portion of the empty fit function which is meant to 
				// describe the omegas coming from the second FDC package
				// from the yield of etas estimated from the empty target
				// runs.
				// 
				continue;
			}
			
			locfEta->SetParameter(ipar,0.0);
			locfEta->SetParError(ipar,0.0);
		}
	}
	
	// the last parameter of the empty target fit function is the bin size:
	locfEta->SetParameter(locfEta->GetNpar()-1, emptyBinSize);
	
	/*
	// Print out parameter list just to make sure everything was zeroed properly:
	printf("\n");
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		printf(" p%d: %s = %f\n", ipar, locParName.Data(), locfEta->GetParameter(ipar));
	}
	*/
	//-----------------------------------------------//
	
	int minMggBin = h_empty->FindBin(minMggCut);
	int maxMggBin = h_empty->FindBin(maxMggCut)-1;
	
	// integrate empty function:
	
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

void MggFitter::GetOmegaYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	if(fitOption_omega==0) return;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_fit->GetParameters());
	
	// zero the parameters not associated with the omega background:
	
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		if((locParName.Contains("omega")) || (locParName.Contains("sigma"))) continue;
		
		locfEta->SetParameter(ipar,0.0);
		locfEta->SetParError(ipar,0.0);
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut);
	int maxMggBin = h_data->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_data->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int omegaYieldPar = f_fit->GetParNumber("N_{#omega}");
	
	double locRelErr = f_fit->GetParError(omegaYieldPar) / f_fit->GetParameter(omegaYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}



void MggFitter::GetOmegaFitPars(double &mu, double &muErr, double &sigma, double &sigmaErr,
	double &alpha, double &alphaErr, double &n, double &nErr)
{
	if(fitOption_omega!=1) {
		mu       = 0.0;
		muErr    = 0.0;
		sigma    = 0.0;
		sigmaErr = 0.0;
		alpha    = 0.0;
		alphaErr = 0.0;
		n        = 0.0;
		nErr     = 0.0;
		return;
	}
	
	int muPar    = f_fit->GetParNumber("#mu_{#omega}");
	mu           = f_fit->GetParameter(muPar);
	muErr        = f_fit->GetParError(muPar);
	
	int sigmaPar = f_fit->GetParNumber("#sigma_{#omega}");
	sigma        = f_fit->GetParameter(sigmaPar);
	sigmaErr     = f_fit->GetParError(sigmaPar);
	
	int alphaPar = f_fit->GetParNumber("#alpha_{#omega}");
	alpha        = f_fit->GetParameter(alphaPar);
	alphaErr     = f_fit->GetParError(alphaPar);
	
	int nPar     = f_fit->GetParNumber("n_{#omega}");
	n            = f_fit->GetParameter(nPar);
	nErr         = f_fit->GetParError(nPar);
}

void MggFitter::GetBkgdYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	if(fitOption_bkgd==4) return;
	
	//-----------------------------------------------//
	
	TF1 *locfBkgd;
	int nParameters = InitializeFitFunction(&locfBkgd, "locfBkgd");
	locfBkgd->SetParameters(f_fit->GetParameters());
	
	// zero the parameters not associated with the background:
	
	int nBkgdPars = 0;
	switch(fitOption_bkgd) {
		case 1:
			nBkgdPars = fitOption_poly + 1;
			break;
		case 2:
			nBkgdPars = 4;
			break;
		case 3:
			nBkgdPars = fitOption_poly + 1;
			break;
		case 4:
			nBkgdPars = 0;
			break;
	}
	
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfBkgd->GetParName(ipar));
		if(locParName.Contains("sigma")) continue;
		
		bool isBkgdPar = false;
		for(int jpar=0; jpar<nBkgdPars; jpar++) {
			TString locBkgdParName = Form("p%d",jpar);
			if(locParName==locBkgdParName) {
				isBkgdPar = true;
				break;
			}
		}
		if(isBkgdPar) continue;
		
		locfBkgd->SetParameter(ipar,0.0);
		locfBkgd->SetParError(ipar,0.0);
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut);
	int maxMggBin = h_data->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfBkgd->Eval(h_data->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// (placeholder) uncertainty:
	
	yieldErr = sqrt(yield);
	
	delete locfBkgd;
	return;
}

void MggFitter::GetHadronicBkgdYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	
	if(fitOption_signal<6) return;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_fit->GetParameters());
	
	// zero the parameters not associated with the eta+pion background:
	
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		
		if(locParName=="N_{#eta}") {
			if(fitOption_signal<11) continue;
		}
		if(locParName=="#Delta#mu_{#eta}") continue;
		if(locParName=="frac_{bkgd}")      continue;
		if(locParName=="A_{#eta#pi#pi}")   continue;
		if(locParName.Contains("sigma"))   continue;
		if(locParName.Contains("alpha"))   continue;
		
		locfEta->SetParameter(ipar,0.0);
		locfEta->SetParError(ipar,0.0);
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut);
	int maxMggBin = h_data->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_data->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int bkgdYieldPar;
	if(fitOption_signal>10) bkgdYieldPar = f_fit->GetParNumber("A_{#eta#pi#pi}");
	else bkgdYieldPar = f_fit->GetParNumber("frac_{bkgd}");
	
	double locRelErr = f_fit->GetParError(bkgdYieldPar) / f_fit->GetParameter(bkgdYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}

void MggFitter::GetEtaPionYield(double &yield, double &yieldErr)
{
	yield    = 0.0;
	yieldErr = 0.0;
	
	if(fitOption_signal<9) return;
	
	//-----------------------------------------------//
	
	TF1 *locfEta;
	int nParameters = InitializeFitFunction(&locfEta, "locfEta");
	locfEta->SetParameters(f_fit->GetParameters());
	
	// zero the parameters not associated with the eta+pion background:
	
	for(int ipar=0; ipar<nParameters; ipar++) {
		TString locParName = Form("%s",locfEta->GetParName(ipar));
		
		if(locParName=="N_{#eta}") {
			if(fitOption_signal<11) continue;
		}
		if(locParName=="#Delta#mu_{#eta}") continue;
		if(locParName=="frac_{#eta#pi}")   continue;
		if(locParName=="A_{#eta#pi}")      continue;
		if(locParName.Contains("sigma"))   continue;
		if(locParName.Contains("alpha"))   continue;
		
		locfEta->SetParameter(ipar,0.0);
		locfEta->SetParError(ipar,0.0);
	}
	
	//-----------------------------------------------//
	
	int minMggBin = h_data->FindBin(minMggCut);
	int maxMggBin = h_data->FindBin(maxMggCut)-1;
	
	// integrate function:
	
	for(int ibin=minMggBin; ibin<=maxMggBin; ibin++) {
		double locYield = locfEta->Eval(h_data->GetBinCenter(ibin));
		yield += locYield;
	}
	
	// estimate uncertainty from fit parameter:
	
	int etaPionYieldPar;
	if(fitOption_signal>10) etaPionYieldPar = f_fit->GetParNumber("A_{#eta#pi}");
	else etaPionYieldPar = f_fit->GetParNumber("frac_{#eta#pi}");
	
	double locRelErr = f_fit->GetParError(etaPionYieldPar) / f_fit->GetParameter(etaPionYieldPar);
	if(locRelErr>2.0) locRelErr = 2.0;
	
	yieldErr = yield * locRelErr;
	
	delete locfEta;
	return;
}

void MggFitter::GetLineshapeShift(double &shift, double &shiftErr)
{
	shift    = 0.0;
	shiftErr = 0.0;
	
	if(fitOption_signal<6) return;
	
	int ipar = f_fit->GetParNumber("#Delta#mu_{#eta}");
	shift    = f_fit->GetParameter(ipar);
	shiftErr = f_fit->GetParError(ipar);
	return;
}

void MggFitter::ZeroSignalPars(TF1 *f1, int subtractHadronicBkgd)
{
	switch(fitOption_signal) {
		case 1:
			f1->SetParameter("N_{#eta}", 0.0);
			break;
		case 2:
			f1->SetParameter("N_{#eta}", 0.0);
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
			f1->SetParameter("N_{#eta}", 0.0);
			if(!subtractHadronicBkgd) f1->SetParameter("N_{#eta,bkgd}", 0.0);
			break;
		case 7:
			if(subtractHadronicBkgd) {
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
		case 8:
			if(subtractHadronicBkgd) {
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
		case 9:
			if(subtractHadronicBkgd) {
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
		case 10:
			if(subtractHadronicBkgd) {
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
		case 11:
			if(subtractHadronicBkgd) {
				// If we are subtracting the eta+pion background, we
				// don't consider it as 'signal' parameter. So when
				// we zero the signal parameters, we only zero the portion
				// due to exclusive etas.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("A_{#eta#pi}", 0.0);
				f1->SetParameter("A_{#eta#pi#pi}",    0.0);
			}
			break;
		case 12:
			if(subtractHadronicBkgd) {
				// If we are subtracting the eta+pion background, we
				// don't consider it as 'signal' parameter. So when
				// we zero the signal parameters, we only zero the portion
				// due to exclusive etas.
				f1->SetParameter("N_{#eta}", 0.0);
			}
			else {
				// By default, zero everything peaking in the eta mass region.
				f1->SetParameter("N_{#eta}",    0.0);
				f1->SetParameter("A_{#eta#pi}", 0.0);
				f1->SetParameter("A_{#eta#pi#pi}",    0.0);
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
	f1->SetParameter("N_{#rho}",  0.0);
	f1->SetParameter("N_{#eta'}",   0.0);
	f1->SetParameter("N_{empty}",   0.0);
	switch(fitOption_signal) {
		case 6:
			f1->SetParameter("N_{#eta,bkgd}",  0.0);
			break;
		case 7:
			f1->SetParameter("frac_{bkgd}", 0.0);
			break;
		case 8:
			f1->SetParameter("frac_{bkgd}", 0.0);
			break;
		case 9:
			f1->SetParameter("frac_{#eta#pi}", 0.0);
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 10:
			f1->SetParameter("frac_{#eta#pi}", 0.0);
			f1->SetParameter("frac_{bkgd}",    0.0);
			break;
		case 11:
			f1->SetParameter("A_{#eta#pi}", 0.0);
			f1->SetParameter("A_{#eta#pi#pi}",    0.0);
			break;
		case 12:
			f1->SetParameter("A_{#eta#pi}", 0.0);
			f1->SetParameter("A_{#eta#pi#pi}",    0.0);
			break;
	}
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

void MggFitter::GetHadronicBkgdFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	// redundant check:
	if(fitOption_signal<7) return;
	
	int fractionPar;
	if((fitOption_signal==11) || (fitOption_signal==12)) fractionPar = f_fit->GetParNumber("A_{#eta#pi#pi}");
	else fractionPar = f_fit->GetParNumber("frac_{bkgd}");
	
	fraction    = f_fit->GetParameter(fractionPar);
	fractionErr = f_fit->GetParError(fractionPar);
	return;
}

void MggFitter::GetEtaPionBkgdFraction(double &fraction, double &fractionErr)
{
	fraction    = 0.0;
	fractionErr = 0.0;
	
	// redundant check:
	if(fitOption_signal<9) return;
	
	int fractionPar;
	if((fitOption_signal==11) || (fitOption_signal==12)) fractionPar = f_fit->GetParNumber("A_{#eta#pi}");
	else fractionPar = f_fit->GetParNumber("frac_{#eta#pi}");
	
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
