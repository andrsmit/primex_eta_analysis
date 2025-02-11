#include "MggFitter.h"

void MggFitter::FitEmpty()
{
	InitializeEmptyFitFunction(&f_empty);
	
	Option_t *locFitOption = "R0QL";
	
	m_nEmptyParameters = InitializeEmptyFitParameters();
	excludeRegions.clear();
	
	excludeRegions.push_back({0.50,0.59});
	excludeRegions.push_back({0.70,0.85});
	GuessEmptyBkgdParameters();
	f_empty->SetRange(minEmptyFitRange, maxEmptyFitRange);
	h_empty->Fit(f_empty, locFitOption);
	
	//------------------------------//
	// now try to fit eta peak:
	
	if(emptyFitOption_eta>0) {
		excludeRegions.clear();
		FixEmptyBkgdParameters();
		
		GuessEmptyEtaParameters();
		f_empty->SetRange(0.5,0.6);
		h_empty->Fit(f_empty, locFitOption);
		
		excludeRegions.push_back({0.70, 0.85});
		ReleaseEmptyBkgdParameters();
		f_empty->SetRange(minEmptyFitRange, maxEmptyFitRange);
		h_empty->Fit(f_empty, locFitOption);
	}
	
	//------------------------------//
	// now try to fit omega peak:
	
	if(emptyFitOption_omega>0) {
		excludeRegions.clear();
		FixEmptyBkgdParameters();
		FixEmptyEtaParameters();
		GuessEmptyOmegaParameters();
		f_empty->SetRange(0.7,0.85);
		h_empty->Fit(f_empty, locFitOption);
		
		ReleaseEmptyBkgdParameters();
		f_empty->SetRange(minEmptyFitRange, maxEmptyFitRange);
		h_empty->Fit(f_empty, locFitOption);
	}
	
	//------------------------------//
	// now try to fit peaking structures from FDC:
	
	if(emptyFitOption_fdc>0) {
		excludeRegions.clear();
		FixEmptyBkgdParameters();
		FixEmptyEtaParameters();
		FixEmptyOmegaParameters();
		GuessEmptyFDCParameters();
		h_empty->Fit(f_empty, locFitOption);
		
		ReleaseEmptyBkgdParameters();
		h_empty->Fit(f_empty, locFitOption);
	}
	
	//------------------------------//
	
	f_empty->SetRange(minEmptyFitRange, maxEmptyFitRange);
	h_empty->Fit(f_empty, locFitOption);
	
	m_nEmptyEtaPar   = f_empty->GetParameter("N_{#eta}");
	m_nEmptyOmegaPar = f_empty->GetParameter("N_{#omega}");
	
	excludeRegions.clear();
	f_empty->SetRange(minEmptyFitRange, maxEmptyFitRange);
	return;
}

void MggFitter::GuessEmptyBkgdParameters()
{
	if(emptyFitOption_bkgd==4) return;
	
	int p0Par = f_empty->GetParNumber("p0");
	int p1Par = f_empty->GetParNumber("p1");
	int p2Par = f_empty->GetParNumber("p2");
	int p3Par = f_empty->GetParNumber("p3");
	
	double p0Guess = 0.;
	double p1Guess = 0., p2Guess = 0., p3Guess = 0., p4Guess = 0.;
	
	if(emptyFitOption_bkgd==2) {
		p0Guess = h_empty->GetBinContent(h_empty->FindBin(minFitRange)) - f_empty->Eval(minFitRange);
		p1Guess =  0.0;
		p2Guess = -1.0;
		p3Guess =  0.0;
		p4Guess = minFitRange;
	}
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->SetParameter(locParIndex, 0.0);
				f_empty->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p4Par = f_empty->GetParNumber("p4");
			
			f_empty->SetParameter(p0Par, p0Guess);
			f_empty->SetParameter(p1Par, p1Guess);
			f_empty->SetParameter(p2Par, p2Guess);
			f_empty->SetParameter(p3Par, p3Guess);
			f_empty->SetParameter(p4Par, p4Guess);
			
			f_empty->SetParLimits(p0Par,  0.00, 1.e4);
			f_empty->SetParLimits(p1Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p2Par, -1.e3, 0.);
			f_empty->SetParLimits(p3Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p4Par,  0.00, 0.50);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->SetParameter(locParIndex, 0.0001);
				f_empty->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
	}
	return;
}

void MggFitter::FixEmptyBkgdParameters() {
	
	if(emptyFitOption_bkgd==4) return;
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->FixParameter(locParIndex, f_empty->GetParameter(locParIndex));
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p4_par = f_empty->GetParNumber("p4");
			for(int ipar=0; ipar<5; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->FixParameter(locParIndex, f_empty->GetParameter(locParIndex));
			}
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->FixParameter(locParIndex, f_empty->GetParameter(locParIndex));
			}
			break;
		}
	}
	return;
}

void MggFitter::ReleaseEmptyBkgdParameters() {
	
	if(emptyFitOption_bkgd==4) return;
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// 3rd-order polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->ReleaseParameter(locParIndex);
				f_empty->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p0Par = f_empty->GetParNumber("p0");
			int p1Par = f_empty->GetParNumber("p1");
			int p2Par = f_empty->GetParNumber("p2");
			int p3Par = f_empty->GetParNumber("p3");
			int p4Par = f_empty->GetParNumber("p4");
			
			f_empty->ReleaseParameter(p0Par);
			f_empty->ReleaseParameter(p1Par);
			f_empty->ReleaseParameter(p2Par);
			f_empty->ReleaseParameter(p3Par);
			f_empty->ReleaseParameter(p4Par);
			
			f_empty->SetParLimits(p0Par,  0.00, 1.e4);
			f_empty->SetParLimits(p1Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p2Par, -1.e3, 0.);
			f_empty->SetParLimits(p3Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p4Par,  0.00, 0.50);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_empty->GetParNumber(Form("p%d",ipar));
				f_empty->ReleaseParameter(locParIndex);
				f_empty->SetParLimits(locParIndex, -1.e4, 1.e4);
			}
			break;
		}
	}
	return;
}

void MggFitter::GuessEmptyEtaParameters() {
	
	// guess number of eta' by integrating histogram and subtracting background:
	
	double minEtaFit = 0.50;
	double maxEtaFit = 0.60;
	
	double     nEtaGuess = h_empty->Integral(h_empty->FindBin(minEtaFit), h_empty->FindBin(maxEtaFit));
	double    muEtaGuess = 0.530;
	double sigmaEtaGuess = 0.015;
	
	double locEtaMax = 0.0;
	for(int ibin=h_empty->FindBin(minEtaFit); ibin<=h_empty->FindBin(maxEtaFit); ibin++) {
		if(h_empty->GetBinContent(ibin) > locEtaMax) {
			locEtaMax  = h_empty->GetBinContent(ibin);
			muEtaGuess = h_empty->GetBinCenter(ibin);
		}
	}
	
	switch(emptyFitOption_eta) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			int     nEtaPar = f_empty->GetParNumber("N_{#eta}");
			int    muEtaPar = f_empty->GetParNumber("#mu_{#eta}");
			int sigmaEtaPar = f_empty->GetParNumber("#sigma_{#eta}");
			
			f_empty->SetParameter(    nEtaPar,     nEtaGuess);
			f_empty->SetParameter(   muEtaPar,    muEtaGuess);
			f_empty->SetParameter(sigmaEtaPar, sigmaEtaGuess);
			
			f_empty->SetParLimits(    nEtaPar, 0.00, 1.e5);
			f_empty->SetParLimits(   muEtaPar, 0.52, 0.57);
			f_empty->SetParLimits(sigmaEtaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   nEtaPar = f_empty->GetParNumber("N_{#eta}");
			int dmuEtaPar = f_empty->GetParNumber("#Delta#mu_{#eta}");
			
			f_empty->SetParameter(  nEtaPar, nEtaGuess);
			f_empty->SetParameter(dmuEtaPar, 0.005);
			
			f_empty->SetParLimits(  nEtaPar,  0.00, 1.e5);
			f_empty->SetParLimits(dmuEtaPar, -0.05, 0.02);
			break;
		}
	}
	return;
}

void MggFitter::FixEmptyEtaParameters() {
	
	switch(emptyFitOption_eta) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			int     nEtaPar = f_empty->GetParNumber("N_{#eta}");
			int    muEtaPar = f_empty->GetParNumber("#mu_{#eta}");
			int sigmaEtaPar = f_empty->GetParNumber("#sigma_{#eta}");
			
			f_empty->FixParameter(    nEtaPar, f_empty->GetParameter(    nEtaPar));
			f_empty->FixParameter(   muEtaPar, f_empty->GetParameter(   muEtaPar));
			f_empty->FixParameter(sigmaEtaPar, f_empty->GetParameter(sigmaEtaPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   nEtaPar = f_empty->GetParNumber("N_{#eta}");
			int dmuEtaPar = f_empty->GetParNumber("#Delta#mu_{#eta}");
			
			f_empty->FixParameter(  nEtaPar, f_empty->GetParameter(  nEtaPar));
			f_empty->FixParameter(dmuEtaPar, f_empty->GetParameter(dmuEtaPar));
			break;
		}
	}
	return;
}

void MggFitter::GuessEmptyOmegaParameters() {
	
	switch(emptyFitOption_omega) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// Crystal Ball function:
			
			int     nOmegaPar = f_empty->GetParNumber("N_{#omega}");
			int    muOmegaPar = f_empty->GetParNumber("#mu_{#omega}");
			int sigmaOmegaPar = f_empty->GetParNumber("#sigma_{#omega}");
			int alphaOmegaPar = f_empty->GetParNumber("#alpha_{#omega}");
			int    nnOmegaPar = f_empty->GetParNumber("n_{#omega}");
			
			f_empty->SetParameter(    nOmegaPar, 0.000);
			f_empty->SetParameter(   muOmegaPar, 0.760);
			f_empty->SetParameter(sigmaOmegaPar, 0.025);
			f_empty->SetParameter(alphaOmegaPar, 1.000);
			f_empty->SetParameter(   nnOmegaPar, 2.000);
			
			f_empty->SetParLimits(    nOmegaPar, 0.000, 1.0e6);
			f_empty->SetParLimits(   muOmegaPar, 0.740, 0.780);
			f_empty->SetParLimits(sigmaOmegaPar, 0.015, 0.050);
			f_empty->SetParLimits(alphaOmegaPar, 0.500, 9.999);
			f_empty->SetParLimits(   nnOmegaPar, 0.500, 9.999);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   nOmegaPar = f_empty->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_empty->GetParNumber("#Delta#mu_{#omega}");
			
			f_empty->SetParameter(  nOmegaPar,  0.00);
			f_empty->SetParameter(dmuOmegaPar, -0.01);
			
			f_empty->SetParLimits(  nOmegaPar,  0.00, 1.e5);
			f_empty->SetParLimits(dmuOmegaPar, -0.04, 0.03);
			break;
		}
	}
	return;
}

void MggFitter::FixEmptyOmegaParameters() {
	
	switch(emptyFitOption_eta) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// Crystal Ball:
			
			int     nOmegaPar = f_empty->GetParNumber("N_{#omega}");
			int    muOmegaPar = f_empty->GetParNumber("#mu_{#omega}");
			int sigmaOmegaPar = f_empty->GetParNumber("#sigma_{#omega}");
			int alphaOmegaPar = f_empty->GetParNumber("#alpha_{#omega}");
			int    nnOmegaPar = f_empty->GetParNumber("n_{#omega}");
			
			f_empty->FixParameter(    nOmegaPar, f_empty->GetParameter(    nOmegaPar));
			f_empty->FixParameter(   muOmegaPar, f_empty->GetParameter(   muOmegaPar));
			f_empty->FixParameter(sigmaOmegaPar, f_empty->GetParameter(sigmaOmegaPar));
			f_empty->FixParameter(alphaOmegaPar, f_empty->GetParameter(alphaOmegaPar));
			f_empty->FixParameter(   nnOmegaPar, f_empty->GetParameter(   nnOmegaPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   nOmegaPar = f_empty->GetParNumber("N_{#omega}");
			int dmuOmegaPar = f_empty->GetParNumber("#Delta#mu_{#omega}");
			
			f_empty->FixParameter(  nOmegaPar, f_empty->GetParameter(  nOmegaPar));
			f_empty->FixParameter(dmuOmegaPar, f_empty->GetParameter(dmuOmegaPar));
			break;
		}
	}
	return;
}

void MggFitter::GuessEmptyFDCParameters() {
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar   = f_empty->GetParNumber(Form("N_{fdc,%d}",i+1));
				int locMuPar  = f_empty->GetParNumber(Form("#mu_{fdc,%d}",i+1));
				int locSigPar = f_empty->GetParNumber(Form("#sigma_{fdc,%d}",i+1));
				f_empty->SetParameter(locNPar, 0.0);
				f_empty->SetParameter(locMuPar, m_muFDC[i]);
				f_empty->SetParameter(locSigPar, 0.02);
				
				f_empty->SetParLimits(locNPar,   0.0, 1.e5);
				f_empty->SetParLimits(locMuPar,  m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				f_empty->SetParLimits(locSigPar, 0.010, 0.035);
			}
			break;
		}
		case 2:
		{
			int locNPar = f_empty->GetParNumber("N_{fdc}");
			f_empty->SetParameter(locNPar, 0.0);
			f_empty->SetParLimits(locNPar, 0.0, 1.e5);
			for(int i=0; i<m_muFDC.size(); i++) {
				int locDeltaMuPar = f_empty->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_empty->SetParameter(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]);
				f_empty->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar = f_empty->GetParNumber(Form("N_{fdc,%d}",i+1));
				f_empty->SetParameter(locNPar, 0.0);
				f_empty->SetParLimits(locNPar, 0.0, 1.e5);
				
				int locDeltaMuPar = f_empty->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_empty->SetParameter(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]);
				f_empty->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
	}
	return;
}

void MggFitter::DumpEmptyFitParameters()
{
	if(f_empty==NULL || m_nEmptyParameters==0) {
		printf("\n\nNo empty target fit was performed.\n");
		return;
	}
	printf("\n\nEmpty-Target Fit parameters:\n");
	for(int ipar=0; ipar<m_nEmptyParameters; ipar++) {
		printf("  p%d (%s): %f +/- %f\n", ipar, f_empty->GetParName(ipar), f_empty->GetParameter(ipar), f_empty->GetParError(ipar));
	}
	//printf("\n\n");
	//printf("excludeRegions.size() = %d\n", excludeRegions.size());
	return;
}
