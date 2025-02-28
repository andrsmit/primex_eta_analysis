#include "MggFitter.h"

void MggFitter::FitEmpty()
{
	InitializeEmptyFitFunction(&f_empty);
	
	emptyBinSize = h_empty->GetBinWidth(1);
	
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
			
			f_empty->SetParLimits(p0Par,  0.00, 1.e5);
			f_empty->SetParLimits(p1Par,  0.00, 1.00);
			f_empty->SetParLimits(p2Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p3Par, -1.e3, 1.e3);
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
			for(int ipar=0; ipar<4; ipar++) {
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
			
			f_empty->ReleaseParameter(p0Par);
			f_empty->ReleaseParameter(p1Par);
			f_empty->ReleaseParameter(p2Par);
			f_empty->ReleaseParameter(p3Par);
			
			f_empty->SetParLimits(p0Par,  0.00, 1.e5);
			f_empty->SetParLimits(p1Par,  0.00, 1.00);
			f_empty->SetParLimits(p2Par, -1.e3, 1.e3);
			f_empty->SetParLimits(p3Par, -1.e3, 1.e3);
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
	
	double minEtaFit  = 0.50;
	double maxEtaFit  = 0.60;
	
	double     NGuess = h_empty->Integral(h_empty->FindBin(minEtaFit), h_empty->FindBin(maxEtaFit));
	double    muGuess = 0.540;
	double sigmaGuess = 0.015;
	
	switch(emptyFitOption_eta) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			int     NPar = f_empty->GetParNumber("N_{#eta}");
			int    muPar = f_empty->GetParNumber("#mu_{#eta}");
			int sigmaPar = f_empty->GetParNumber("#sigma_{#eta}");
			
			f_empty->SetParameter(    NPar,     NGuess);
			f_empty->SetParameter(   muPar,    muGuess);
			f_empty->SetParameter(sigmaPar, sigmaGuess);
			
			f_empty->SetParLimits(    NPar, 0.00, 1.e5);
			f_empty->SetParLimits(   muPar, 0.52, 0.57);
			f_empty->SetParLimits(sigmaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_empty->GetParNumber("N_{#eta}");
			int dmuPar = f_empty->GetParNumber("#Delta#mu_{#eta}");
			
			f_empty->SetParameter(  NPar, NGuess);
			f_empty->SetParameter(dmuPar, 0.005);
			
			f_empty->SetParLimits(  NPar,  0.00, 1.e5);
			f_empty->SetParLimits(dmuPar, -0.05, 0.02);
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
			
			int     NPar = f_empty->GetParNumber("N_{#eta}");
			int    muPar = f_empty->GetParNumber("#mu_{#eta}");
			int sigmaPar = f_empty->GetParNumber("#sigma_{#eta}");
			
			f_empty->FixParameter(    NPar, f_empty->GetParameter(    NPar));
			f_empty->FixParameter(   muPar, f_empty->GetParameter(   muPar));
			f_empty->FixParameter(sigmaPar, f_empty->GetParameter(sigmaPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_empty->GetParNumber("N_{#eta}");
			int dmuPar = f_empty->GetParNumber("#Delta#mu_{#eta}");
			
			f_empty->FixParameter(  NPar, f_empty->GetParameter(  NPar));
			f_empty->FixParameter(dmuPar, f_empty->GetParameter(dmuPar));
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
			
			int     NPar = f_empty->GetParNumber("N_{#omega}");
			int    muPar = f_empty->GetParNumber("#mu_{#omega}");
			int sigmaPar = f_empty->GetParNumber("#sigma_{#omega}");
			int alphaPar = f_empty->GetParNumber("#alpha_{#omega}");
			int     nPar = f_empty->GetParNumber("n_{#omega}");
			
			f_empty->SetParameter(    NPar, 0.000);
			f_empty->SetParameter(   muPar, 0.760);
			f_empty->SetParameter(sigmaPar, 0.025);
			f_empty->SetParameter(alphaPar, 1.000);
			f_empty->SetParameter(    nPar, 2.000);
			
			f_empty->SetParLimits(    NPar, 0.000,  1.0e6);
			f_empty->SetParLimits(   muPar, 0.740,  0.780);
			f_empty->SetParLimits(sigmaPar, 0.015,  0.050);
			f_empty->SetParLimits(alphaPar, 0.200,  9.999);
			f_empty->SetParLimits(    nPar, 1.100, 49.999);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_empty->GetParNumber("N_{#omega}");
			int dmuPar = f_empty->GetParNumber("#Delta#mu_{#omega}");
			
			f_empty->SetParameter(  NPar,  0.00);
			f_empty->SetParameter(dmuPar, -0.01);
			
			f_empty->SetParLimits(  NPar,  0.00, 1.e5);
			f_empty->SetParLimits(dmuPar, -0.04, 0.03);
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
			
			int     NPar = f_empty->GetParNumber("N_{#omega}");
			int    muPar = f_empty->GetParNumber("#mu_{#omega}");
			int sigmaPar = f_empty->GetParNumber("#sigma_{#omega}");
			int alphaPar = f_empty->GetParNumber("#alpha_{#omega}");
			int     nPar = f_empty->GetParNumber("n_{#omega}");
			
			f_empty->FixParameter(    NPar, f_empty->GetParameter(    NPar));
			f_empty->FixParameter(   muPar, f_empty->GetParameter(   muPar));
			f_empty->FixParameter(sigmaPar, f_empty->GetParameter(sigmaPar));
			f_empty->FixParameter(alphaPar, f_empty->GetParameter(alphaPar));
			f_empty->FixParameter(    nPar, f_empty->GetParameter(    nPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_empty->GetParNumber("N_{#omega}");
			int dmuPar = f_empty->GetParNumber("#Delta#mu_{#omega}");
			
			f_empty->FixParameter(  NPar, f_empty->GetParameter(  NPar));
			f_empty->FixParameter(dmuPar, f_empty->GetParameter(dmuPar));
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
