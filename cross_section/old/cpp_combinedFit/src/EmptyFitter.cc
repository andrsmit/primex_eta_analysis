#include "MggFitter.h"

void MggFitter::FitEmptyWide()
{
	InitializeEmptyFitFunction(&f_emptyWide);
	
	emptyBinSize = h_emptyWide->GetBinWidth(1);
	
	Option_t *locFitOption = "R0QL";
	
	InitializeEmptyWideFitParameters();
	excludeRegions.clear();
	
	excludeRegions.push_back({0.50,0.59});
	excludeRegions.push_back({0.70,0.85});
	GuessEmptyBkgdParameters();
	
	f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
	h_emptyWide->Fit(f_emptyWide, locFitOption);
	
	excludeRegions.clear();
	
	//------------------------------//
	// now try to fit peaking structures from FDC:
	
	if(emptyFitOption_fdc>0) {
		excludeRegions.clear();
		GuessEmptyFDCParameters();
		h_emptyWide->Fit(f_emptyWide, locFitOption);
	}
	
	//------------------------------//
	// now try to fit eta peak:
	
	if(emptyFitOption_eta>0) {
		excludeRegions.clear();
		
		GuessEmptyEtaParameters();
		f_emptyWide->SetRange(0.5,0.6);
		h_emptyWide->Fit(f_emptyWide, locFitOption);
		
		excludeRegions.push_back({0.70, 0.85});
		ReleaseEmptyBkgdParameters();
		f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
		h_emptyWide->Fit(f_emptyWide, locFitOption);
	}
	
	//------------------------------//
	// now try to fit omega peak:
	
	if(emptyFitOption_omega>0) {
		excludeRegions.clear();
		FixEmptyBkgdParameters();
		FixEmptyEtaParameters();
		FixEmptyFDCParameters();
		GuessEmptyOmegaParameters();
		f_emptyWide->SetRange(0.7,0.85);
		h_emptyWide->Fit(f_emptyWide, locFitOption);
		
		ReleaseEmptyBkgdParameters();
		f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
		h_emptyWide->Fit(f_emptyWide, locFitOption);
	}
	
	//------------------------------//
	ReleaseEmptyFDCParameters();
	ReleaseEmptyEtaParameters();
	
	f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
	h_emptyWide->Fit(f_emptyWide, locFitOption);
	
	excludeRegions.clear();
	f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
	
	/*
	// test
	f_emptyWide->SetParameter(1, 103.853507);
	f_emptyWide->SetParameter(3, 943.221777);
	f_emptyWide->SetParameter(4, 0.774695);
	f_emptyWide->SetParameter(5, 0.015000);
	f_emptyWide->SetParameter(6, 0.200000);
	f_emptyWide->SetParameter(7, 1.100000);
	f_emptyWide->SetParameter(8, -485442.086493);
	f_emptyWide->SetParameter(9, 816447.594345);
	f_emptyWide->SetParameter(10, -470457.637431);
	f_emptyWide->SetParameter(11, 166833.295377);
	f_emptyWide->SetParameter(12, -27107.413507);
	f_emptyWide->SetParameter(13, 89.571410);
	f_emptyWide->SetParameter(14, 0.615000);
	f_emptyWide->SetParameter(16, 112.291070);
	f_emptyWide->SetParameter(17, 0.531256);
	f_emptyWide->SetParameter(19, 326.685666);
	f_emptyWide->SetParameter(20, 0.451555);
	f_emptyWide->SetParameter(22, 106.037975);
	f_emptyWide->SetParameter(23, 0.405420);
	f_emptyWide->SetParameter(25, 82.103802);
	f_emptyWide->SetParameter(26, 0.331875);
	*/
	
	return;
}

void MggFitter::GuessEmptyBkgdParameters()
{
	if(emptyFitOption_bkgd==4) return;
	
	int p0Par = f_emptyWide->GetParNumber("p0,empty");
	int p1Par = f_emptyWide->GetParNumber("p1,empty");
	int p2Par = f_emptyWide->GetParNumber("p2,empty");
	int p3Par = f_emptyWide->GetParNumber("p3,empty");
	
	double p0Guess = 0., p1Guess = 0., p2Guess = 0., p3Guess = 0.;
	
	if(emptyFitOption_bkgd==2) {
		p0Guess =  (h_emptyWide->GetBinContent(h_emptyWide->FindBin(minEmptyFitRange)) - f_emptyWide->Eval(minEmptyFitRange))/emptyBinSize;
		p1Guess =  minEmptyFitRange;
		p2Guess =  0.0;
		p3Guess =  0.0;
	}
	
	switch(emptyFitOption_bkgd) {
		case 1:
		{
			// polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->SetParameter(locParIndex, 0.0);
				f_emptyWide->SetParLimits(locParIndex, -1.e7, 1.e7);
			}
			break;
		}
		case 2:
		{
			// exponential:
			
			f_emptyWide->SetParameter(p0Par, p0Guess);
			f_emptyWide->FixParameter(p1Par, p1Guess);
			f_emptyWide->SetParameter(p2Par, p2Guess);
			f_emptyWide->SetParameter(p3Par, p3Guess);
			
			f_emptyWide->SetParLimits(p0Par,  0.00, 1.e6);
			//f_emptyWide->SetParLimits(p1Par,  0.00, 1.00);
			f_emptyWide->SetParLimits(p2Par, -1.e3, 1.e3);
			f_emptyWide->SetParLimits(p3Par, -1.e3, 1.e3);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->SetParameter(locParIndex, 0.0001);
				f_emptyWide->SetParLimits(locParIndex, -1.e7, 1.e7);
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
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->FixParameter(locParIndex, f_emptyWide->GetParameter(locParIndex));
			}
			break;
		}
		case 2:
		{
			// Exponential:
			for(int ipar=0; ipar<4; ipar++) {
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->FixParameter(locParIndex, f_emptyWide->GetParameter(locParIndex));
			}
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->FixParameter(locParIndex, f_emptyWide->GetParameter(locParIndex));
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
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->ReleaseParameter(locParIndex);
				f_emptyWide->SetParLimits(locParIndex, -1.e7, 1.e7);
			}
			break;
		}
		case 2:
		{
			// Exponential:
			int p0Par = f_emptyWide->GetParNumber("p0,empty");
			int p1Par = f_emptyWide->GetParNumber("p1,empty");
			int p2Par = f_emptyWide->GetParNumber("p2,empty");
			int p3Par = f_emptyWide->GetParNumber("p3,empty");
			
			f_emptyWide->ReleaseParameter(p0Par);
			//f_emptyWide->ReleaseParameter(p1Par);
			f_emptyWide->ReleaseParameter(p2Par);
			f_emptyWide->ReleaseParameter(p3Par);
			
			f_emptyWide->SetParLimits(p0Par,  0.00, 1.e6);
			//f_emptyWide->SetParLimits(p1Par,  0.00, 1.00);
			f_emptyWide->SetParLimits(p2Par, -1.e3, 1.e3);
			f_emptyWide->SetParLimits(p3Par, -1.e3, 1.e3);
			break;
		}
		case 3:
		{
			// Chebyshev polynomial:
			for(int ipar=0; ipar<=emptyFitOption_poly; ipar++) {
				int locParIndex = f_emptyWide->GetParNumber(Form("p%d,empty",ipar));
				f_emptyWide->ReleaseParameter(locParIndex);
				f_emptyWide->SetParLimits(locParIndex, -1.e7, 1.e7);
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
	
	double     NGuess = h_emptyWide->Integral(h_emptyWide->FindBin(minEtaFit), h_emptyWide->FindBin(maxEtaFit));
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
			
			int     NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int    muPar = f_emptyWide->GetParNumber("#mu_{#eta,empty}");
			int sigmaPar = f_emptyWide->GetParNumber("#sigma_{#eta,empty}");
			
			f_emptyWide->SetParameter(    NPar,     NGuess);
			f_emptyWide->SetParameter(   muPar,    muGuess);
			f_emptyWide->SetParameter(sigmaPar, sigmaGuess);
			
			f_emptyWide->SetParLimits(    NPar, 0.00, 1.e5);
			f_emptyWide->SetParLimits(   muPar, 0.52, 0.57);
			f_emptyWide->SetParLimits(sigmaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int dmuPar = f_emptyWide->GetParNumber("#Delta#mu_{#eta,empty}");
			
			f_emptyWide->SetParameter(  NPar, NGuess);
			f_emptyWide->FixParameter(dmuPar, 0.0025);
			
			f_emptyWide->SetParLimits(  NPar,  0.000, 1.e5);
			//f_emptyWide->SetParLimits(dmuPar,  0.002, 0.01);
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
			
			int     NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int    muPar = f_emptyWide->GetParNumber("#mu_{#eta,empty}");
			int sigmaPar = f_emptyWide->GetParNumber("#sigma_{#eta,empty}");
			
			f_emptyWide->FixParameter(    NPar, f_emptyWide->GetParameter(    NPar));
			f_emptyWide->FixParameter(   muPar, f_emptyWide->GetParameter(   muPar));
			f_emptyWide->FixParameter(sigmaPar, f_emptyWide->GetParameter(sigmaPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int dmuPar = f_emptyWide->GetParNumber("#Delta#mu_{#eta,empty}");
			
			f_emptyWide->FixParameter(  NPar, f_emptyWide->GetParameter(  NPar));
			f_emptyWide->FixParameter(dmuPar, f_emptyWide->GetParameter(dmuPar));
			break;
		}
	}
	return;
}

void MggFitter::ReleaseEmptyEtaParameters() {
	
	switch(emptyFitOption_eta) {
		case 0:
		{
			break;
		}
		case 1:
		{
			// single Gaussian:
			
			int     NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int    muPar = f_emptyWide->GetParNumber("#mu_{#eta,empty}");
			int sigmaPar = f_emptyWide->GetParNumber("#sigma_{#eta,empty}");
			
			f_emptyWide->SetParameter(    NPar, f_emptyWide->GetParameter(    NPar));
			f_emptyWide->SetParameter(   muPar, f_emptyWide->GetParameter(   muPar));
			f_emptyWide->SetParameter(sigmaPar, f_emptyWide->GetParameter(sigmaPar));
			
			f_emptyWide->SetParLimits(    NPar, 0.00, 1.e5);
			f_emptyWide->SetParLimits(   muPar, 0.52, 0.57);
			f_emptyWide->SetParLimits(sigmaPar, 0.01, 0.03);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_emptyWide->GetParNumber("N_{#eta,empty}");
			int dmuPar = f_emptyWide->GetParNumber("#Delta#mu_{#eta,empty}");
			
			f_emptyWide->SetParameter(NPar, f_emptyWide->GetParameter(NPar));
			f_emptyWide->SetParLimits(NPar, 0.00, 1.e5);
			
			//f_emptyWide->SetParameter(dmuPar, f_emptyWide->GetParameter(dmuPar));
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
			
			int     NPar = f_emptyWide->GetParNumber("N_{#omega,empty}");
			int    muPar = f_emptyWide->GetParNumber("#mu_{#omega,empty}");
			int sigmaPar = f_emptyWide->GetParNumber("#sigma_{#omega,empty}");
			int alphaPar = f_emptyWide->GetParNumber("#alpha_{#omega,empty}");
			int     nPar = f_emptyWide->GetParNumber("n_{#omega,empty}");
			
			f_emptyWide->SetParameter(    NPar, 0.000);
			f_emptyWide->SetParameter(   muPar, 0.760);
			f_emptyWide->SetParameter(sigmaPar, 0.025);
			f_emptyWide->SetParameter(alphaPar, 1.000);
			f_emptyWide->SetParameter(    nPar, 2.000);
			
			f_emptyWide->SetParLimits(    NPar, 0.000,  1.0e6);
			f_emptyWide->SetParLimits(   muPar, 0.740,  0.780);
			f_emptyWide->SetParLimits(sigmaPar, 0.015,  0.050);
			f_emptyWide->SetParLimits(alphaPar, 0.200,  9.999);
			f_emptyWide->SetParLimits(    nPar, 1.100, 49.999);
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_emptyWide->GetParNumber("N_{#omega,empty}");
			int dmuPar = f_emptyWide->GetParNumber("#Delta#mu_{#omega,empty}");
			
			f_emptyWide->SetParameter(  NPar,  0.00);
			f_emptyWide->SetParameter(dmuPar, -0.01);
			
			f_emptyWide->SetParLimits(  NPar,  0.00, 1.e5);
			f_emptyWide->SetParLimits(dmuPar, -0.04, 0.03);
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
			
			int     NPar = f_emptyWide->GetParNumber("N_{#omega,empty}");
			int    muPar = f_emptyWide->GetParNumber("#mu_{#omega,empty}");
			int sigmaPar = f_emptyWide->GetParNumber("#sigma_{#omega,empty}");
			int alphaPar = f_emptyWide->GetParNumber("#alpha_{#omega,empty}");
			int     nPar = f_emptyWide->GetParNumber("n_{#omega,empty}");
			
			f_emptyWide->FixParameter(    NPar, f_emptyWide->GetParameter(    NPar));
			f_emptyWide->FixParameter(   muPar, f_emptyWide->GetParameter(   muPar));
			f_emptyWide->FixParameter(sigmaPar, f_emptyWide->GetParameter(sigmaPar));
			f_emptyWide->FixParameter(alphaPar, f_emptyWide->GetParameter(alphaPar));
			f_emptyWide->FixParameter(    nPar, f_emptyWide->GetParameter(    nPar));
			break;
		}
		case 2:
		{
			// Lineshape fit:
			
			int   NPar = f_emptyWide->GetParNumber("N_{#omega,empty}");
			int dmuPar = f_emptyWide->GetParNumber("#Delta#mu_{#omega,empty}");
			
			f_emptyWide->FixParameter(  NPar, f_emptyWide->GetParameter(  NPar));
			f_emptyWide->FixParameter(dmuPar, f_emptyWide->GetParameter(dmuPar));
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
				int locNPar   = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				int locMuPar  = f_emptyWide->GetParNumber(Form("#mu_{fdc,%d}",i+1));
				int locSigPar = f_emptyWide->GetParNumber(Form("#sigma_{fdc,%d}",i+1));
				
				// only let the parameters float if they fall within our fitting range:
				if(m_muFDC[i] < minEmptyFitRange) continue;
				
				f_emptyWide->SetParameter(locNPar, 0.0);
				f_emptyWide->SetParameter(locMuPar, m_muFDC[i]);
				f_emptyWide->SetParameter(locSigPar, 0.015);
				
				f_emptyWide->SetParLimits(locNPar,   0.0, 1.e5);
				f_emptyWide->SetParLimits(locMuPar,  m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				f_emptyWide->SetParLimits(locSigPar, 0.0075, 0.025);
			}
			
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#mu_{fdc,1}"),    0.619179);
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#sigma_{fdc,1}"), 0.027);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#mu_{fdc,2}"),    0.531688);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#sigma_{fdc,2}"), 0.027);
			
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#mu_{fdc,3}"),    0.453992);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#sigma_{fdc,3}"), 0.017433);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#mu_{fdc,4}"),    0.398570);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#sigma_{fdc,4}"), 0.008740);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#mu_{fdc,5}"),    0.338665);
			f_emptyWide->SetParameter(f_emptyWide->GetParNumber("#sigma_{fdc,5}"), 0.008873);
			
			break;
		}
		case 2:
		{
			int locNPar = f_emptyWide->GetParNumber("N_{fdc}");
			f_emptyWide->SetParameter(locNPar, 0.0);
			f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
			for(int i=0; i<m_muFDC.size(); i++) {
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]);
				f_emptyWide->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locNPar, 0.0);
				f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
				
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]);
				f_emptyWide->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
		case 4:
		{
			for(int i=0; i<m_muFDC_omega.size(); i++) {
				int locNPar = f_emptyWide->GetParNumber(Form("N_{fdc,#omega,%d}",i+1));
				f_emptyWide->SetParameter(locNPar, 0.0);
				f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
				
				int locMuPar = f_emptyWide->GetParNumber(Form("#mu_{fdc,#omega,%d}",i+1));
				f_emptyWide->SetParameter(locMuPar, m_muFDC_omega[i]);
				f_emptyWide->SetParLimits(locMuPar, m_muFDC_omega[i]-0.03, m_muFDC_omega[i]+0.03);
			}
			for(int i=0; i<m_muFDC_eta.size(); i++) {
				int locNPar = f_emptyWide->GetParNumber(Form("N_{fdc,#eta,%d}",i+1));
				f_emptyWide->SetParameter(locNPar, 0.0);
				f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
				
				int locMuPar = f_emptyWide->GetParNumber(Form("#mu_{fdc,#eta,%d}",i+1));
				f_emptyWide->SetParameter(locMuPar, m_muFDC_eta[i]);
				f_emptyWide->SetParLimits(locMuPar, m_muFDC_eta[i]-0.03, m_muFDC_eta[i]+0.03);
			}
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#eta,1}"), m_muFDC_eta[0]);
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#eta,2}"), m_muFDC_eta[1]);
			//f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#eta,3}"), m_muFDC_eta[2]);
			
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#omega,1}"), m_muFDC_omega[0]);
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#omega,2}"), m_muFDC_omega[1]);
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#mu_{fdc,#omega,3}"), m_muFDC_omega[2]);
			break;
		}
	}
	return;
}

void MggFitter::FixEmptyFDCParameters() {
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar   = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				int locMuPar  = f_emptyWide->GetParNumber(Form("#mu_{fdc,%d}",i+1));
				int locSigPar = f_emptyWide->GetParNumber(Form("#sigma_{fdc,%d}",i+1));
				
				f_emptyWide->FixParameter(locNPar,   f_emptyWide->GetParameter(locNPar));
				f_emptyWide->FixParameter(locMuPar,  f_emptyWide->GetParameter(locMuPar));
				f_emptyWide->FixParameter(locSigPar, f_emptyWide->GetParameter(locSigPar));
			}
			break;
		}
		case 2:
		{
			int locNPar = f_emptyWide->GetParNumber("N_{fdc}");
			f_emptyWide->FixParameter(locNPar, f_emptyWide->GetParameter(locNPar));
			for(int i=0; i<m_muFDC.size(); i++) {
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->FixParameter(locDeltaMuPar, f_emptyWide->GetParameter(locDeltaMuPar));
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				f_emptyWide->FixParameter(locNPar, f_emptyWide->GetParameter(locNPar));
				
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->FixParameter(locDeltaMuPar, f_emptyWide->GetParameter(locDeltaMuPar));
			}
			break;
		}
	}
	return;
}

void MggFitter::ReleaseEmptyFDCParameters() {
	
	switch(emptyFitOption_fdc) {
		case 0:
			break;
		case 1:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar   = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				int locMuPar  = f_emptyWide->GetParNumber(Form("#mu_{fdc,%d}",i+1));
				int locSigPar = f_emptyWide->GetParNumber(Form("#sigma_{fdc,%d}",i+1));
				
				f_emptyWide->SetParameter(locNPar,   f_emptyWide->GetParameter(locNPar));
				f_emptyWide->SetParameter(locMuPar,  f_emptyWide->GetParameter(locMuPar));
				
				f_emptyWide->SetParLimits(locNPar,   0.0, 1.e5);
				f_emptyWide->SetParLimits(locMuPar,  m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				
				if(i>1) {
					f_emptyWide->SetParameter(locSigPar, f_emptyWide->GetParameter(locSigPar));
					f_emptyWide->SetParLimits(locSigPar, 0.0075, 0.025);
				}
			}
			break;
		}
		case 2:
		{
			int locNPar = f_emptyWide->GetParNumber("N_{fdc}");
			f_emptyWide->SetParameter(locNPar, f_emptyWide->GetParameter(locNPar));
			f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
			for(int i=0; i<m_muFDC.size(); i++) {
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locDeltaMuPar, f_emptyWide->GetParameter(locDeltaMuPar));
				f_emptyWide->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
		case 3:
		{
			for(int i=0; i<m_muFDC.size(); i++) {
				int locNPar = f_emptyWide->GetParNumber(Form("N_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locNPar, f_emptyWide->GetParameter(locNPar));
				f_emptyWide->SetParLimits(locNPar, 0.0, 1.e5);
				
				int locDeltaMuPar = f_emptyWide->GetParNumber(Form("#Delta#mu_{fdc,%d}",i+1));
				f_emptyWide->SetParameter(locDeltaMuPar, f_emptyWide->GetParameter(locDeltaMuPar));
				f_emptyWide->SetParLimits(locDeltaMuPar, m_muFDC[i]-m_muFDC[0]-0.02, m_muFDC[i]-m_muFDC[0]+0.02);
			}
			break;
		}
	}
	return;
}

void MggFitter::DumpEmptyFitParameters()
{
	if(f_emptyWide==nullptr) {
		printf("\n\nNo empty target fit was performed.\n");
		return;
	}
	printf("\n\nEmpty-Target Fit parameters:\n");
	for(int ipar=0; ipar<f_emptyWide->GetNpar(); ipar++) {
		printf("  p%d (%s): %f +/- %f\n", ipar, f_emptyWide->GetParName(ipar), f_emptyWide->GetParameter(ipar), f_emptyWide->GetParError(ipar));
	}
	//printf("\n\n");
	//printf("excludeRegions.size() = %d\n", excludeRegions.size());
	return;
}
