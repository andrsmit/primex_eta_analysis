#include "EtaAnalyzer.h"
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
	}
	
	f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
	h_emptyWide->Fit(f_emptyWide, locFitOption);
	
	excludeRegions.clear();
	f_emptyWide->SetRange(minEmptyFitRange, maxEmptyFitRange);
	
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
				f_emptyWide->SetParLimits(locSigPar, 0.0075, 0.030);
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
			f_emptyWide->FixParameter(f_emptyWide->GetParNumber("#sigma_{fdc,5}"), 0.008873);
			
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
				
				f_emptyWide->ReleaseParameter(locNPar);
				f_emptyWide->ReleaseParameter(locMuPar);
				
				f_emptyWide->SetParLimits(locNPar,  0.0, 1.e5);
				f_emptyWide->SetParLimits(locMuPar, m_muFDC[i]-0.02, m_muFDC[i]+0.02);
				
				if(i>1) {
					f_emptyWide->ReleaseParameter(locSigPar);
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
