#include "EtaAnalyzer.h"
#include "CrossSection.h"

//---------------------------------------------------------//
// Setter functions for run-specific member variables:

void EtaAnalyzer::SetPhase(int phase)
{
	if((phase<1) || (phase>3)) {
		cout << "\nUnsupported PrimEx-eta phase number provided.\n" << endl;
		exit(1);
	}
	m_phase = phase;
	return;
}

void EtaAnalyzer::SetAnalysisOption(int option)
{
	if((option<0) || (option>8)) {
		cout << "\nUnsupported analysis option provided.\n" << endl;
		exit(1);
	}
	m_analysisOption = option;
	return;
}

void EtaAnalyzer::SetVetoOption(int option)
{
	if((option<0) || (option>8)) {
		cout << "\nUnsupported veto option provided.\n" << endl;
		exit(1);
	}
	m_vetoOption = option;
	return;
}

void EtaAnalyzer::SetMggHistName(TString name)
{
	m_mggHistName = name;
	return;
}
void EtaAnalyzer::SetMatrixHistName(TString name)
{
	m_matrixHistName = name;
	return;
}

//---------------------------------------------------------//
// Setter functions for binning-related member variables:

void EtaAnalyzer::SetRebinsMgg(int rebins) 
{
	m_rebinsMgg  = rebins;
	m_mggBinSize = (1.e-3) * (double)m_rebinsMgg;
	return;
}
void EtaAnalyzer::SetRebinsTheta(int rebins) 
{
	m_rebinsTheta  = rebins;
	m_reconAngleBinSize = 0.01 * (double)m_rebinsTheta;
	return;
}
void EtaAnalyzer::SetRebinsEmptyMgg(int rebins) 
{
	m_rebinsEmptyMgg = rebins;
	m_emptyMggBinSize = (1.e-3) * (double)m_rebinsEmptyMgg;
	return;
}
void EtaAnalyzer::SetBeamEnergy(double min, double max)
{
	m_minBeamEnergy = min;
	m_maxBeamEnergy = max;
	return;
}
void EtaAnalyzer::SetBeamEnergy()
{
	if(m_phase==1) {
		m_minBeamEnergy =  8.0;
		m_maxBeamEnergy = 10.9;
	}
	else if(m_phase==2) {
		m_minBeamEnergy =  8.0;
		m_maxBeamEnergy = 10.0;
	}
	else if(m_phase==3) {
		m_minBeamEnergy =  8.0;
		m_maxBeamEnergy = 11.3;
	}
	return;
}

//---------------------------------------------------------//
// Getter functions for binning-related member variables:

void EtaAnalyzer::GetBeamEnergyBinning(double &binSize, double &min, double &max)
{
	binSize = m_beamEnergyBinSize;
	min     = m_minBeamEnergy;
	max     = m_maxBeamEnergy;
	return;
}
void EtaAnalyzer::GetReconAngleBinning(double &binSize, double &min, double &max)
{
	binSize = m_reconAngleBinSize;
	min     = m_minReconAngle;
	max     = m_maxReconAngle;
	return;
}
void EtaAnalyzer::GetThrownAngleBinning(double &binSize, double &min, double &max)
{
	binSize = m_thrownAngleBinSize;
	min     = m_minThrownAngle;
	max     = m_maxThrownAngle;
	return;
}

//---------------------------------------------------------//
// Setter functions for fitting-related member variables:

void EtaAnalyzer::SetSubtractEmptyTarget(int option)
{
	if((option!=0) && (option!=1)) {
		printf("\nUnsupported empty target subtraction option provided (should be 0 or 1).");
		printf(" Empty target will be subtracted prior to fit.\n");
	}
	else {
		m_subtractEmpty = option;
	}
	return;
}

void EtaAnalyzer::SetFitEmptyTarget(int option)
{
	if((option!=0) && (option!=1)) {
		printf("\nUnsupported empty target fit option provided (should be 0 or 1).");
		printf(" Empty target will not be fit to fitting full target data.\n");
	}
	else {
		m_fitOption_empty = option;
	}
	return;
}

//---------------//

void EtaAnalyzer::SetFitOption_signal(int option) 
{
	if((option<1) || (option>2)) {
		printf("\nUnsupported signal fit option provided (should be 1-2). Using default option: %d\n", m_fitOption_signal);
	}
	else {
		m_fitOption_signal = option;
	}
	return;
}
void EtaAnalyzer::SetFitOption_bkgd(int option) 
{
	if((option<0) || (option>4)) {
		printf("\nUnsupported background fit option provided (should be 1-5). Using default option: %d\n", m_fitOption_bkgd);
	}
	else {
		m_fitOption_bkgd = option;
	}
	return;
}
void EtaAnalyzer::SetFitOption_bkgd(int option, int order) 
{
	SetFitOption_bkgd(option);
	
	int locOrder = order;
	if(order<0)      locOrder = 0;
	else if(order>5) locOrder = 5;
	m_fitOption_poly = locOrder;
	
	return;
}
void EtaAnalyzer::SetFitOption_omega(int option)
{
	if((option<0) || (option>5)) {
		printf("\nUnsupported omega fit option provided (should be 1-3). Using default option: %d\n", m_fitOption_omega);
	}
	else {
		m_fitOption_omega = option;
	}
	return;
}
void EtaAnalyzer::SetFitOption_etap(int option)
{
	if((option<0) || (option>1)) {
		printf("\nUnsupported eta-prime fit option provided (should be 0-1). Using default option: %d\n", m_fitOption_etap);
	}
	else {
		m_fitOption_etap = option;
	}
	return;
}

//---------------//

void EtaAnalyzer::SetEmptyFitOption_eta(int option) 
{
	if((option<0) || (option>2)) {
		printf("\nUnsupported gas eta fit option provided (should be 0-2). Using default option: %d\n", m_emptyFitOption_eta);
	}
	else {
		m_emptyFitOption_eta = option;
	}
	return;
}
void EtaAnalyzer::SetEmptyFitOption_omega(int option) 
{
	if((option<0) || (option>2)) {
		printf("\nUnsupported gas omega fit option provided (should be 0-2). Using default option: %d\n", m_emptyFitOption_omega);
	}
	else {
		m_emptyFitOption_omega = option;
	}
	return;
}
void EtaAnalyzer::SetEmptyFitOption_fdc(int option) 
{
	if((option<0) || (option>4)) {
		printf("\nUnsupported fdc fit option provided (should be 0-3). Using default option: %d\n", m_emptyFitOption_fdc);
	}
	else {
		m_emptyFitOption_fdc = option;
	}
	return;
}
void EtaAnalyzer::SetEmptyFitOption_bkgd(int option, int order) 
{
	SetEmptyFitOption_bkgd(option);
	
	int locOrder = order;
	if(order<0)      locOrder = 0;
	else if(order>5) locOrder = 5;
	m_emptyFitOption_poly = locOrder;
	
	return;
}
void EtaAnalyzer::SetEmptyFitOption_bkgd(int option)
{
	if((option<1) || (option>3)) {
		printf("\nUnsupported empty bkgd fit option provided (should be 1-3). Using default option: %d\n", m_emptyFitOption_bkgd);
	}
	else {
		m_emptyFitOption_bkgd = option;
	}
	return;
}

//---------------//

void EtaAnalyzer::SetFitRange(double min, double max)
{
	double locMin = 0.30;
	//if(m_phase==1) locMin = 0.30;
	//else           locMin = 0.40;
	
	m_minFitRange = min < locMin ? locMin : min;
	m_maxFitRange = max;
	
	if((m_maxFitRange>0.72) && (m_fitOption_omega==0)) m_fitOption_omega = 1;
	return;
}

void EtaAnalyzer::SetEmptyFitRange(double min, double max)
{
	double locMin = 0.30;
	//if(m_phase==1) locMin = 0.30;
	//else           locMin = 0.40;
	
	m_minEmptyFitRange = min < locMin ? locMin : min;
	m_maxEmptyFitRange = max;
	return;
}

//---------------------------------------------------------//

int EtaAnalyzer::GetFitOption(int opt)
{
	switch(opt) {
		case 1:
			return m_fitOption_signal;
		case 2:
			return m_fitOption_bkgd;
		case 3:
			return m_fitOption_poly;
		case 4:
			return m_fitOption_omega;
		case 5:
			return m_fitOption_etap;
		default:
			return 0;
	}
}

int EtaAnalyzer::GetEmptyFitOption(int opt)
{
	switch(opt) {
		case 0:
			return m_fitOption_empty;
		case 1:
			return m_emptyFitOption_eta;
		case 2:
			return m_emptyFitOption_omega;
		case 3:
			return m_emptyFitOption_fdc;
		case 4:
			return m_emptyFitOption_bkgd;
		case 5:
			return m_emptyFitOption_poly;
		default:
			return 0;
	}
}

void EtaAnalyzer::GetFitRange(double &min, double &max)
{
	min = m_minFitRange;
	max = m_maxFitRange;
	return;
}

void EtaAnalyzer::GetEmptyFitRange(double &min, double &max)
{
	min = m_minEmptyFitRange;
	max = m_maxEmptyFitRange;
	return;
}

//---------------------------------------------------------//

TString EtaAnalyzer::GetFitOptionStr(int option)
{
	TString optString = "";
	switch(option) {
		case 0:
			switch(m_fitOption_signal) {
				case 1:
					optString = "histogram as lineshape";
					break;
				case 2:
					optString = "function-parameterized lineshape";
					break;
			}
			break;
		case 1:
			switch(m_fitOption_bkgd) {
				case 0:
					optString = "no background";
					break;
				case 1:
					optString = Form("polynomial (order %d)", m_fitOption_poly);
					break;
				case 2:
					optString = "exponential";
					break;
				case 3:
					optString = Form("Chebyshev polynomial (order %d)", m_fitOption_poly);
					break;
			}
			break;
		case 2:
			switch(m_fitOption_omega) {
				case 0:
					optString = "No fit to omega peak";
					break;
				case 1:
					optString = "Crystal Ball (floating parameters)";
					break;
				case 2:
					optString = "Simulated Lineshape (with double Crystal Ball parameterization)";
					break;
				case 3:
					optString = "Simulated Lineshape";
					break;
				case 4:
					optString = "Simulated Lineshape for omega and rho";
					break;
			}
			break;
		case 3:
			switch(m_fitOption_etap) {
				case 0:
					optString = "No fit attempted";
					break;
				case 1:
					optString = "Gaussian";
					break;
			}
			break;
	}
	return optString;
}

TString EtaAnalyzer::GetBkgdFitName() 
{
	switch(m_fitOption_bkgd) {
		case 0:
			return "";
		case 1:
			switch(m_fitOption_poly) {
				case 1:
					return "1st order poly";
				case 2:
					return "2nd order poly";
				case 3:
					return "3rd order poly";
				default:
					return Form("%dth order poly", m_fitOption_poly);
			}
			break;
		case 2:
			return "exponential";
		case 3:
			switch(m_fitOption_poly) {
				case 1:
					return "1st order poly";
				case 2:
					return "2nd order poly";
				case 3:
					return "3rd order poly";
				default:
					return Form("%dth order poly", m_fitOption_poly);
			}
			break;
	}
	TString locString = "";
	return locString;
}

TString EtaAnalyzer::GetEmptyFitOptionStr()
{
	TString locString = m_subtractEmpty ? "subtracted before fit" : "Fit performed";
	return locString;
}

//---------------------------------------------------------//

void EtaAnalyzer::InitializeFitCanvas()
{
	cFit = new TCanvas("cFit", "Mgg Fit", 1200, 900);
	cFit->cd();
	
	pFit = new TPad("padFit", "Mgg Fit", 0.005, 0.3025, 0.995, 0.9950);
	pFit->SetLeftMargin(0.10);
	pFit->SetRightMargin(0.02);
	pFit->SetTopMargin(0.075);
	pFit->SetBottomMargin(0.015);
	pFit->SetTickx(); pFit->SetTicky();
	pFit->SetFrameLineWidth(2);
	pFit->SetGrid();
	
	pRes = new TPad("padRes", "Mgg Res", 0.005, 0.0050, 0.995, 0.2975);
	pRes->SetLeftMargin(0.10);
	pRes->SetRightMargin(0.02);
	pRes->SetTopMargin(0.005);
	pRes->SetBottomMargin(0.325);
	pRes->SetTickx(); pRes->SetTicky();
	pRes->SetFrameLineWidth(2);
	pRes->SetGrid();
	pFit->Draw();
	pRes->Draw();
	return;
}

void EtaAnalyzer::InitializeEmptyCanvas()
{
	cEmpty = new TCanvas("cEmpty", "Empty Mgg Fit", 1200, 650);
	styleCanvas(cEmpty);
	cEmpty->SetLeftMargin(0.10);
	cEmpty->SetRightMargin(0.02);
	cEmpty->SetTopMargin(0.075);
	cEmpty->SetBottomMargin(0.13);
	cEmpty->SetGrid();
	return;
}

void EtaAnalyzer::DrawInvariantMass(double minAngle, double maxAngle)
{
	if(cFit==NULL) InitializeFitCanvas();
	
	int minAngleBin = h_mggVsThetaFull->GetXaxis()->FindBin(minAngle);
	int maxAngleBin = h_mggVsThetaEmpty->GetXaxis()->FindBin(maxAngle)-1;
	
	TH1F *locHistFull  = (TH1F*)h_mggVsThetaFull->ProjectionY("locHistFull", minAngleBin, maxAngleBin);
	TH1F *locHistEmpty = (TH1F*)h_mggVsThetaEmpty->ProjectionY("locHistEmpty", minAngleBin, maxAngleBin);
	
	locHistFull->Rebin(m_rebinsMgg);
	locHistEmpty->Rebin(m_rebinsMgg);
	
	styleMggHistogram(locHistFull);
	styleMggHistogram(locHistEmpty, kBlue);
	
	cFit->cd();
	locHistFull->Draw("PE1");
	locHistEmpty->Draw("PE1 same");
	cFit->Update();
	cFit->Modified();
	return;
}

void EtaAnalyzer::DrawFitResult(TH1F *h1, TH1F *hPull, TF1 *fFit, TF1 *fSignal, TF1 *fBkgd, TF1 *fOmega, TF1 *fRho,
	TF1 *fHadronicBkgd, TF1 *fEtaPi, TF1 *fEmpty, double minAngle, double maxAngle)
{
	if(cFit==NULL) InitializeFitCanvas();
	
	if(l0==NULL) {
		l0 = new TLine(m_minFitRange, 0.0, m_maxFitRange, 0.0);
		l0->SetLineColor(kBlack);
	}
	if(lm==NULL) {
		lm = new TLine(m_minFitRange, -2.0, m_maxFitRange, -2.0);
		lm->SetLineColor(kRed);
	}
	if(lp==NULL) {
		lp = new TLine(m_minFitRange, +2.0, m_maxFitRange, +2.0);
		lp->SetLineColor(kRed);
	}
	if(lx1==NULL) {
		lx1 = new TLine(0.5, 0.0, 0.5, 1.0);
		lx1->SetLineColor(kRed);
		lx1->SetLineStyle(4);
	}
	if(lx2==NULL) {
		lx2 = new TLine(0.6, 0.0, 0.6, 1.0);
		lx2->SetLineColor(kRed);
		lx2->SetLineStyle(4);
	}
	
	styleMggHistogram(h1);
	
	//h1->SetMinimum(-100.);
	h1->GetXaxis()->SetRangeUser(m_minFitRange,m_maxFitRange);
	hPull->GetXaxis()->SetRangeUser(m_minFitRange, m_maxFitRange);
	
	//h1->SetMinimum(-20.);
	//h1->GetXaxis()->SetRangeUser(0.45,0.65);
	//hPull->GetXaxis()->SetRangeUser(0.45,0.65);
	
	//gStyle->SetOptFit(0);
	
	pFit->cd();
	h1->Draw("PE1");
	fFit->Draw("same");
	fSignal->Draw("same");
	fBkgd->Draw("same");
	fOmega->Draw("same");
	if(fRho) fRho->Draw("same");
	if(fHadronicBkgd) fHadronicBkgd->Draw("same");
	if(fEtaPi) fEtaPi->Draw("same");
	if(fEmpty) fEmpty->Draw("same");
	
	pFit->Update();
	
	l0->Draw("same");
	lx1->SetY1(gPad->GetUymin());
	lx1->SetY2(gPad->GetUymax());
	lx1->Draw("same");
	lx2->SetY1(gPad->GetUymin());
	lx2->SetY2(gPad->GetUymax());
	lx2->Draw("same");
	
	TLegend *locLeg = new TLegend(0.655, 0.550, 0.950, 0.900);
	locLeg->AddEntry(fFit,    "Full Fit",         "l");
	locLeg->AddEntry(fSignal, "Signal Lineshape", "l");
	if(fEmpty) locLeg->AddEntry(fEmpty, "Empty-Target Bkgd", "l");
	if(fEtaPi) locLeg->AddEntry(fEtaPi, "Eta+Pion Bkgd", "l");
	if(fHadronicBkgd) locLeg->AddEntry(fHadronicBkgd, "Other Hadronic Bkgd", "l");
	locLeg->AddEntry(fBkgd, Form("Additional Bkgd (%s)", GetBkgdFitName().Data()), "l");
	if(fRho) locLeg->AddEntry(fRho, "#rho bkgd", "l");
	//locLeg->Draw();
	
	TLatex locLatex;
	locLatex.DrawLatexNDC(0.137, 0.834, 
		Form("#scale[1.0]{#theta_{#gamma#gamma}: %.2f#circ - %.2f#circ}", 
		minAngle, maxAngle));
	
	pRes->cd();
	
	hPull->Draw("PE1");
	l0->Draw("same");
	lm->Draw("same");
	lp->Draw("same");
	
	pRes->Update();
	cFit->Update();
	cFit->Modified();
	//cFit->SaveAs(Form("fit_result_%.2fdeg_%.2fdeg.pdf",minAngle,maxAngle));
	
	return;
}

void EtaAnalyzer::ExtractAngularYield(MggFitter& locFitter, int drawOption)
{
	if(m_binningSet==false) InitializeBinning();
	
	for(int iThetaBin=0; iThetaBin<m_angularBin.size(); iThetaBin++) {
		
		double locAngle    = m_angularBin[iThetaBin].first;
		double locMinAngle = locAngle - m_angularBin[iThetaBin].second;
		double locMaxAngle = locAngle + m_angularBin[iThetaBin].second;
		
		//if(locAngle<0.5 || locAngle>0.38) continue;
		//if(locAngle<4.0) continue;
		
		//----------------------------------------------//
		// Get 1-d projection of invariant mass spectrum:
		
		int minAngleBin = m_rebinsTheta*(iThetaBin) + 1;
		int maxAngleBin = m_rebinsTheta*(iThetaBin+1);
		
		// verify that minAngleBin and maxAngleBin correctly correspond to bin edges:
		
		double locMinAngleHist = h_mggVsThetaFull->GetXaxis()->GetBinCenter(minAngleBin) 
			- 0.5*h_mggVsThetaFull->GetXaxis()->GetBinWidth(1);
		double locMaxAngleHist = h_mggVsThetaFull->GetXaxis()->GetBinCenter(maxAngleBin) 
			+ 0.5*h_mggVsThetaFull->GetXaxis()->GetBinWidth(1);
		if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
			printf("\nWarning: Histogram bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
			printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
			printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
		}
		
		TH1F *locHistFull  = (TH1F*)h_mggVsThetaFull->ProjectionY("locHistFull", minAngleBin, maxAngleBin, "e");
		TH1F *locHistEmpty = (TH1F*)h_mggVsThetaEmpty->ProjectionY("locHistEmpty", minAngleBin, maxAngleBin, "e");
		
		TH1F *locHistFull_acc  = (TH1F*)h_mggVsThetaFull_acc->ProjectionY("locHistFull_acc", minAngleBin, maxAngleBin, "e");
		TH1F *locHistEmpty_acc = (TH1F*)h_mggVsThetaEmpty_acc->ProjectionY("locHistEmpty_acc", minAngleBin, maxAngleBin, "e");
		
		if(m_subtractEmpty) {
			locHistFull->Add(locHistEmpty,-1.0);
			locHistFull_acc->Add(locHistEmpty_acc,-1.0);
		}
		
		locHistFull->Rebin(m_rebinsMgg);
		locHistFull_acc->Rebin(m_rebinsMgg);
		
		// set bin errors for empty bins to 1, otherwise they are ignored in a chi-squared minimization:
		for(int ibin=1; ibin<=locHistFull->GetXaxis()->GetNbins(); ibin++) {
			double locBinError = locHistFull->GetBinError(ibin);
			if((locBinError>0.0) && (locBinError<1.0)) {
				locHistFull->SetBinError(ibin, 1.0);
			}
		}
		
		//----------------------------------------------//
		// Set up fitter object:
		
		InitializeMggFitter(locFitter, this, m_angularBin[iThetaBin].first, m_angularBin[iThetaBin].second);
		locFitter.SetData(locHistFull, locHistFull_acc);
		
		//=============================//
		// Offset between simulated lineshape and data:
		// 8/26/25: applied manual corrections in the analysis code, so this should now be fixed to 0:
		
		locFitter.lineshapeOffset = 0.0;
		/*
		double p0 =  0.00165029;
		double p1 = -0.00429059;
		double p2 =  0.00403517;
		double p3 = -0.00121199;
		double p4 = 0.000119845;
		locFitter.lineshapeOffset = p0 + p1*locAngle + p2*pow(locAngle,2.0) 
			+ p3*pow(locAngle,3.0) + p4*pow(locAngle,4.0);
		*/
		//=============================//
		
		double locBinSize = 0.0;
		double lineshapeWindowSize, lineshapeAngleLow, lineshapeAngleHigh;
		double binSizeRatio = 1.0;
		
		//----------------------------------------------//
		// Signal MC Lineshape
		
		// fraction of incoherent cross section in this angular bin:
		locFitter.incFraction_theory = h_incFraction->GetBinContent(h_incFraction->FindBin(locAngle+1.e-6));
		
		// Coherent MC:
		
		locBinSize  = h_etaLineshapeCoh->GetXaxis()->GetBinWidth(1);
		minAngleBin = h_etaLineshapeCoh->GetXaxis()->FindBin(locMinAngle + 0.5*locBinSize);
		maxAngleBin = h_etaLineshapeCoh->GetXaxis()->FindBin(locMaxAngle - 0.5*locBinSize)-1;
		
		locMinAngleHist = h_etaLineshapeCoh->GetXaxis()->GetBinCenter(minAngleBin) 
			- 0.5*h_etaLineshapeCoh->GetXaxis()->GetBinWidth(1);
		locMaxAngleHist = h_etaLineshapeCoh->GetXaxis()->GetBinCenter(maxAngleBin) 
			+ 0.5*h_etaLineshapeCoh->GetXaxis()->GetBinWidth(1);
		if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
			printf("\nWarning: Eta lineshape bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
			printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
			printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
		}
		
		TH1F *hEtaCoh = (TH1F*)h_etaLineshapeCoh->ProjectionY("hEtaCoh", minAngleBin, maxAngleBin, "e");
		locFitter.SetCohLineshape(hEtaCoh);
		
		// Quasifree MC:
		
		/*TH1F *hEtaQF  = (TH1F*)h_etaLineshapeQF->ProjectionY("hEtaQF", minAngleBin, maxAngleBin, "e");
		locFitter.SetQFLineshape(hEtaQF, 1);*/
		
		lineshapeWindowSize = 0.1;
		
		lineshapeAngleLow  = locAngle - lineshapeWindowSize;
		lineshapeAngleHigh = locAngle + lineshapeWindowSize;
		if(lineshapeAngleLow < 0.0) {
			lineshapeAngleLow  = 0.00;
			lineshapeAngleHigh = 2.0*lineshapeWindowSize;
		}
		TH1F *hEtaQF = (TH1F*)h_etaLineshapeQF->ProjectionY("hEtaQF",
			h_etaLineshapeQF->GetXaxis()->FindBin(lineshapeAngleLow),
			h_etaLineshapeQF->GetXaxis()->FindBin(lineshapeAngleHigh)-1, "e");
		hEtaQF->Scale((locMaxAngle-locMinAngle)/(2.0*lineshapeWindowSize));
		locFitter.SetQFLineshape(hEtaQF);
		
		//----------------------------------------------//
		// Omega Lineshape
		
		lineshapeWindowSize = 0.5;
		
		lineshapeAngleLow  = locAngle - lineshapeWindowSize;
		lineshapeAngleHigh = locAngle + lineshapeWindowSize;
		if(lineshapeAngleLow < 0.0) {
			lineshapeAngleLow  = 0.00;
			lineshapeAngleHigh = 2.0*lineshapeWindowSize;
		}
		TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega",
			h_omegaLineshape->GetXaxis()->FindBin(lineshapeAngleLow),
			h_omegaLineshape->GetXaxis()->FindBin(lineshapeAngleHigh)-1, "e");
		hOmega->Rebin(m_rebinsMgg);
		hOmega->Scale((locMaxAngle-locMinAngle)/(2.0*lineshapeWindowSize));
		locFitter.SetOmegaLineshape(hOmega);
		
		if(m_fitOption_omega==4) {
			lineshapeWindowSize = 2.0;
			
			lineshapeAngleLow  = locAngle - lineshapeWindowSize;
			lineshapeAngleHigh = locAngle + lineshapeWindowSize;
			if(lineshapeAngleLow < 0.0) {
				lineshapeAngleLow  = 0.00;
				lineshapeAngleHigh = 2.0*lineshapeWindowSize;
			}
			TH1F *hRho = (TH1F*)h_rhoLineshape->ProjectionY("hRho",
				h_rhoLineshape->GetXaxis()->FindBin(lineshapeAngleLow),
				h_rhoLineshape->GetXaxis()->FindBin(lineshapeAngleHigh)-1, "e");
			hRho->Rebin(m_rebinsMgg);
			hRho->Scale((locMaxAngle-locMinAngle)/(2.0*lineshapeWindowSize));
			locFitter.SetRhoLineshape(hRho);
		}
		
		//----------------------------------------------//
		// Hadronic Background Lineshapes:
		
		lineshapeWindowSize = 2.0;
		
		lineshapeAngleLow  = locAngle - lineshapeWindowSize;
		lineshapeAngleHigh = locAngle + lineshapeWindowSize;
		if(lineshapeAngleLow < 0.0) {
			lineshapeAngleLow  = 0.00;
			lineshapeAngleHigh = 2.0*lineshapeWindowSize;
		}
		TH1F *hHadronicBkgd = (TH1F*)h_hadronicBkgdLineshape->ProjectionY("hHadronicBkgd",
			h_hadronicBkgdLineshape->GetXaxis()->FindBin(lineshapeAngleLow),
			h_hadronicBkgdLineshape->GetXaxis()->FindBin(lineshapeAngleHigh)-1, "e");
		hHadronicBkgd->Rebin(m_rebinsMgg);
		
		// assume flat distribution in theta and normalize by bin size ratio:
		
		TH1F *lochHadronicBkgd = (TH1F*)h_hadronicBkgdLineshape->ProjectionY("hNarrowHadronicBkgd",
			minAngleBin, maxAngleBin, "e");
		double nHadronicNarrow = lochHadronicBkgd->Integral();
		double nHadronicWide   = hHadronicBkgd->Integral();
		if(nHadronicWide>=1.0) {
			hHadronicBkgd->Scale(nHadronicNarrow/nHadronicWide);
		}
		locFitter.SetHadronicBkgdLineshape(hHadronicBkgd);
		
		double locHadronicBkgdFrac    = h_HadronicBkgdFraction_bggen->GetBinContent(h_HadronicBkgdFraction_bggen->FindBin(locAngle));
		double locHadronicBkgdFracErr = h_HadronicBkgdFraction_bggen->GetBinError(h_HadronicBkgdFraction_bggen->FindBin(locAngle));
		locFitter.SetHadronicBkgdFraction(locHadronicBkgdFrac, locHadronicBkgdFracErr);
		
		TH1F *hEtaPionBkgd = (TH1F*)h_eta1PionLineshape->ProjectionY("hEtaPionBkgd",
			h_eta1PionLineshape->GetXaxis()->FindBin(lineshapeAngleLow),
			h_eta1PionLineshape->GetXaxis()->FindBin(lineshapeAngleHigh)-1, "e");
		hEtaPionBkgd->Rebin(m_rebinsMgg);
		
		TH1F *lochEtaPionBkgd = (TH1F*)h_eta1PionLineshape->ProjectionY("hNarrowEtaPionBkgd",
			minAngleBin, maxAngleBin, "e");
		double nEtaPiNarrow = lochEtaPionBkgd->Integral();
		double nEtaPiWide   = hEtaPionBkgd->Integral();
		if(nEtaPiWide>0.0) hEtaPionBkgd->Scale(nEtaPiNarrow/nEtaPiWide);
		
		locFitter.SetEtaPionLineshape(hEtaPionBkgd);
		
		double locEtaPionFrac    = h_EtaPionBkgdFraction_bggen->GetBinContent(h_EtaPionBkgdFraction_bggen->FindBin(locAngle));
		double locEtaPionFracErr = h_EtaPionBkgdFraction_bggen->GetBinError(h_EtaPionBkgdFraction_bggen->FindBin(locAngle));
		locFitter.SetEtaPionBkgdFraction(locEtaPionFrac, locEtaPionFracErr);
		
		//----------------------------------------------//
		// Do fit and extract yield:
		
		TH1F *locHistEmptyWide, *locHistEmptyWide_acc;
		double emptyAngleLow  = locMinAngle;
		double emptyAngleHigh = locMaxAngle;
		
		if(m_fitOption_empty==1) {
			
			// To get the pdf of the empty target background, we need to combine a wider angular range.
			
			double emptyWindowSize = m_phase==1 ? 0.225 : 0.125;
			
			if((locMaxAngle-locMinAngle) > (2.0*emptyWindowSize)) emptyWindowSize = 0.5*(locMaxAngle-locMinAngle);
			
			emptyAngleLow  = locAngle - emptyWindowSize;
			emptyAngleHigh = locAngle + emptyWindowSize;
			if(emptyAngleLow < 0.0) {
				emptyAngleLow  = 0.00;
				emptyAngleHigh = 2.0*emptyWindowSize;
			}
			
			locHistEmptyWide     = (TH1F*)h_mggVsThetaEmpty->ProjectionY("EmptyHistWide",
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleLow),
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleHigh)-1, "e");
			locHistEmptyWide_acc = (TH1F*)h_mggVsThetaEmpty_acc->ProjectionY("EmptyHistWide_acc",
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleLow),
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleHigh)-1, "e");
			
			locHistEmptyWide->Rebin(m_rebinsEmptyMgg);
			locHistEmpty->Rebin(m_rebinsEmptyMgg);
			locHistEmptyWide_acc->Rebin(m_rebinsEmptyMgg);
			locHistEmpty_acc->Rebin(m_rebinsEmptyMgg);
			
			double nEmptyNarrow  = locHistEmpty->Integral(
				locHistEmpty->GetXaxis()->FindBin(m_minEmptyFitRange), locHistEmpty->GetXaxis()->FindBin(m_maxEmptyFitRange));
			double nEmptyWide    = locHistEmptyWide->Integral(
				locHistEmptyWide->GetXaxis()->FindBin(m_minEmptyFitRange), locHistEmptyWide->GetXaxis()->FindBin(m_maxEmptyFitRange));
			double locEmptyRatio = nEmptyNarrow / nEmptyWide;
			
			double locEmptyRatioErr = sqrt(nEmptyNarrow)/nEmptyNarrow;
			
			locHistEmptyWide->Scale(locEmptyRatio);
			locHistEmptyWide_acc->Scale(locEmptyRatio);
			
			locHistEmpty->Add(locHistEmpty_acc, -0.1);
			locHistEmptyWide->Add(locHistEmptyWide_acc, -0.1);
			
			locFitter.SetEmpty(locHistEmpty, 1.0, 0.05);
			locFitter.SetEmptyWide(locHistEmptyWide);
		}
		
		//----------------------------------------------//
		// Do the fit:
		
		locFitter.FitData();
		locFitter.DumpFitParameters();
		
		//----------------------------------------------//
		// Save Results:
		
		// Integrated counts for full and empty target histograms:
		
		double locCounts, locCountsErr;
		IntegrateHistogram(locHistFull, locCounts, locCountsErr, 0.5, 0.6);
		
		double locCountsEmpty, locCountsEmptyErr;
		IntegrateHistogram(locHistEmpty, locCountsEmpty, locCountsEmptyErr, 0.5, 0.6);
		
		m_angularCounts[iThetaBin]      = {locCounts,      locCountsErr     };
		m_angularCountsEmpty[iThetaBin] = {locCountsEmpty, locCountsEmptyErr};
		
		// Yield of exclusive eta's estimated from histogram counts minus background fit functions:
		
		double locYield, locYieldErr;
		locFitter.GetYield(locYield, locYieldErr, 0, 1, 1);
		
		m_angularYield[iThetaBin] = {locYield, locYieldErr};
		
		// Yield of exclusive eta's estimated from integrating signal fit function
		
		double locYieldFit, locYieldFitErr;
		locFitter.GetYield(locYieldFit, locYieldFitErr, 1, 1);
		
		m_angularYieldFit[iThetaBin] = {locYieldFit, locYieldFitErr};
		
		// Yield of inclusive eta's estimated from integrating signal fit function
		
		double locYieldInc, locYieldIncErr;
		locFitter.GetYield(locYieldInc, locYieldIncErr, 1, 0);
		
		m_angularYieldInc[iThetaBin] = {locYieldInc, locYieldIncErr};
		
		// Empty target background integrated between mgg cut (from fit function):
		
		if(m_fitOption_empty) {
			double emptyYield, emptyYieldErr;
			locFitter.GetEmptyYield(emptyYield, emptyYieldErr);
			m_angularYieldEmpty[iThetaBin] = {emptyYield, emptyYieldErr};
		}
		
		// Only the peaking part of the empty target background integrated between mgg cut:
		
		if(m_fitOption_empty && (m_emptyFitOption_eta || m_emptyFitOption_fdc)) {
			double emptyYield, emptyYieldErr;
			locFitter.GetEmptyYield(emptyYield, emptyYieldErr, 1);
			m_angularYieldEmptyPeaking[iThetaBin] = {emptyYield, emptyYieldErr};
		}
		
		// Fraction of hadronic background integrated over all mgg:
		
		double hadronicBkgdFrac, hadronicBkgdFracErr;
		locFitter.GetHadronicBkgdFraction(hadronicBkgdFrac, hadronicBkgdFracErr);
		m_angularHadronicBkgdFraction[iThetaBin] = {hadronicBkgdFrac, hadronicBkgdFracErr};
		
		// eta+pion:
		double etaPionFrac, etaPionFracErr;
		locFitter.GetEtaPionBkgdFraction(etaPionFrac, etaPionFracErr);
		m_angularEtaPionBkgdFraction[iThetaBin] = {etaPionFrac, etaPionFracErr};
		
		// Yield of eta+pion background, integrated between mgg cut:
		
		double hadronicBkgdYield, hadronicBkgdYieldErr;
		locFitter.GetHadronicBkgdYield(hadronicBkgdYield, hadronicBkgdYieldErr);
		m_angularYieldHadronicBkgd[iThetaBin] = {hadronicBkgdYield, hadronicBkgdYieldErr};
		
		// eta+pion:
		double etaPionYield, etaPionYieldErr;
		locFitter.GetEtaPionYield(etaPionYield, etaPionYieldErr);
		m_angularYieldEtaPion[iThetaBin] = {etaPionYield, etaPionYieldErr};
		
		// Yield of omega background, integrated between mgg cut:
		
		double locYieldOmega, locYieldOmegaErr;
		locFitter.GetOmegaYield(locYieldOmega, locYieldOmegaErr);
		m_angularYieldOmega[iThetaBin] = {locYieldOmega, locYieldOmegaErr};
		
		// Get Omega fit parameters to monitor variation across all angular bins:
		
		double locOmegaMu,    locOmegaMuErr,    locOmegaSigma, locOmegaSigmaErr;
		double locOmegaAlpha, locOmegaAlphaErr, locOmegaN,     locOmegaNErr;
		locFitter.GetOmegaFitPars(locOmegaMu, locOmegaMuErr, locOmegaSigma, locOmegaSigmaErr, 
			locOmegaAlpha, locOmegaAlphaErr, locOmegaN, locOmegaNErr);
		m_omegaFitMu[iThetaBin]    = {locOmegaMu,    locOmegaMuErr};
		m_omegaFitSigma[iThetaBin] = {locOmegaSigma, locOmegaSigmaErr};
		m_omegaFitAlpha[iThetaBin] = {locOmegaAlpha, locOmegaAlphaErr};
		m_omegaFitN[iThetaBin]     = {locOmegaN,     locOmegaNErr};
		
		// Yield of remaining background (polynomial, exponential, chebyshev,...), integrated between mgg cut:
		
		double locYieldBkgd, locYieldBkgdErr;
		locFitter.GetBkgdYield(locYieldBkgd, locYieldBkgdErr);
		m_angularYieldBkgd[iThetaBin] = {locYieldBkgd, locYieldBkgdErr};
		
		// Shift of invariant mass lineshape fitted to the data:
		double locShift, locShiftErr;
		locFitter.GetLineshapeShift(locShift, locShiftErr);
		m_fitResult_shift[iThetaBin] = {locShift, locShiftErr};
		
		// Relative fraction of coherent and quasifree lineshapes from fit to data:
		double locZqf, locZqfErr;
		locFitter.GetQFFraction(locZqf, locZqfErr);
		m_fitResult_zqf[iThetaBin] = {locZqf, locZqfErr};
		
		//----------------------------------------------//
		
		//printf(" angle, yield1,  yield2 = %.2f,  %f,  %f\n", m_angularBin[iThetaBin].first, locYield, locYieldFit);
		printf(" angle: %f\n", m_angularBin[iThetaBin].first);
		printf("   integrated counts (full): %f\n", locCounts);
		printf("   integrated counts (empty): %f\n", locCountsEmpty);
		printf("   integrated counts (sub): %f\n", (locCounts-locCountsEmpty));
		printf("   yield from signal fit function integration: %f\n", locYieldFit);
		printf("   yield from background fit function integration: %f\n", locYield);
		printf("   hadronic background yield: %f\n", (m_angularYieldHadronicBkgd[iThetaBin].first+m_angularYieldEtaPion[iThetaBin].first));
		printf("   Uncertainty on eta yield: %.1f\%\n", 1.e2*locYieldFitErr/locYieldFit); 
		printf("\n");
		
		if(drawOption)
		{
			TH1F *locPull = (TH1F*)locHistFull->Clone("lochPull");
			styleMggHistogram(locPull);
			locFitter.FillPull(locPull);
			
			TF1 *locfFit    = (TF1*)locFitter.GetFitFunction()->Clone("locfFit");
			
			TF1 *locfSignal, *locfBkgd;
			TF1 *locfOmega = nullptr;
			TF1 *locfRho   = nullptr;
			TF1 *locfEmpty = nullptr;
			TF1 *locfHadronicBkgd = nullptr, *locfEtaPion = nullptr;
			
			locFitter.GetSignalFunction(&locfSignal, "locfSignal");
			locFitter.GetOmegaFunction(&locfOmega, "locfOmega");
			locFitter.GetBkgdFunction(&locfBkgd, "locfBkgd");
			
			if(m_fitOption_omega==4) locFitter.GetRhoFunction(&locfRho, "locfRho");
			
			locFitter.GetHadronicBkgdFunction(&locfHadronicBkgd, "locfHadronicBkgd");
			locFitter.GetEtaPionFunction(&locfEtaPion, "locfEtaPion");
			
			if(m_fitOption_empty==1) {
				// Change the binWidth parameter in the empty target fit function to match
				// the full target data we're plotting it over:
				locFitter.GetEmptyFitFunction(&locfEmpty, "locfEmpty", locHistFull->GetBinWidth(1));
			}
			DrawFitResult(locHistFull, locPull, locfFit, locfSignal, locfBkgd, locfOmega, locfRho, 
				locfHadronicBkgd, locfEtaPion, locfEmpty, locMinAngle, locMaxAngle);
			
			// Draw the empty target fit results separately:
			
			if(m_fitOption_empty==1) {
				
				locHistEmpty->GetXaxis()->SetRangeUser(m_minEmptyFitRange, m_maxEmptyFitRange);
				locHistEmptyWide->GetXaxis()->SetRangeUser(m_minEmptyFitRange, m_maxEmptyFitRange);
				
				if(cEmpty==NULL) InitializeEmptyCanvas();
				
				TF1 *locfEmptyDraw = nullptr;
				locFitter.GetEmptyFitFunction(&locfEmptyDraw, "locfEmptyDraw", locHistEmptyWide->GetBinWidth(1), 0);
				// the last parameter ('0') tells it to not exclude the accidental background from the returned function
				
				TF1 *locfEmptyBkgdDraw = nullptr;
				locFitter.GetEmptyBGFitFunction(&locfEmptyBkgdDraw, "locfEmptyBkgdDraw", locHistEmptyWide->GetBinWidth(1), 0);
				locfEmptyBkgdDraw->SetLineColor(kMagenta);
				
				styleMggHistogram(locHistEmpty, kBlue);
				styleMggHistogram(locHistEmptyWide, kCyan);
				
				cEmpty->cd();
				locHistEmptyWide->Draw("PE");
				locfEmptyDraw->Draw("same");
				locfEmptyBkgdDraw->Draw("same");
				//locHistEmpty->Draw("PE same");
				
				TLegend *locLeg = new TLegend(0.60, 0.60, 0.95, 0.89);
				locLeg->AddEntry(locHistEmpty,     Form("Empty Bkgd from %.2f#circ - %.2f#circ", locMinAngle,   locMaxAngle));
				locLeg->AddEntry(locHistEmptyWide, Form("Empty Bkgd from %.2f#circ - %.2f#circ", emptyAngleLow, emptyAngleHigh));
				//locLeg->Draw();
				
				TLatex latEmpty;
				latEmpty.SetTextColor(kBlack);
				latEmpty.SetTextFont(52);
				
				double emptyWindowSize = m_phase==1 ? 0.225 : 0.125;
				double minEmptyFitAngle = locAngle - emptyWindowSize;
				double maxEmptyFitAngle = locAngle + emptyWindowSize;
				if(minEmptyFitAngle < 0.0) {
					minEmptyFitAngle = 0.00;
					maxEmptyFitAngle = 2.0*emptyWindowSize;
				}
				latEmpty.DrawLatexNDC(0.635,0.85,Form("%.2f#circ < #theta_{#gamma#gamma}< %.2f#circ", minEmptyFitAngle, maxEmptyFitAngle));
				
				cEmpty->Update();
				cEmpty->Modified();
				//cEmpty->SaveAs("empty_fit.pdf");
				
				//cFit->SaveAs("full_target_fit.pdf");
				//cEmpty->SaveAs("empty_target_fit.pdf");
			}
			//getchar();
			//int character = getchar();
			//printf("The entered character is : %c", character);
			//int input;
			//cout << "Enter 'p' to print fit results to canvas... ";
			//cin  >> input;
			//cout << "Entered character is: " << input << endl;
		}
	}
	return;
}

void EtaAnalyzer::PlotAngularYield()
{
	h_Yield          = new TH1F("AngularYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_YieldFit       = new TH1F("AngularYieldFit", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_YieldInclusive = new TH1F("AngularYieldInclusive", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	h_Yield->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_Yield->GetXaxis()->SetTitleSize(0.05);
	h_Yield->GetXaxis()->SetTitleOffset(1.0);
	h_Yield->GetXaxis()->CenterTitle("");
	h_Yield->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_reconAngleBinSize));
	h_Yield->GetYaxis()->SetTitleSize(0.05);
	h_Yield->GetYaxis()->SetTitleOffset(1.0);
	h_Yield->GetYaxis()->CenterTitle("");
	h_Yield->SetTitle("");
	h_Yield->SetMarkerStyle(4);
	h_Yield->SetMarkerSize(0.7);
	h_Yield->SetMarkerColor(kBlue);
	h_Yield->SetLineColor(kBlue);
	h_Yield->SetLineWidth(2);
	
	h_YieldFit->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_YieldFit->GetXaxis()->SetTitleSize(0.05);
	h_YieldFit->GetXaxis()->SetTitleOffset(1.0);
	h_YieldFit->GetXaxis()->CenterTitle("");
	h_YieldFit->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_reconAngleBinSize));
	h_YieldFit->GetYaxis()->SetTitleSize(0.05);
	h_YieldFit->GetYaxis()->SetTitleOffset(1.0);
	h_YieldFit->GetYaxis()->CenterTitle("");
	h_YieldFit->SetTitle("");
	h_YieldFit->SetMarkerStyle(4);
	h_YieldFit->SetMarkerSize(0.7);
	h_YieldFit->SetMarkerColor(kRed);
	h_YieldFit->SetLineColor(kRed);
	h_YieldFit->SetLineWidth(2);
	
	h_YieldInclusive->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_YieldInclusive->GetXaxis()->SetTitleSize(0.05);
	h_YieldInclusive->GetXaxis()->SetTitleOffset(1.0);
	h_YieldInclusive->GetXaxis()->CenterTitle("");
	h_YieldInclusive->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_reconAngleBinSize));
	h_YieldInclusive->GetYaxis()->SetTitleSize(0.05);
	h_YieldInclusive->GetYaxis()->SetTitleOffset(1.0);
	h_YieldInclusive->GetYaxis()->CenterTitle("");
	h_YieldInclusive->SetTitle("");
	h_YieldInclusive->SetMarkerStyle(4);
	h_YieldInclusive->SetMarkerSize(0.7);
	h_YieldInclusive->SetMarkerColor(kBlack);
	h_YieldInclusive->SetLineColor(kBlack);
	h_YieldInclusive->SetLineWidth(2);
	
	double locMaxYield = 0.0;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		
		h_Yield->SetBinContent(ibin+1, m_angularYield[ibin].first);
		h_Yield->SetBinError(ibin+1, m_angularYield[ibin].second);
		
		h_YieldFit->SetBinContent(ibin+1, m_angularYieldFit[ibin].first);
		h_YieldFit->SetBinError(ibin+1, m_angularYieldFit[ibin].second);
		
		/*
		double locYieldInclusive = m_angularYieldFit[ibin].first 
			+ m_angularYieldEtaPion[ibin].first 
			+ m_angularYieldHadronicBkgd[ibin].first;
		double locYieldInclusiveErr = sqrt(pow(m_angularYieldFit[ibin].second,2.0)
			+ pow(m_angularYieldEtaPion[ibin].second,2.0)
			+ pow(m_angularYieldHadronicBkgd[ibin].second,2.0));
		
		h_YieldInclusive->SetBinContent(ibin+1, locYieldInclusive);
		h_YieldInclusive->SetBinError(ibin+1, locYieldInclusiveErr);
		*/
		h_YieldInclusive->SetBinContent(ibin+1, m_angularYieldInc[ibin].first);
		h_YieldInclusive->SetBinError(ibin+1, m_angularYieldInc[ibin].second);
		
		if(m_angularYieldInc[ibin].first>locMaxYield) {
			locMaxYield = m_angularYieldInc[ibin].first;
		}
	}
	
	h_Yield->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxYield);
	h_YieldFit->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxYield);
	h_YieldInclusive->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxYield);
	
	if(cYield==NULL) {
		cYield = new TCanvas("cYield", "Angular Yield", 950, 700);
		styleCanvas(cYield);
	}
	
	cYield->cd();
	h_YieldInclusive->Draw("PE1X0");
	h_Yield->Draw("PE1X0 same");
	h_YieldFit->Draw("PE1X0 same");
	cYield->Update();
	cYield->Modified();
	
	return;
}

void EtaAnalyzer::PlotCrossSection()
{
	h_CrossSection          = new TH1F("CrossSection", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_CrossSectionFit       = new TH1F("CrossSectionFit", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_CrossSectionInclusive = new TH1F("CrossSectionInclusive", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	h_CrossSection->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_CrossSection->GetXaxis()->SetTitleSize(0.05);
	h_CrossSection->GetXaxis()->SetTitleOffset(1.0);
	h_CrossSection->GetXaxis()->CenterTitle("");
	h_CrossSection->GetYaxis()->SetTitle("d#sigma/d#theta [#mub / rad.]");
	h_CrossSection->GetYaxis()->SetTitleSize(0.05);
	h_CrossSection->GetYaxis()->SetTitleOffset(1.0);
	h_CrossSection->GetYaxis()->CenterTitle("");
	h_CrossSection->SetTitle("");
	h_CrossSection->SetMarkerStyle(4);
	h_CrossSection->SetMarkerSize(0.7);
	h_CrossSection->SetMarkerColor(kBlue);
	h_CrossSection->SetLineColor(kBlue);
	h_CrossSection->SetLineWidth(2);
	
	h_CrossSectionFit->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_CrossSectionFit->GetXaxis()->SetTitleSize(0.05);
	h_CrossSectionFit->GetXaxis()->SetTitleOffset(1.0);
	h_CrossSectionFit->GetXaxis()->CenterTitle("");
	h_CrossSectionFit->GetYaxis()->SetTitle("d#sigma/d#theta [#mub / rad.]");
	h_CrossSectionFit->GetYaxis()->SetTitleSize(0.05);
	h_CrossSectionFit->GetYaxis()->SetTitleOffset(1.0);
	h_CrossSectionFit->GetYaxis()->CenterTitle("");
	h_CrossSectionFit->SetTitle("");
	h_CrossSectionFit->SetMarkerStyle(4);
	h_CrossSectionFit->SetMarkerSize(0.7);
	h_CrossSectionFit->SetMarkerColor(kRed);
	h_CrossSectionFit->SetLineColor(kRed);
	h_CrossSectionFit->SetLineWidth(2);
	
	h_CrossSectionInclusive->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_CrossSectionInclusive->GetXaxis()->SetTitleSize(0.05);
	h_CrossSectionInclusive->GetXaxis()->SetTitleOffset(1.0);
	h_CrossSectionInclusive->GetXaxis()->CenterTitle("");
	h_CrossSectionInclusive->GetYaxis()->SetTitle("d#sigma/d#theta [#mub / rad.]");
	h_CrossSectionInclusive->GetYaxis()->SetTitleSize(0.05);
	h_CrossSectionInclusive->GetYaxis()->SetTitleOffset(1.0);
	h_CrossSectionInclusive->GetYaxis()->CenterTitle("");
	h_CrossSectionInclusive->SetTitle("");
	h_CrossSectionInclusive->SetMarkerStyle(4);
	h_CrossSectionInclusive->SetMarkerSize(0.7);
	h_CrossSectionInclusive->SetMarkerColor(kBlack);
	h_CrossSectionInclusive->SetLineColor(kBlack);
	h_CrossSectionInclusive->SetLineWidth(2);
	
	double locMaxCS = 1.65;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		
		// Yield from integrated counts minus background:
		double locYield    = h_Yield->GetBinContent(ibin+1);
		double locYieldErr = h_Yield->GetBinError(ibin+1);
		
		// Yield from integrating signal component of fit function:
		double locYieldFit    = h_YieldFit->GetBinContent(ibin+1);
		double locYieldFitErr = h_YieldFit->GetBinError(ibin+1);
		
		// Yield from integration of signal+etapi+eta2pi components of fit function
		double locYieldInclusive    = h_YieldInclusive->GetBinContent(ibin+1);
		double locYieldInclusiveErr = h_YieldInclusive->GetBinError(ibin+1);
		
		double locAcc      = h_Acceptance->GetBinContent(ibin+1);
		double locAccErr   = h_Acceptance->GetBinError(ibin+1);
		
		if(locAcc <= 0.0) continue;
		double locBinSize = m_reconAngleBinSize * TMath::DegToRad(); // bin size in rad.
		
		//-------------------//
		
		double locCS    = locYield    / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		double locCSErr = locYieldErr / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		
		double locCS_upper = (locYield + locYieldErr) / ((locAcc-locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCS_lower = (locYield - locYieldErr) / ((locAcc+locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCS_error = (locCS_upper - locCS_lower) / 2.0;
		
		//-------------------//
		
		double locCSFit    = locYieldFit    / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFitErr = locYieldFitErr / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		
		double locCSFit_upper = (locYieldFit + locYieldFitErr) / ((locAcc-locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFit_lower = (locYieldFit - locYieldFitErr) / ((locAcc+locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFit_error = (locCSFit_upper - locCSFit_lower) / 2.0;
		
		//-------------------//
		
		double locCSInclusive    = locYieldInclusive    / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		double locCSInclusiveErr = locYieldInclusiveErr / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		
		double locCSInclusive_upper = (locYieldInclusive + locYieldInclusiveErr) 
			/ ((locAcc-locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSInclusive_lower = (locYieldInclusive - locYieldInclusiveErr) 
			/ ((locAcc+locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSInclusive_error = (locCSInclusive_upper - locCSInclusive_lower) / 2.0;
		
		//-------------------//
		
		if(0) {
			double locTheta = h_CrossSection->GetBinCenter(ibin+1);
			locCS          /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			locCS_error    /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			
			locCSFit       /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			locCSFit_error /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			
			locCSInclusive       /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			locCSInclusive_error /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
		}
		
		h_CrossSection->SetBinContent(ibin+1, locCS);
		h_CrossSection->SetBinError(ibin+1, locCS_error);
		
		h_CrossSectionFit->SetBinContent(ibin+1, locCSFit);
		h_CrossSectionFit->SetBinError(ibin+1, locCSFit_error);
		
		h_CrossSectionInclusive->SetBinContent(ibin+1, locCSInclusive);
		h_CrossSectionInclusive->SetBinError(ibin+1, locCSInclusive_error);
		
		if(locCSInclusive>locMaxCS) {
			locMaxCS = locCSInclusive;
		}
	}
	
	h_CrossSection->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxCS);
	h_CrossSectionFit->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxCS);
	h_CrossSectionInclusive->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxCS);
	
	if(cCrossSection==NULL) {
		cCrossSection = new TCanvas("cCrossSection", "Cross Section", 950, 700);
		styleCanvas(cCrossSection);
	}
	
	cCrossSection->cd();
	cCrossSection->SetGrid();
	h_CrossSectionInclusive->Draw("PE1X0");
	h_CrossSection->Draw("PE1X0 same");
	h_CrossSectionFit->Draw("PE1X0 same");
	cCrossSection->Update();
	cCrossSection->Modified();
	
	return;
}

void EtaAnalyzer::PlotEmptyEtaRatio()
{
	h_EmptyEtaRatio = new TH1F("EmptyEtaRatio", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	h_EmptyEtaRatio->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_EmptyEtaRatio->GetXaxis()->SetTitleSize(0.05);
	h_EmptyEtaRatio->GetXaxis()->SetTitleOffset(1.0);
	h_EmptyEtaRatio->GetXaxis()->CenterTitle("");
	h_EmptyEtaRatio->GetYaxis()->SetTitle("N_{#eta,gas} / N_{#eta}");
	h_EmptyEtaRatio->GetYaxis()->SetTitleSize(0.05);
	h_EmptyEtaRatio->GetYaxis()->SetTitleOffset(1.0);
	h_EmptyEtaRatio->GetYaxis()->CenterTitle("");
	h_EmptyEtaRatio->SetTitle("");
	h_EmptyEtaRatio->SetMarkerStyle(4);
	h_EmptyEtaRatio->SetMarkerSize(0.7);
	h_EmptyEtaRatio->SetMarkerColor(kBlue);
	h_EmptyEtaRatio->SetLineColor(kBlue);
	h_EmptyEtaRatio->SetLineWidth(2);
	
	double locMaxRatio = 0.0;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		
		//double num    = m_angularYieldEmpty[ibin].first;
		//double numErr = m_angularYieldEmpty[ibin].second;
		double num    = m_angularYieldEmptyPeaking[ibin].first;
		double numErr = m_angularYieldEmptyPeaking[ibin].second;
		
		double den    = (m_angularYieldFit[ibin].first + m_angularYieldHadronicBkgd[ibin].first + m_angularYieldEtaPion[ibin].first);
		double denErr = m_angularYield[ibin].second * (den/m_angularYield[ibin].first);
		
		double locRatio    = num/den;
		double locRatioErr = sqrt(pow(numErr/den,2.0) + pow(denErr*num/den/den,2.0));
		
		if(locRatio>locMaxRatio) locMaxRatio = locRatio;
		
		h_EmptyEtaRatio->SetBinContent(ibin+1, locRatio);
		h_EmptyEtaRatio->SetBinError(ibin+1, locRatioErr);
	}
	
	h_EmptyEtaRatio->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxRatio);
	
	if(cEmptyRatio==NULL) {
		cEmptyRatio = new TCanvas("cEmptyRatio", "Empty Ratio", 950, 700);
		styleCanvas(cEmptyRatio);
	}
	
	cEmptyRatio->cd();
	h_EmptyEtaRatio->Draw("PE1X0");
	h_EmptyEtaRatio->Draw("PE1X0 same");
	cEmptyRatio->Update();
	cEmptyRatio->Modified();
	
	return;
}

void EtaAnalyzer::PlotHadronicBkgdFraction()
{
	h_HadronicBkgdFraction = new TH1F("HadronicBkgdFraction", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	h_HadronicBkgdFraction->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_HadronicBkgdFraction->GetXaxis()->SetTitleSize(0.05);
	h_HadronicBkgdFraction->GetXaxis()->SetTitleOffset(1.0);
	h_HadronicBkgdFraction->GetXaxis()->CenterTitle("");
	h_HadronicBkgdFraction->GetYaxis()->SetTitle("N_{bkgd} / N_{#eta}");
	h_HadronicBkgdFraction->GetYaxis()->SetTitleSize(0.05);
	h_HadronicBkgdFraction->GetYaxis()->SetTitleOffset(1.0);
	h_HadronicBkgdFraction->GetYaxis()->CenterTitle("");
	h_HadronicBkgdFraction->SetTitle("");
	h_HadronicBkgdFraction->SetMarkerColor(kBlue);
	h_HadronicBkgdFraction->SetLineColor(kBlue);
	h_HadronicBkgdFraction->SetMarkerStyle(4);
	h_HadronicBkgdFraction->SetMarkerSize(0.7);
	h_HadronicBkgdFraction->SetLineWidth(2);
	
	h_HadronicBkgdFraction_bggen->SetLineColor(kMagenta);
	h_HadronicBkgdFraction_bggen->SetMarkerColor(kMagenta);
	h_HadronicBkgdFraction_bggen->SetMarkerStyle(4);
	h_HadronicBkgdFraction_bggen->SetMarkerSize(0.7);
	h_HadronicBkgdFraction_bggen->SetLineWidth(2);
	
	for(int ibin=1; ibin<=h_HadronicBkgdFraction_bggen->GetXaxis()->GetNbins(); ibin++) {
		h_HadronicBkgdFraction_bggen->SetBinContent(ibin, 1.0);
		h_HadronicBkgdFraction_bggen->SetBinError(ibin, 0.0);
	}
	
	double locMaxRatio = 0.0;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		
		double locRatio    = m_angularHadronicBkgdFraction[ibin].first;
		double locRatioErr = m_angularHadronicBkgdFraction[ibin].second;
		
		if(locRatio>locMaxRatio) locMaxRatio = locRatio;
		
		h_HadronicBkgdFraction->SetBinContent(ibin+1, locRatio);
		h_HadronicBkgdFraction->SetBinError(ibin+1, locRatioErr);
	}
	
	h_HadronicBkgdFraction->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxRatio);
	
	if(cHadronicBkgdFraction==NULL) {
		cHadronicBkgdFraction = new TCanvas("cHadronicBkgdFraction", "Hadronic Bkgd Fraction", 950, 700);
		styleCanvas(cHadronicBkgdFraction);
	}
	
	h_HadronicBkgdFraction->GetYaxis()->SetRangeUser(0.0,2.0);
	
	cHadronicBkgdFraction->cd();
	h_HadronicBkgdFraction->Draw("PE1X0");
	h_HadronicBkgdFraction_bggen->Draw("PE1X0 same");
	
	TLegend *locLeg = new TLegend(0.143, 0.774, 0.367, 0.916);
	locLeg->AddEntry(h_HadronicBkgdFraction_bggen, "BGGEN", "PE1X0");
	locLeg->AddEntry(h_HadronicBkgdFraction, "Fit Result", "PE1X0");
	locLeg->Draw();
	
	cHadronicBkgdFraction->Update();
	cHadronicBkgdFraction->Modified();
	return;
}

void EtaAnalyzer::PlotEtaPionFraction()
{
	h_EtaPionBkgdFraction = new TH1F("EtaPionBkgdFraction", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	h_EtaPionBkgdFraction->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_EtaPionBkgdFraction->GetXaxis()->SetTitleSize(0.05);
	h_EtaPionBkgdFraction->GetXaxis()->SetTitleOffset(1.0);
	h_EtaPionBkgdFraction->GetXaxis()->CenterTitle("");
	h_EtaPionBkgdFraction->GetYaxis()->SetTitle("N_{#eta+#pi} / N_{#eta}");
	h_EtaPionBkgdFraction->GetYaxis()->SetTitleSize(0.05);
	h_EtaPionBkgdFraction->GetYaxis()->SetTitleOffset(1.0);
	h_EtaPionBkgdFraction->GetYaxis()->CenterTitle("");
	h_EtaPionBkgdFraction->SetTitle("");
	h_EtaPionBkgdFraction->SetMarkerColor(kBlue);
	h_EtaPionBkgdFraction->SetLineColor(kBlue);
	h_EtaPionBkgdFraction->SetMarkerStyle(4);
	h_EtaPionBkgdFraction->SetMarkerSize(0.7);
	h_EtaPionBkgdFraction->SetLineWidth(2);
	
	h_EtaPionBkgdFraction_bggen->SetLineColor(kMagenta);
	h_EtaPionBkgdFraction_bggen->SetMarkerColor(kMagenta);
	h_EtaPionBkgdFraction_bggen->SetMarkerStyle(4);
	h_EtaPionBkgdFraction_bggen->SetMarkerSize(0.7);
	h_EtaPionBkgdFraction_bggen->SetLineWidth(2);
	
	double locMaxRatio = 0.0;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		
		double locRatio    = m_angularEtaPionBkgdFraction[ibin].first;
		double locRatioErr = m_angularEtaPionBkgdFraction[ibin].second;
		
		if(locRatio>locMaxRatio) locMaxRatio = locRatio;
		
		h_EtaPionBkgdFraction->SetBinContent(ibin+1, locRatio);
		h_EtaPionBkgdFraction->SetBinError(ibin+1, locRatioErr);
	}
	
	h_EtaPionBkgdFraction->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxRatio);
	
	if(cEtaPionFraction==NULL) {
		cEtaPionFraction = new TCanvas("cEtaPionFraction", "Eta+Pion Fraction", 950, 700);
		styleCanvas(cEtaPionFraction);
	}
	
	h_EtaPionBkgdFraction->GetYaxis()->SetRangeUser(0.0,2.0);
	
	cEtaPionFraction->cd();
	h_EtaPionBkgdFraction->Draw("PE1X0");
	h_EtaPionBkgdFraction_bggen->Draw("PE1X0 same");
	
	TLegend *locLeg = new TLegend(0.143, 0.774, 0.367, 0.916);
	locLeg->AddEntry(h_EtaPionBkgdFraction_bggen, "BGGEN", "PE1X0");
	locLeg->AddEntry(h_EtaPionBkgdFraction, "Fit Result", "PE1X0");
	locLeg->Draw();
	
	cEtaPionFraction->Update();
	cEtaPionFraction->Modified();
	
	return;
}

void EtaAnalyzer::StyleYieldHistogram(TH1F *h1, int markerStyle, int markerColor)
{
	h1->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitleOffset(1.0);
	h1->GetXaxis()->CenterTitle(true);
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	h1->SetTitle("");
	h1->SetMarkerStyle(markerStyle);
	h1->SetMarkerSize(0.7);
	h1->SetMarkerColor(markerColor);
	h1->SetLineColor(markerColor);
	h1->SetLineWidth(2);
}

void EtaAnalyzer::PlotOmegaFitPars()
{
	if(m_fitOption_omega!=1) return;
	
	h_omega_mu_fit    = new TH1F("omega_mu_fit", "",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_omega_sigma_fit = new TH1F("omega_sigma_fit", "",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_omega_alpha_fit = new TH1F("omega_alpha_fit", "",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_omega_n_fit     = new TH1F("omega_n_fit", "",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	StyleYieldHistogram(h_omega_mu_fit,    4);
	StyleYieldHistogram(h_omega_sigma_fit, 4);
	StyleYieldHistogram(h_omega_alpha_fit, 4);
	StyleYieldHistogram(h_omega_n_fit,     4);
	
	h_omega_mu_fit->GetYaxis()->SetTitle("#mu_{#omega} [GeV/c^{2}]");
	h_omega_sigma_fit->GetYaxis()->SetTitle("#sigma_{#omega} [GeV/c^{2}]");
	h_omega_alpha_fit->GetYaxis()->SetTitle("#alpha_{#omega}");
	h_omega_n_fit->GetYaxis()->SetTitle("n_{#omega}");
	
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) 
	{
		h_omega_mu_fit->SetBinContent(ibin+1, m_omegaFitMu[ibin].first);
		h_omega_mu_fit->SetBinError(ibin+1, m_omegaFitMu[ibin].second);
		
		h_omega_sigma_fit->SetBinContent(ibin+1, m_omegaFitSigma[ibin].first);
		h_omega_sigma_fit->SetBinError(ibin+1, m_omegaFitSigma[ibin].second);
		
		h_omega_alpha_fit->SetBinContent(ibin+1, m_omegaFitAlpha[ibin].first);
		h_omega_alpha_fit->SetBinError(ibin+1, m_omegaFitAlpha[ibin].second);
		
		h_omega_n_fit->SetBinContent(ibin+1, m_omegaFitN[ibin].first);
		h_omega_n_fit->SetBinError(ibin+1, m_omegaFitN[ibin].second);
	}
	
	h_omega_mu_fit->GetYaxis()->SetRangeUser(0.750, 0.790);
	h_omega_sigma_fit->GetYaxis()->SetRangeUser(0.024, 0.040);
	h_omega_alpha_fit->GetYaxis()->SetRangeUser(0.0, 2.0);
	h_omega_n_fit->GetYaxis()->SetRangeUser(0.0, 10.0);
	
	cOmegaFitPars = new TCanvas("cOmegaFitPars", "Omega Fit Pars", 1000, 750);
	cOmegaFitPars->Divide(2,2);
	
	cOmegaFitPars->cd(1);
	h_omega_mu_fit->Draw("PE1X0");
	
	cOmegaFitPars->cd(2);
	h_omega_sigma_fit->Draw("PE1X0");
	
	cOmegaFitPars->cd(3);
	h_omega_alpha_fit->Draw("PE1X0");
	
	cOmegaFitPars->cd(4);
	h_omega_n_fit->Draw("PE1X0");
	
	cOmegaFitPars->Update();
	cOmegaFitPars->Modified();
}

void EtaAnalyzer::PlotBackgrounds()
{
	h_Counts      = new TH1F("Counts", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_EmptyCounts = new TH1F("EmptyCounts", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	StyleYieldHistogram(h_Counts, 4, kBlack);
	StyleYieldHistogram(h_EmptyCounts, 4, kRed);
	
	h_Counts->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	h_EmptyCounts->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	
	//----------------------------------//
	
	h_EmptyYield = new TH1F("EmptyYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_OmegaYield = new TH1F("OmegaYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_BkgdYield = new TH1F("BkgdYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_HadronicBkgdYield = new TH1F("HadronicBkgdYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_EtaPionYield = new TH1F("EtaPionYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	
	StyleYieldHistogram(h_EmptyYield,        4, kRed);
	StyleYieldHistogram(h_OmegaYield,        4, kGreen);
	StyleYieldHistogram(h_BkgdYield,         4, kMagenta);
	StyleYieldHistogram(h_HadronicBkgdYield, 4, kTeal);
	StyleYieldHistogram(h_EtaPionYield,      4, kCyan);
	
	h_EmptyYield->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	h_OmegaYield->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	h_BkgdYield->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	h_HadronicBkgdYield->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	h_EtaPionYield->GetYaxis()->SetTitle(Form("counts / %.2f#circ", 2.0*m_angularBin[0].second));
	
	//----------------------------------//
	
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) 
	{
		h_Counts->SetBinContent(ibin+1, m_angularCounts[ibin].first);
		h_Counts->SetBinError(ibin+1, m_angularCounts[ibin].second);
		
		h_EmptyCounts->SetBinContent(ibin+1, m_angularCountsEmpty[ibin].first);
		h_EmptyCounts->SetBinError(ibin+1, m_angularCountsEmpty[ibin].second);
		
		h_EmptyYield->SetBinContent(ibin+1, m_angularYieldEmpty[ibin].first);
		h_EmptyYield->SetBinError(ibin+1, m_angularYieldEmpty[ibin].second);
		
		h_OmegaYield->SetBinContent(ibin+1, m_angularYieldOmega[ibin].first);
		h_OmegaYield->SetBinError(ibin+1, m_angularYieldOmega[ibin].second);
		
		h_BkgdYield->SetBinContent(ibin+1, m_angularYieldBkgd[ibin].first);
		h_BkgdYield->SetBinError(ibin+1, m_angularYieldBkgd[ibin].second);
		
		h_HadronicBkgdYield->SetBinContent(ibin+1, m_angularYieldHadronicBkgd[ibin].first);
		h_HadronicBkgdYield->SetBinError(ibin+1, m_angularYieldHadronicBkgd[ibin].second);
		
		h_EtaPionYield->SetBinContent(ibin+1, m_angularYieldEtaPion[ibin].first);
		h_EtaPionYield->SetBinError(ibin+1, m_angularYieldEtaPion[ibin].second);
	}
	
	if(cCounts==NULL) {
		cCounts = new TCanvas("cCounts", "Counts", 950, 700);
		styleCanvas(cCounts);
	}
	
	cCounts->cd();
	h_Counts->Draw("PE1X0");
	h_EmptyCounts->Draw("PE1X0 same");
	cCounts->Update();
	cCounts->Modified();
	
	if(cBackgrounds==NULL) {
		cBackgrounds = new TCanvas("cBackgrounds", "Backgrounds", 950, 700);
		styleCanvas(cBackgrounds);
	}
	
	TH1F *locCounts = (TH1F*)h_Counts->Clone("locCounts");
	locCounts->Add(h_EmptyCounts,-1.0);
	
	h_HadronicBkgdYield->SetLineColor(kRed);
	h_HadronicBkgdYield->SetMarkerColor(kRed);
	
	cBackgrounds->cd();
	locCounts->Draw("PE1X0");
	if(h_Yield) h_Yield->Draw("PE1X0 same");
	//h_EmptyYield->Draw("PE1X0 same");
	h_OmegaYield->Draw("PE1X0 same");
	h_BkgdYield->Draw("PE1X0 same");
	h_HadronicBkgdYield->Draw("PE1X0 same");
	h_EtaPionYield->Draw("PE1X0 same");
	cBackgrounds->Update();
	cBackgrounds->Modified();
	
	TLegend *legBg = new TLegend(0.6, 0.6, 0.8, 0.9);
	legBg->AddEntry(locCounts, "Integrated Counts (full-empty)", "PE1X0");
	legBg->AddEntry(h_Yield, "Extracted Yield", "PE1X0");
	legBg->AddEntry(h_OmegaYield, "Omega Bkgd", "PE1X0");
	legBg->AddEntry(h_BkgdYield, "Smooth Bkgd", "PE1X0");
	legBg->AddEntry(h_EtaPionYield, "Eta+Pion Bkgd", "PE1X0");
	legBg->AddEntry(h_HadronicBkgdYield, "Other Hadronic Bkgd", "PE1X0");
	legBg->Draw();
	
	return;
}

void EtaAnalyzer::PlotLineshapeShift()
{
	TH1F *h_LineshapeShift = new TH1F("LineshapeShift", "Shift of Lineshape Fitted to Data",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_LineshapeShift->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_LineshapeShift->GetXaxis()->SetTitleSize(0.05);
	h_LineshapeShift->GetXaxis()->SetTitleOffset(1.0);
	h_LineshapeShift->GetXaxis()->CenterTitle(true);
	h_LineshapeShift->GetYaxis()->SetTitle("#Delta#mu_{#eta} [GeV/c^{2}]");
	h_LineshapeShift->GetYaxis()->SetTitleSize(0.05);
	h_LineshapeShift->GetYaxis()->SetTitleOffset(1.0);
	h_LineshapeShift->GetYaxis()->CenterTitle(true);
	h_LineshapeShift->SetMarkerStyle(4);
	h_LineshapeShift->SetMarkerSize(0.8);
	h_LineshapeShift->SetMarkerColor(kBlue);
	h_LineshapeShift->SetLineColor(kBlue);
	
	h_LineshapeShift->GetYaxis()->SetRangeUser(-0.01,0.03);
	
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) 
	{
		h_LineshapeShift->SetBinContent(ibin+1, m_fitResult_shift[ibin].first);
		double locShiftErr = m_fitResult_shift[ibin].second;
		if(locShiftErr == 0.0) {
			locShiftErr = 0.002;
		}
		h_LineshapeShift->SetBinError(ibin+1, locShiftErr);
	}
	
	TCanvas *cShift = new TCanvas("cShift", "cShift", 950, 700);
	styleCanvas(cShift);
	h_LineshapeShift->Draw("PE1X0");
	cShift->Update();
	cShift->Modified();
	
	return;
}

void EtaAnalyzer::PlotQFFraction()
{
	TH1F *h_qfFraction = new TH1F("qfFraction", "Relative Fraction of QF Lineshape in Fit to Data",
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_qfFraction->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_qfFraction->GetXaxis()->SetTitleSize(0.05);
	h_qfFraction->GetXaxis()->SetTitleOffset(1.0);
	h_qfFraction->GetXaxis()->CenterTitle(true);
	h_qfFraction->GetYaxis()->SetTitle("z_{qf}");
	h_qfFraction->GetYaxis()->SetTitleSize(0.05);
	h_qfFraction->GetYaxis()->SetTitleOffset(1.0);
	h_qfFraction->GetYaxis()->CenterTitle(true);
	h_qfFraction->SetMarkerStyle(4);
	h_qfFraction->SetMarkerSize(0.8);
	h_qfFraction->SetMarkerColor(kBlue);
	h_qfFraction->SetLineColor(kBlue);
	
	h_qfFraction->GetYaxis()->SetRangeUser(0.0,1.0);
	
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) 
	{
		h_qfFraction->SetBinContent(ibin+1, m_fitResult_zqf[ibin].first);
		double locFracErr = m_fitResult_zqf[ibin].second;
		if(locFracErr == 0.0) {
			locFracErr = 0.01;
		}
		h_qfFraction->SetBinError(ibin+1, locFracErr);
	}
	
	TCanvas *cQFFrac = new TCanvas("cQFFrac", "cQFFrac", 950, 700);
	styleCanvas(cQFFrac);
	h_qfFraction->Draw("PE1X0");
	cQFFrac->Update();
	cQFFrac->Modified();
	
	return;
}

void EtaAnalyzer::IntegrateHistogram(TH1F *h1, double& counts, double& countsErr, double minCut, double maxCut)
{
	counts    = 0.0;
	countsErr = 0.0;
	
	// First check that the histogram binnning properly aligns with specified cut range, and provide warning if not:
	
	int locMinBin = h1->FindBin(minCut);
	int locMaxBin = h1->FindBin(maxCut)-1;
	
	double locBinSize = h1->GetXaxis()->GetBinWidth(1);
	double locMinCut  = h1->GetBinCenter(locMinBin) - 0.5*locBinSize;
	double locMaxCut  = h1->GetBinCenter(locMaxBin) + 0.5*locBinSize;
	
	if((fabs(locMinCut-minCut)>1.e-6) || (fabs(locMaxCut-maxCut)>1.e-6)) {
		printf("\n\nWarning inside 'IntegrateHistogram' function call: Cut ranges don't overlap with histogram bin edges.\n");
		printf("  Cut range: %f - %f\n", minCut, maxCut);
		printf("  Bin edges used in signal integration: %f - %f\n\n", locMinCut, locMaxCut);
	}
	
	for(int ibin=locMinBin; ibin<=locMaxBin; ibin++) {
		counts    += h1->GetBinContent(ibin);
		countsErr += pow(h1->GetBinError(ibin),2.0);
	}
	countsErr = sqrt(countsErr);
	
	return;
}

void EtaAnalyzer::InitializeBinning()
{
	double locAngle = m_minReconAngle + 0.5*m_reconAngleBinSize;
	
	while(locAngle < m_maxReconAngle) 
	{
		m_angularBin.push_back({locAngle, 0.5*m_reconAngleBinSize});
		
		m_angularCounts.push_back({0.0, 0.0});
		m_angularCountsEmpty.push_back({0.0, 0.0});
		
		m_angularYield.push_back({0.0, 0.0});
		m_angularYieldFit.push_back({0.0, 0.0});
		m_angularYieldInc.push_back({0.0, 0.0});
		
		m_angularYieldEmpty.push_back({0.0, 0.0});
		m_angularYieldEmptyPeaking.push_back({0.0, 0.0});
		
		m_angularHadronicBkgdFraction.push_back({0.0, 0.0});
		m_angularYieldHadronicBkgd.push_back({0.0, 0.0});
		
		m_angularEtaPionBkgdFraction.push_back({0.0, 0.0});
		m_angularYieldEtaPion.push_back({0.0, 0.0});
		
		m_angularYieldOmega.push_back({0.0, 0.0});
		m_angularYieldBkgd.push_back({0.0, 0.0});
		
		m_omegaFitMu.push_back({0.0, 0.0});
		m_omegaFitSigma.push_back({0.0, 0.0});
		m_omegaFitAlpha.push_back({0.0, 0.0});
		m_omegaFitN.push_back({0.0, 0.0});
		
		m_fitResult_shift.push_back({0.0, 0.0});
		m_fitResult_bkgdShift.push_back({0.0, 0.0});
		m_fitResult_zqf.push_back({0.0, 0.0});
		
		locAngle += m_reconAngleBinSize;
	}
	m_binningSet = true;
	return;
}

void InitializeMggFitter(MggFitter &fitter, EtaAnalyzer *anaObj, double angle, double angleWidth)
{
	fitter.angle            = angle;
	fitter.angleWidth       = angleWidth;
	fitter.empty_flux_ratio = anaObj->GetEmptyFluxRatio();
	fitter.binSize          = anaObj->GetMggBinSize();
	fitter.emptyBinSize     = anaObj->GetEmptyMggBinSize();
	fitter.fitOption_signal = anaObj->GetFitOption(1);
	fitter.fitOption_bkgd   = anaObj->GetFitOption(2);
	fitter.fitOption_poly   = anaObj->GetFitOption(3);
	fitter.fitOption_omega  = anaObj->GetFitOption(4);
	fitter.fitOption_etap   = anaObj->GetFitOption(5);
	anaObj->GetFitRange(fitter.minFitRange, fitter.maxFitRange);
	
	fitter.useRawMass = anaObj->GetRawMassOption();
	fitter.vetoOption = anaObj->GetVetoOption();
	
	// only fit the empty target if we don't subtract it first:
	if(anaObj->GetEmptySubtractOption()==0) {
		fitter.fitOption_empty  = anaObj->GetEmptyFitOption(0);
	}
	
	fitter.emptyFitOption_eta   = anaObj->GetEmptyFitOption(1);
	fitter.emptyFitOption_omega = anaObj->GetEmptyFitOption(2);
	fitter.emptyFitOption_fdc   = anaObj->GetEmptyFitOption(3);
	fitter.emptyFitOption_bkgd  = anaObj->GetEmptyFitOption(4);
	fitter.emptyFitOption_poly  = anaObj->GetEmptyFitOption(5);
	anaObj->GetEmptyFitRange(fitter.minEmptyFitRange, fitter.maxEmptyFitRange);
	return;
}

void EtaAnalyzer::WriteROOTFile(TString fileName)
{
	TFile *fOut = new TFile(fileName.Data(), "RECREATE");
	fOut->cd();
	
	// Angular Yield:
	if(h_Yield) h_Yield->Write();
	if(h_YieldFit) h_YieldFit->Write();
	if(h_YieldInclusive) h_YieldInclusive->Write();
	
	// Cross Section:
	if(h_CrossSection) h_CrossSection->Write();
	if(h_CrossSectionFit) h_CrossSectionFit->Write();
	if(h_CrossSectionInclusive) h_CrossSectionInclusive->Write();
	
	// Counts:
	if(h_Counts) h_Counts->Write();
	if(h_EmptyCounts) h_EmptyCounts->Write();
	
	// Fitted empty yield:
	if(h_EmptyYield) h_EmptyYield->Write();
	
	// Fitted background Yields:
	if(h_OmegaYield) h_OmegaYield->Write(); // omega peaking background
	if(h_BkgdYield) h_BkgdYield->Write(); // smooth (e.g. polynomial, expo, cheb)
	if(h_EtaPionYield) h_EtaPionYield->Write();
	if(h_HadronicBkgdYield) h_HadronicBkgdYield->Write();
	
	if(h_Acceptance) h_Acceptance->Write();
	
	fOut->Write();
	fOut->Close();
	
	return;
}

void EtaAnalyzer::DumpSettings() 
{
	printf("\n\n");
	printf("=================================================\n");
	printf("Extracting eta->gg angular yield from phase %d.\n\n", m_phase);
	printf("  Beam Energy Range: %.2f GeV - %.2f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
	printf("  Angular binning: %.2f degrees\n", m_reconAngleBinSize);
	printf("  Fit functions:\n");
	printf("    Signal: %s\n", GetFitOptionStr(0).Data());
	printf("    Empty: %s\n", GetEmptyFitOptionStr().Data());
	printf("    Background: %s\n", GetFitOptionStr(1).Data());
	printf("    Omega: %s\n", GetFitOptionStr(2).Data());
	printf("    Eta prime: %s\n", GetFitOptionStr(3).Data());
	printf("    Fitting range: %.3f GeV - %.3f GeV\n", m_minFitRange, m_maxFitRange);
	printf("\n=================================================\n");
	
	return;
}
