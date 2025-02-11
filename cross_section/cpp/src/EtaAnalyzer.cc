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
	if((option<0) || (option>5)) {
		cout << "\nUnsupported analysis option provided.\n" << endl;
		exit(1);
	}
	m_analysisOption = option;
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

void EtaAnalyzer::SetBeamEnergy(double min, double max)
{
	m_minBeamEnergy = min;
	m_maxBeamEnergy = max;
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
	if((option<1) || (option>7)) {
		printf("\nUnsupported signal fit option provided (should be 1-6). Using default option: %d\n", m_fitOption_signal);
	}
	else {
		m_fitOption_signal = option;
	}
	return;
}
void EtaAnalyzer::SetFitOption_bkgd(int option) 
{
	if((option<1) || (option>5)) {
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
	if((option<1) || (option>3)) {
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
	if((option<0) || (option>3)) {
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
	double locMin = 0.40;
	if(m_phase==1) locMin = 0.30;
	else           locMin = 0.40;
	
	m_minFitRange = min < locMin ? locMin : min;
	m_maxFitRange = max;
	return;
}

void EtaAnalyzer::SetEmptyFitRange(double min, double max)
{
	double locMin = 0.40;
	if(m_phase==1) locMin = 0.30;
	else           locMin = 0.40;
	
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
					optString = "single Gaussian";
					break;
				case 2:
					optString = "double Gaussian";
					break;
				case 3:
					optString = "Crystal Ball";
					break;
				case 4:
					optString = "Crystal Ball + Gaussian";
					break;
				case 5:
					optString = "Simulated Lineshape";
					break;
				case 6:
					optString = "Simulated Lineshape + Gaussian";
					break;
			}
			break;
		case 1:
			switch(m_fitOption_bkgd) {
				case 1:
					optString = Form("polynomial (order %d)", m_fitOption_poly);
					break;
				case 2:
					optString = "exponential";
					break;
				case 3:
					optString = Form("Chebyshev polynomial (order %d)", m_fitOption_poly);
					break;
				case 4:
					optString = "no background";
					break;
				case 5:
					optString = "empty target lineshape";
					break;
			}
			break;
		case 2:
			switch(m_fitOption_omega) {
				case 1:
					optString = "Crystal Ball (floating parameters)";
					break;
				case 2:
					optString = "Simulated Lineshape (with double Crystal Ball parameterization)";
					break;
				case 3:
					optString = "Simulated Lineshape";
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

//---------------------------------------------------------//

void EtaAnalyzer::InitializeFitCanvas()
{
	cFit = new TCanvas("cFit", "Mgg Fit", 1200, 900);
	
	pFit = new TPad("padFit", "Mgg Fit", 0.005, 0.3025, 0.995, 0.995);
	pRes = new TPad("padRes", "Mgg Res", 0.005, 0.005,  0.995, 0.2975);
	
	pFit->SetLeftMargin(0.10);
	pFit->SetRightMargin(0.02);
	pFit->SetTopMargin(0.075);
	pFit->SetBottomMargin(0.015);
	pFit->SetTickx(); pFit->SetTicky();
	pFit->SetFrameLineWidth(2);
	
	pRes->SetLeftMargin(0.10);
	pRes->SetRightMargin(0.02);
	pRes->SetTopMargin(0.005);
	pRes->SetBottomMargin(0.325);
	pRes->SetTickx(); pRes->SetTicky();
	pRes->SetFrameLineWidth(2);
	
	cFit->cd();
	pFit->Draw();
	pRes->Draw();
	
	return;
}

void EtaAnalyzer::InitializeEmptyCanvas()
{
	cEmpty = new TCanvas("cEmpty", "Empty Mgg Fit", 950, 700);
	styleCanvas(cEmpty);
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

void EtaAnalyzer::FitInvariantMass(double minAngle, double maxAngle, int drawFitResult, int drawOmegaFit)
{
	int minAngleBin = h_mggVsThetaFull->GetXaxis()->FindBin(minAngle);
	int maxAngleBin = h_mggVsThetaEmpty->GetXaxis()->FindBin(maxAngle)-1;
	
	TH1F *locHistFull  = (TH1F*)h_mggVsThetaFull->ProjectionY("locHistFull", minAngleBin, maxAngleBin);
	TH1F *locHistEmpty = (TH1F*)h_mggVsThetaEmpty->ProjectionY("locHistEmpty", minAngleBin, maxAngleBin);
	
	locHistFull->Rebin(m_rebinsMgg);
	locHistEmpty->Rebin(m_rebinsMgg);
	
	styleMggHistogram(locHistFull);
	styleMggHistogram(locHistEmpty, kBlue);
	
	
	MggFitter locFitter;
	InitializeMggFitter(locFitter, this, 0.5*(minAngle+maxAngle));
	
	// Load 1-D lineshape projections:
	
	if(h_etaLineshape) {
		minAngleBin = h_etaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_etaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hEta;
		if(m_fitOption_signal==7) {
			hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta");
			hEta->Rebin(m_rebinsMgg);
		}
		else {
			hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
		}
		locFitter.SetEtaLineshape(hEta);
	}
	
	if(h_omegaLineshape) {
		minAngleBin = h_omegaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_omegaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega");//, minAngleBin, maxAngleBin);
		//hOmega->Rebin(m_rebinsMgg);
		locFitter.SetOmegaLineshape(hOmega, drawOmegaFit);
	}
	
	if(h_fdcOmegaLineshape) {
		minAngleBin = h_fdcOmegaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_fdcOmegaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hFDCOmega = (TH1F*)h_fdcOmegaLineshape->ProjectionY("hFDCOmega");//, minAngleBin, maxAngleBin);
		//hFDCOmega->Rebin(m_rebinsMgg);
		locFitter.SetFDCOmegaLineshape(hFDCOmega, drawOmegaFit);
	}
	
	if(h_etaPionLineshape && (m_fitOption_signal==7)) {
		TH1F *hEtaPion = (TH1F*)h_etaPionLineshape->ProjectionY("hEtaPion");
		hEtaPion->Rebin(m_rebinsMgg);
		locFitter.SetEtaPionLineshape(hEtaPion);
	}
	
	locFitter.SetData(locHistFull);
	locFitter.SetEmpty(locHistEmpty);
	
	locFitter.FitEmpty();
	
	locHistFull->GetXaxis()->SetRangeUser(m_minFitRange, m_maxFitRange);
	
	if(cFit==NULL) InitializeFitCanvas();
	
	TH1F *locEmptyPull = (TH1F*)locHistEmpty->Clone("lochEmptyPull");
	locFitter.FillEmptyPull(locEmptyPull);
	TF1 *locfEmpty = (TF1*)locFitter.GetEmptyFitFunction()->Clone("locfEmpty");
	
	pFit->cd();
	locHistFull->Draw("PE1");
	//locHistEmpty->Draw("PE1");
	locfEmpty->Draw("same");
	pRes->cd();
	locEmptyPull->Draw("PE1");
	cFit->Update();
	cFit->Modified();
	locFitter.DumpEmptyFitParameters();
	getchar();
	/*
	locHistFull->Add(locHistEmpty, -1.0);
	
	MggFitter locFitter;
	InitializeMggFitter(locFitter, this, 0.5*(minAngle+maxAngle));
	
	// Load 1-D lineshape projections:
	
	if(h_etaLineshape && (m_fitOption_signal>=5)) {
		minAngleBin = h_etaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_etaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hEta;
		if(m_fitOption_signal==7) {
			hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta");
			hEta->Rebin(m_rebinsMgg);
		}
		else {
			hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
		}
		locFitter.SetEtaLineshape(hEta);
	}
	
	if(h_omegaLineshape) {
		minAngleBin = h_omegaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_omegaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega");//, minAngleBin, maxAngleBin);
		//hOmega->Rebin(m_rebinsMgg);
		locFitter.SetOmegaLineshape(hOmega, drawOmegaFit);
	}
	
	if(h_etaPionLineshape && (m_fitOption_signal==7)) {
		TH1F *hEtaPion = (TH1F*)h_etaPionLineshape->ProjectionY("hEtaPion");
		hEtaPion->Rebin(m_rebinsMgg);
		locFitter.SetEtaPionLineshape(hEtaPion);
	}
	
	locFitter.SetData(locHistFull);
	locFitter.FitData();
	
	locFitter.DumpFitParameters();
	
	double yield, yieldErr;
	locFitter.GetYield(yield, yieldErr);
	
	double yieldFit, yieldFitErr;
	locFitter.GetYield(yieldFit, yieldFitErr, 1);
	
	printf("\n");
	printf("Yield (method 1) = %f +/- %f\n", yield,    yieldErr);
	printf("Yield (method 2) = %f +/- %f\n", yieldFit, yieldFitErr);
	
	if(drawFitResult)
	{
		TH1F *locPull = (TH1F*)locHistFull->Clone("lochPull");
		locFitter.FillPull(locPull);
		TF1 *locfFit     = (TF1*)locFitter.GetFitFunction()->Clone("locfFit");
		TF1 *locfSignal  = (TF1*)locFitter.GetSignalFunction()->Clone("locfSignal");
		TF1 *locfBkgd    = (TF1*)locFitter.GetBkgdFunction()->Clone("locfBkgd");
		TF1 *locfEtaPion = (TF1*)locFitter.GetEtaPionFunction()->Clone("locfEtaPion");
		
		double locDiff = 0.0;
		for(int ibin=1; ibin<=locHistFull->GetXaxis()->GetNbins(); ibin++) {
			double locCounts = locHistFull->GetBinContent(ibin);
			double locBkgd   = locfFit->Eval(locHistFull->GetXaxis()->GetBinCenter(ibin));
			double locX = locHistFull->GetXaxis()->GetBinCenter(ibin);
			if(locX>0.5 && locX<0.6) {
				locDiff += (locCounts - locBkgd);
			}
			//locHistFull->SetBinContent(ibin, locCounts-locBkgd);
		}
		cout << "Differnence between data and fit: " << locDiff << endl;
		
		DrawFitResult(locHistFull, locPull, locfFit, locfSignal, locfBkgd, locfEtaPion);
	}
	*/
	return;
}

void EtaAnalyzer::DrawFitResult(TH1F *h1, TH1F *hPull, TF1 *fFit, TF1 *fSignal, TF1 *fBkgd, TF1 *fEtaPi, TF1 *fEmpty, 
	double minAngle, double maxAngle)
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
	
	pFit->cd();
	h1->Draw("PE1");
	fFit->Draw("same");
	fSignal->Draw("same");
	fBkgd->Draw("same");
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
	
	TLegend *locLeg = new TLegend(0.135, 0.550, 0.350, 0.800);
	locLeg->AddEntry(fFit,    "Full Fit",         "l");
	locLeg->AddEntry(fSignal, "Signal Lineshape", "l");
	if(fEmpty) locLeg->AddEntry(fEmpty, "Empty-Target Bkgd", "l");
	if(fEtaPi) locLeg->AddEntry(fEtaPi, "Eta+Pion Lineshape", "l");
	locLeg->AddEntry(fBkgd,   Form("Additional Bkgd (%s)", GetBkgdFitName().Data()), "l");
	locLeg->Draw();
	
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
	return;
}

void EtaAnalyzer::ExtractAngularYield(int drawOption)
{
	if(m_binningSet==false) InitializeBinning();
	
	for(int iThetaBin=0; iThetaBin<m_angularBin.size(); iThetaBin++) {
		
		double locAngle    = m_angularBin[iThetaBin].first;
		double locMinAngle = locAngle - m_angularBin[iThetaBin].second;
		double locMaxAngle = locAngle + m_angularBin[iThetaBin].second;
		
		//if(locAngle<0.26 || locAngle>0.28) continue;
		
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
		
		TH1F *locHistFull  = (TH1F*)h_mggVsThetaFull->ProjectionY("locHistFull", minAngleBin, maxAngleBin);
		locHistFull->Rebin(m_rebinsMgg);
		styleMggHistogram(locHistFull);
		
		TH1F *locHistEmpty = (TH1F*)h_mggVsThetaEmpty->ProjectionY("locHistEmpty", minAngleBin, maxAngleBin);
		locHistEmpty->Rebin(m_rebinsMgg);
		styleMggHistogram(locHistEmpty, kBlue);
		
		if(m_subtractEmpty) {
			locHistFull->Add(locHistEmpty,-1.0);
		}
		
		// set bin errors for empty bins to 1, otherwise they are ignored in a chi-squared minimization:
		for(int ibin=1; ibin<=locHistFull->GetXaxis()->GetNbins(); ibin++) {
			if(locHistFull->GetBinContent(ibin)==0.0) {
				if(locHistFull->GetBinError(ibin)==0.0) {
					locHistFull->SetBinError(ibin, 1.0);
				}
			}
		}
		
		//----------------------------------------------//
		// Set up fitter object:
		
		MggFitter locFitter;
		InitializeMggFitter(locFitter, this, m_angularBin[iThetaBin].first);
		locFitter.SetData(locHistFull);
		
		double locBinSize = 0.0;
		if((h_etaLineshape) && (m_fitOption_signal>=5)) {
			
			locBinSize  = h_etaLineshape->GetXaxis()->GetBinWidth(1);
			minAngleBin = h_etaLineshape->GetXaxis()->FindBin(locMinAngle + 0.5*locBinSize);
			maxAngleBin = h_etaLineshape->GetXaxis()->FindBin(locMaxAngle - 0.5*locBinSize)-1;
			
			locMinAngleHist = h_etaLineshape->GetXaxis()->GetBinCenter(minAngleBin) 
				- 0.5*h_etaLineshape->GetXaxis()->GetBinWidth(1);
			locMaxAngleHist = h_etaLineshape->GetXaxis()->GetBinCenter(maxAngleBin) 
				+ 0.5*h_etaLineshape->GetXaxis()->GetBinWidth(1);
			if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
				printf("\nWarning: Eta lineshape bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
				printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
				printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
			}
			
			/*
			TH1F *hEta;
			if(m_fitOption_signal==7) {
				hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta");
				hEta->Rebin(m_rebinsMgg);
			}
			else {
				hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
			}
			*/
			TH1F *hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
			locFitter.SetEtaLineshape(hEta);
		}
		if(h_omegaLineshape) {
			
			locBinSize  = h_omegaLineshape->GetXaxis()->GetBinWidth(1);
			minAngleBin = h_omegaLineshape->GetXaxis()->FindBin(locMinAngle + 0.5*locBinSize);
			maxAngleBin = h_omegaLineshape->GetXaxis()->FindBin(locMaxAngle - 0.5*locBinSize)-1;
			
			locMinAngleHist = h_omegaLineshape->GetXaxis()->GetBinCenter(minAngleBin) 
				- 0.5*h_omegaLineshape->GetXaxis()->GetBinWidth(1);
			locMaxAngleHist = h_omegaLineshape->GetXaxis()->GetBinCenter(maxAngleBin) 
				+ 0.5*h_omegaLineshape->GetXaxis()->GetBinWidth(1);
			if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
				printf("\nWarning: Omega lineshape bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
				printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
				printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
			}
			
			TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega");//, minAngleBin, maxAngleBin);
			locFitter.SetOmegaLineshape(hOmega);
		}
		if((h_etaPionLineshape) && (m_fitOption_signal==7)) {
			
			locBinSize  = h_etaPionLineshape->GetXaxis()->GetBinWidth(1);
			minAngleBin = h_etaPionLineshape->GetXaxis()->FindBin(locMinAngle + 0.5*locBinSize);
			maxAngleBin = h_etaPionLineshape->GetXaxis()->FindBin(locMaxAngle - 0.5*locBinSize)-1;
			
			locMinAngleHist = h_etaPionLineshape->GetXaxis()->GetBinCenter(minAngleBin) 
				- 0.5*h_etaPionLineshape->GetXaxis()->GetBinWidth(1);
			locMaxAngleHist = h_etaPionLineshape->GetXaxis()->GetBinCenter(maxAngleBin) 
				+ 0.5*h_etaPionLineshape->GetXaxis()->GetBinWidth(1);
			if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
				printf("\nWarning: Eta+Pion lineshape bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
				printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
				printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
			}
			
			TH1F *hEtaPion = (TH1F*)h_etaPionLineshape->ProjectionY("hEtaPion");
			hEtaPion->Rebin(m_rebinsMgg);
			locFitter.SetEtaPionLineshape(hEtaPion);
		}
		
		if((h_fdcOmegaLineshape) && (m_subtractEmpty==0) && (m_fitOption_empty==1) && (m_emptyFitOption_fdc>=2)) {
			
			locBinSize  = h_fdcOmegaLineshape->GetXaxis()->GetBinWidth(1);
			minAngleBin = h_fdcOmegaLineshape->GetXaxis()->FindBin(locMinAngle + 0.5*locBinSize);
			maxAngleBin = h_fdcOmegaLineshape->GetXaxis()->FindBin(locMaxAngle - 0.5*locBinSize)-1;
			
			locMinAngleHist = h_fdcOmegaLineshape->GetXaxis()->GetBinCenter(minAngleBin) 
				- 0.5*h_fdcOmegaLineshape->GetXaxis()->GetBinWidth(1);
			locMaxAngleHist = h_fdcOmegaLineshape->GetXaxis()->GetBinCenter(maxAngleBin) 
				+ 0.5*h_fdcOmegaLineshape->GetXaxis()->GetBinWidth(1);
			if((fabs(locMinAngle-locMinAngleHist)>1.e-6) || (fabs(locMinAngle-locMinAngleHist)>1.e-6)) {
				printf("\nWarning: FDC lineshape bin edges do not overlap with minimum and/or maximum angular bin ranges.\n");
				printf("  Desired bin range: %f-%f\n", locMinAngle, locMaxAngle);
				printf("  Histogram bin range: %f-%f\n", locMinAngleHist, locMaxAngleHist);
			}
			
			TH1F *hFDCOmega = (TH1F*)h_fdcOmegaLineshape->ProjectionY("hFDCOmega");//, minAngleBin, maxAngleBin);
			//hFDCOmega->Rebin(m_rebinsMgg);
			locFitter.SetFDCOmegaLineshape(hFDCOmega);
		}
		
		//----------------------------------------------//
		// Do fit and extract yield:
		
		TH1F *locHistEmptyWide;
		double emptyAngleLow  = locMinAngle;
		double emptyAngleHigh = locMaxAngle;
		
		if((m_subtractEmpty==0) && (m_fitOption_empty==1)) {
			
			// To get the pdf of the empty target background, we need to combine a wider angular range.
			// For now, we'll use +/-0.25 degrees from the central bin:
			
			double emptyWindowSize = m_phase==1 ? 0.25 : 0.15;
			
			emptyAngleLow  = locAngle - emptyWindowSize;
			emptyAngleHigh = locAngle + emptyWindowSize;
			if(emptyAngleLow < 0.0) {
				emptyAngleLow  = 0.00;
				emptyAngleHigh = 2.0*emptyWindowSize;
			}
			
			locHistEmptyWide = (TH1F*)h_mggVsThetaEmpty->ProjectionY("EmptyHistWide",
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleLow),
				h_mggVsThetaEmpty->GetXaxis()->FindBin(emptyAngleHigh)-1, "e");
			locHistEmptyWide->Rebin(m_rebinsMgg);
			styleMggHistogram(locHistEmptyWide, kCyan);
			
			double nEmptyNarrow  = locHistEmpty->Integral(
				locHistEmpty->GetXaxis()->FindBin(m_minEmptyFitRange), locHistEmpty->GetXaxis()->FindBin(m_maxEmptyFitRange));
			double nEmptyWide    = locHistEmptyWide->Integral(
				locHistEmptyWide->GetXaxis()->FindBin(m_minEmptyFitRange), locHistEmptyWide->GetXaxis()->FindBin(m_maxEmptyFitRange));
			double locEmptyRatio = nEmptyNarrow / nEmptyWide;
			
			double locEmptyRatioErr = sqrt(nEmptyNarrow)/nEmptyNarrow;
			
			locHistEmptyWide->Scale(locEmptyRatio);
			
			locFitter.SetEmpty(locHistEmptyWide, 1.0, locEmptyRatioErr);
			locFitter.FitEmpty();
		}
		
		locFitter.FitData();
		
		double locYield, locYieldErr;
		locFitter.GetYield(locYield, locYieldErr);
		
		double locYieldFit, locYieldFitErr;
		locFitter.GetYield(locYieldFit, locYieldFitErr, 1);
		
		// Inflate error bars from empty target background:
		double nEmpty = locHistEmpty->Integral(locHistEmpty->FindBin(0.5), locHistEmpty->FindBin(0.6)) / m_emptyTargetFluxRatio;
		locYieldErr = sqrt(locYield + pow(m_emptyTargetFluxRatio,2.0)*nEmpty);
		
		m_angularYield[iThetaBin]    = {locYield, locYieldErr};
		m_angularYieldFit[iThetaBin] = {locYieldFit, locYieldErr};
		
		if(drawOption)
		{
			TH1F *locPull = (TH1F*)locHistFull->Clone("lochPull");
			locFitter.FillPull(locPull);
			TF1 *locfFit    = (TF1*)locFitter.GetFitFunction()->Clone("locfFit");
			TF1 *locfSignal = (TF1*)locFitter.GetSignalFunction()->Clone("locfSignal");
			TF1 *locfBkgd   = (TF1*)locFitter.GetBkgdFunction()->Clone("locfBkgd");
			TF1 *locfEtaPi  = NULL;
			TF1 *locfEmpty  = NULL;
			
			if(m_fitOption_signal==7) {
				locfEtaPi = (TF1*)locFitter.GetEtaPionFunction()->Clone("locfEtaPi");
			}
			if((m_fitOption_empty==1) && (m_subtractEmpty==0)) {
				locfEmpty = (TF1*)locFitter.GetEmptyFitFunction()->Clone("locfEmpty");
			}
			DrawFitResult(locHistFull, locPull, locfFit, locfSignal, locfBkgd, locfEtaPi, locfEmpty, locMinAngle, locMaxAngle);
			
			// Draw the empty target fit results separately:
			
			if((m_fitOption_empty==1) && (m_subtractEmpty==0)) {
				
				locHistEmpty->GetXaxis()->SetRangeUser(m_minEmptyFitRange, m_maxEmptyFitRange);
				
				if(cEmpty==NULL) InitializeEmptyCanvas();
				
				cEmpty->cd();
				locHistEmpty->Draw("PE");
				locHistEmptyWide->Draw("PE same");
				locfEmpty->Draw("same");
				
				TLegend *locLeg = new TLegend(0.60, 0.60, 0.95, 0.89);
				locLeg->AddEntry(locHistEmpty,     Form("Empty Bkgd from %.2f#circ - %.2f#circ", locMinAngle,   locMaxAngle));
				locLeg->AddEntry(locHistEmptyWide, Form("Empty Bkgd from %.2f#circ - %.2f#circ", emptyAngleLow, emptyAngleHigh));
				locLeg->Draw();
				
				cEmpty->Update();
				cEmpty->Modified();
			}
		}
		printf(" angle, yield1,  yield2 = %f,  %f,  %f\n", m_angularBin[iThetaBin].first, locYield, locYieldFit);
	}
	
	return;
}

void EtaAnalyzer::PlotAngularYield()
{
	h_Yield    = new TH1F("AngularYield", "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_YieldFit = new TH1F("AngularYieldFit", "", 
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
	
	double locMaxYield = 0.0;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		h_Yield->SetBinContent(ibin+1, m_angularYield[ibin].first);
		h_Yield->SetBinError(ibin+1, m_angularYield[ibin].second);
		
		h_YieldFit->SetBinContent(ibin+1, m_angularYieldFit[ibin].first);
		h_YieldFit->SetBinError(ibin+1, m_angularYieldFit[ibin].second);
		
		if(m_angularYield[ibin].first>locMaxYield) {
			locMaxYield = m_angularYield[ibin].first;
		}
		if(m_angularYieldFit[ibin].first>locMaxYield) {
			locMaxYield = m_angularYieldFit[ibin].first;
		}
	}
	
	h_Yield->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxYield);
	h_YieldFit->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxYield);
	
	if(cYield==NULL) {
		cYield = new TCanvas("cYield", "Angular Yield", 950, 700);
		styleCanvas(cYield);
	}
	
	cYield->cd();
	h_Yield->Draw("PE1X0");
	h_YieldFit->Draw("PE1X0 same");
	cYield->Update();
	cYield->Modified();
	
	return;
}

void EtaAnalyzer::PlotCrossSection()
{
	h_CrossSection    = new TH1F("CrossSection",    "", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_CrossSectionFit = new TH1F("CrossSectionFit", "", 
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
	
	double locMaxCS = 1.65;
	for(int ibin=0; ibin<m_angularBin.size(); ibin++) {
		double locYield    = h_Yield->GetBinContent(ibin+1);
		double locYieldErr = h_Yield->GetBinError(ibin+1);
		
		double locYieldFit    = h_YieldFit->GetBinContent(ibin+1);
		double locYieldFitErr = h_YieldFit->GetBinError(ibin+1);
		
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
		
		if(0) {
			double locTheta = h_CrossSection->GetBinCenter(ibin+1);
			locCS          /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			locCS_error    /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			
			locCSFit       /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
			locCSFit_error /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
		}
		
		h_CrossSection->SetBinContent(ibin+1, locCS);
		h_CrossSection->SetBinError(ibin+1, locCS_error);
		
		h_CrossSectionFit->SetBinContent(ibin+1, locCSFit);
		h_CrossSectionFit->SetBinError(ibin+1, locCSFit_error);
		
		if(locCS>locMaxCS) {
			locMaxCS = locCS;
		}
		if(locCSFit>locMaxCS) {
			locMaxCS = locCSFit;
		}
	}
	
	h_CrossSection->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxCS);
	h_CrossSectionFit->GetYaxis()->SetRangeUser(0.0, 1.2*locMaxCS);
	
	if(cCrossSection==NULL) {
		cCrossSection = new TCanvas("cCrossSection", "Cross Section", 950, 700);
		styleCanvas(cCrossSection);
	}
	
	cCrossSection->cd();
	h_CrossSection->Draw("PE1X0");
	h_CrossSectionFit->Draw("PE1X0 same");
	cCrossSection->Update();
	cCrossSection->Modified();
	
	return;
}

void EtaAnalyzer::InitializeBinning()
{
	double locAngle = m_minReconAngle + 0.5*m_reconAngleBinSize;
	
	while(locAngle < m_maxReconAngle) 
	{
		m_angularBin.push_back({locAngle, 0.5*m_reconAngleBinSize});
		m_angularYield.push_back({0.0, 0.0});
		m_angularYieldFit.push_back({0.0, 0.0});
		m_angularYieldEmpty.push_back({0.0, 0.0});
		m_angularSBR.push_back({0.0, 0.0});
		
		locAngle += m_reconAngleBinSize;
	}
	m_binningSet = true;
	return;
}

void InitializeMggFitter(MggFitter &fitter, EtaAnalyzer *anaObj, double angle)
{
	fitter.angle            = angle;
	fitter.binSize          = anaObj->GetMggBinSize();
	fitter.fitOption_signal = anaObj->GetFitOption(1);
	fitter.fitOption_bkgd   = anaObj->GetFitOption(2);
	fitter.fitOption_poly   = anaObj->GetFitOption(3);
	fitter.fitOption_omega  = anaObj->GetFitOption(4);
	fitter.fitOption_etap   = anaObj->GetFitOption(5);
	anaObj->GetFitRange(fitter.minFitRange, fitter.maxFitRange);
	
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
	if(h_Yield) h_Yield->Write();
	if(h_YieldFit) h_YieldFit->Write();
	if(h_CrossSection) h_CrossSection->Write();
	if(h_CrossSectionFit) h_CrossSectionFit->Write();
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
	printf("    Background: %s\n", GetFitOptionStr(1).Data());
	printf("    Omega: %s\n", GetFitOptionStr(2).Data());
	printf("    Eta prime: %s\n", GetFitOptionStr(3).Data());
	printf("    Fitting range: %.3f GeV - %.3f GeV\n", m_minFitRange, m_maxFitRange);
	printf("\n=================================================\n");
	printf("\n\n");
	
	return;
}
