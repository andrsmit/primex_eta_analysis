#include "EtaAnalyzer.h"
#include "CrossSection.h"

EtaAnalyzer::EtaAnalyzer(int phase) 
{
	if((phase<1) || (phase>3)) {
		cout << "\nUnsupported PrimEx-eta phase number provided.\n" << endl;
		exit(1);
	}
	m_phase = phase;
	
	h_etaLineshape    = NULL;
	h_omegaLineshape  = NULL;
	
	h_mggVsThetaFull  = NULL;
	h_mggVsThetaEmpty = NULL;
	
	cFit              = NULL;
	cYield            = NULL;
	cCrossSection     = NULL;
	cAcceptance       = NULL;
	
	l0  = NULL;
	lp  = NULL;
	lm  = NULL;
	lx1 = NULL;
	lx2 = NULL;
}

void EtaAnalyzer::InitializeFitCanvas()
{
	cFit = new TCanvas("cFit", "Mgg Fit", 700, 500);
	
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

//---------------------------------------------------------//

void EtaAnalyzer::SetRebinsMgg(int rebins) 
{
	m_rebinsMgg  = rebins;
	m_mggBinSize = (1.e-3) * (double)m_rebinsMgg;
	return;
}
void EtaAnalyzer::SetRebinsTheta(int rebins) 
{
	m_rebinsTheta  = rebins;
	m_thetaBinSize = 0.01 * (double)m_rebinsTheta;
	return;
}

//---------------------------------------------------------//

void EtaAnalyzer::SetFitOption_signal(int option) 
{
	if((option<1) || (option>6)) {
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

void EtaAnalyzer::SetFitRange(double min, double max)
{
	m_minFitRange = min;
	m_maxFitRange = max;
	return;
}

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

void EtaAnalyzer::GetFitRange(double &min, double &max)
{
	min = m_minFitRange;
	max = m_maxFitRange;
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
	
	locHistFull->Add(locHistEmpty, -1.0);
	
	MggFitter locFitter;
	InitializeFitterSettings(locFitter, this);
	
	// Load 1-D lineshape projections:
	
	if(h_etaLineshape) {
		minAngleBin = h_etaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_etaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
		hEta->Rebin(m_rebinsMgg);
		locFitter.SetEtaLineshape(hEta);
	}
	
	if(h_omegaLineshape) {
		minAngleBin = h_omegaLineshape->GetXaxis()->FindBin(minAngle);
		maxAngleBin = h_omegaLineshape->GetXaxis()->FindBin(maxAngle)-1;
		TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega", minAngleBin, maxAngleBin);
		hOmega->Rebin(m_rebinsMgg);
		locFitter.SetOmegaLineshape(hOmega, drawOmegaFit);
	}
	
	locFitter.SetData(locHistFull);
	locFitter.FitData();
	
	double yield, yieldErr;
	locFitter.GetYield(yield, yieldErr);
	printf("\nYield = %f +/- %f\n", yield, yieldErr);
	
	if(drawFitResult)
	{
		TH1F *locPull = (TH1F*)locHistFull->Clone("lochPull");
		locFitter.FillPull(locPull);
		TF1 *locfFit    = (TF1*)locFitter.GetFitFunction()->Clone("locfFit");
		TF1 *locfSignal = (TF1*)locFitter.GetSignalFunction()->Clone("locfSignal");
		TF1 *locfBkgd   = (TF1*)locFitter.GetBkgdFunction()->Clone("locfBkgd");
		DrawFitResult(locHistFull, locPull, locfFit, locfSignal, locfBkgd);
	}
	return;
}

void EtaAnalyzer::DrawFitResult(TH1F *h1, TH1F *hPull, TF1 *fFit, TF1 *fSignal, TF1 *fBkgd)
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
	
	pFit->Update();
	
	l0->Draw("same");
	lx1->SetY1(gPad->GetUymin());
	lx1->SetY2(gPad->GetUymax());
	lx1->Draw("same");
	lx2->SetY1(gPad->GetUymin());
	lx2->SetY2(gPad->GetUymax());
	lx2->Draw("same");
	
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
		
		//----------------------------------------------//
		// Get 1-d projection of invariant mass spectrum:
		
		int minAngleBin = m_rebinsTheta*(iThetaBin);
		int maxAngleBin = m_rebinsTheta*(iThetaBin+1);
		
		TH1F *locHistFull  = (TH1F*)h_mggVsThetaFull->ProjectionY("locHistFull", minAngleBin, maxAngleBin);
		locHistFull->Rebin(m_rebinsMgg);
		styleMggHistogram(locHistFull);
		
		TH1F *locHistEmpty = (TH1F*)h_mggVsThetaEmpty->ProjectionY("locHistEmpty", minAngleBin, maxAngleBin);
		locHistEmpty->Rebin(m_rebinsMgg);
		styleMggHistogram(locHistFull);
		
		locHistFull->Add(locHistEmpty,-1.0);
		
		//----------------------------------------------//
		// Set up fitter object:
		
		MggFitter locFitter;
		InitializeFitterSettings(locFitter, this);
		
		locFitter.SetData(locHistFull);
		if(h_etaLineshape) {
			minAngleBin = h_etaLineshape->GetXaxis()->FindBin(m_angularBin[iThetaBin].first - m_angularBin[iThetaBin].second);
			maxAngleBin = h_etaLineshape->GetXaxis()->FindBin(m_angularBin[iThetaBin].first + m_angularBin[iThetaBin].second)-1;
			TH1F *hEta = (TH1F*)h_etaLineshape->ProjectionY("hEta", minAngleBin, maxAngleBin);
			locFitter.SetEtaLineshape(hEta);
		}
		if(h_omegaLineshape) {
			minAngleBin = h_omegaLineshape->GetXaxis()->FindBin(m_angularBin[iThetaBin].first - m_angularBin[iThetaBin].second);
			maxAngleBin = h_omegaLineshape->GetXaxis()->FindBin(m_angularBin[iThetaBin].first + m_angularBin[iThetaBin].second)-1;
			TH1F *hOmega = (TH1F*)h_omegaLineshape->ProjectionY("hOmega");//, minAngleBin, maxAngleBin);
			locFitter.SetOmegaLineshape(hOmega);
		}
		
		//----------------------------------------------//
		// Do fit and extract yield:
		
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
			DrawFitResult(locHistFull, locPull, locfFit, locfSignal, locfBkgd);
		}
		printf(" angle, yield1,  yield2 = %f,  %f,  %f\n", m_angularBin[iThetaBin].first, locYield, locYieldFit);
	}
	
	return;
}

void EtaAnalyzer::PlotAngularYield()
{
	h_Yield    = new TH1F("AngularYield",    "", m_angularBin.size(), 0.0, m_thetaBinSize*(double)(m_angularBin.size()));
	h_YieldFit = new TH1F("AngularYieldFit", "", m_angularBin.size(), 0.0, m_thetaBinSize*(double)(m_angularBin.size()));
	
	h_Yield->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_Yield->GetXaxis()->SetTitleSize(0.05);
	h_Yield->GetXaxis()->SetTitleOffset(1.0);
	h_Yield->GetXaxis()->CenterTitle("");
	h_Yield->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_thetaBinSize));
	h_Yield->GetYaxis()->SetTitleSize(0.05);
	h_Yield->GetYaxis()->SetTitleOffset(1.0);
	h_Yield->GetYaxis()->CenterTitle("");
	h_Yield->SetTitle("");
	h_Yield->SetMarkerStyle(4);
	h_Yield->SetMarkerSize(0.9);
	h_Yield->SetMarkerColor(kBlue);
	h_Yield->SetLineColor(kBlue);
	h_Yield->SetLineWidth(2);
	
	h_YieldFit->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_YieldFit->GetXaxis()->SetTitleSize(0.05);
	h_YieldFit->GetXaxis()->SetTitleOffset(1.0);
	h_YieldFit->GetXaxis()->CenterTitle("");
	h_YieldFit->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_thetaBinSize));
	h_YieldFit->GetYaxis()->SetTitleSize(0.05);
	h_YieldFit->GetYaxis()->SetTitleOffset(1.0);
	h_YieldFit->GetYaxis()->CenterTitle("");
	h_YieldFit->SetTitle("");
	h_YieldFit->SetMarkerStyle(4);
	h_YieldFit->SetMarkerSize(0.9);
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
	
	h_Yield->GetYaxis()->SetRangeUser(0.0, 1.4*locMaxYield);
	h_YieldFit->GetYaxis()->SetRangeUser(0.0, 1.4*locMaxYield);
	
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
	h_CrossSection    = new TH1F("CrossSection",    "", m_angularBin.size(), 0.0, m_thetaBinSize*(double)(m_angularBin.size()));
	h_CrossSectionFit = new TH1F("CrossSectionFit", "", m_angularBin.size(), 0.0, m_thetaBinSize*(double)(m_angularBin.size()));
	
	h_CrossSection->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_CrossSection->GetXaxis()->SetTitleSize(0.05);
	h_CrossSection->GetXaxis()->SetTitleOffset(1.0);
	h_CrossSection->GetXaxis()->CenterTitle("");
	h_CrossSection->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_thetaBinSize));
	h_CrossSection->GetYaxis()->SetTitleSize(0.05);
	h_CrossSection->GetYaxis()->SetTitleOffset(1.0);
	h_CrossSection->GetYaxis()->CenterTitle("");
	h_CrossSection->SetTitle("");
	h_CrossSection->SetMarkerStyle(4);
	h_CrossSection->SetMarkerSize(0.9);
	h_CrossSection->SetMarkerColor(kBlue);
	h_CrossSection->SetLineColor(kBlue);
	h_CrossSection->SetLineWidth(2);
	
	h_CrossSectionFit->GetXaxis()->SetTitle("Polar Angle, #theta_{#eta} [#circ]");
	h_CrossSectionFit->GetXaxis()->SetTitleSize(0.05);
	h_CrossSectionFit->GetXaxis()->SetTitleOffset(1.0);
	h_CrossSectionFit->GetXaxis()->CenterTitle("");
	h_CrossSectionFit->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_thetaBinSize));
	h_CrossSectionFit->GetYaxis()->SetTitleSize(0.05);
	h_CrossSectionFit->GetYaxis()->SetTitleOffset(1.0);
	h_CrossSectionFit->GetYaxis()->CenterTitle("");
	h_CrossSectionFit->SetTitle("");
	h_CrossSectionFit->SetMarkerStyle(4);
	h_CrossSectionFit->SetMarkerSize(0.9);
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
		double locBinSize = m_thetaBinSize * TMath::DegToRad(); // bin size in rad.
		
		double locCS    = locYield    / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		double locCSErr = locYieldErr / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		
		double locCS_upper = (locYield + locYieldErr) / ((locAcc-locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCS_lower = (locYield - locYieldErr) / ((locAcc+locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCS_error = (locCS_upper - locCS_lower) / 2.0;
		
		//double locTheta = h_CrossSection->GetBinCenter(ibin+1);
		//locCS       /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
		//locCS_error /= (2.0*TMath::Pi()*sin(locTheta*TMath::DegToRad()));
		
		h_CrossSection->SetBinContent(ibin+1, locCS);
		h_CrossSection->SetBinError(ibin+1, locCS_error);
		
		double locCSFit    = locYieldFit    / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFitErr = locYieldFitErr / (locAcc * m_luminosity * m_branchingRatio * locBinSize);
		
		double locCSFit_upper = (locYieldFit + locYieldFitErr) / ((locAcc-locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFit_lower = (locYieldFit - locYieldFitErr) / ((locAcc+locAccErr) * m_luminosity * m_branchingRatio * locBinSize);
		double locCSFit_error = (locCSFit_upper - locCSFit_lower) / 2.0;
		
		h_CrossSectionFit->SetBinContent(ibin+1, locCSFit);
		h_CrossSectionFit->SetBinError(ibin+1, locCSFit_error);
		/*
		if(locCS>locMaxCS) {
			locMaxCS = locCS;
		}
		if(locCSFit>locMaxCS) {
			locMaxCS = locCSFit;
		}
		*/
	}
	
	h_CrossSection->GetYaxis()->SetRangeUser(0.0, 1.4*locMaxCS);
	h_CrossSectionFit->GetYaxis()->SetRangeUser(0.0, 1.4*locMaxCS);
	
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
	double maxAngleFit = 4.5;
	
	int nBinsFit = (int)(1.e2*maxAngleFit/m_rebinsTheta);
	for(int ibin=0; ibin<nBinsFit; ibin++) {
		int locBinLow  = m_rebinsTheta*(ibin);
		int locBinHigh = m_rebinsTheta*(ibin+1);
		double locMinAngle = 0.01*(double)locBinLow;
		double locMaxAngle = 0.01*(double)locBinHigh;
		double locAngle    = 0.5*(locMinAngle+locMaxAngle);
		double locAngleErr = 0.5*(locMaxAngle-locMinAngle);
		
		m_angularBin.push_back({locAngle, locAngleErr});
		m_angularYield.push_back({0.0, 0.0});
		m_angularYieldFit.push_back({0.0, 0.0});
		m_angularYieldEmpty.push_back({0.0, 0.0});
		m_angularSBR.push_back({0.0, 0.0});
	}
	m_binningSet = true;
	return;
}

void InitializeFitterSettings(MggFitter &fitter, EtaAnalyzer *anaObj)
{
	fitter.binSize          = anaObj->GetMggBinSize();
	fitter.fitOption_signal = anaObj->GetFitOption(1);
	fitter.fitOption_bkgd   = anaObj->GetFitOption(2);
	fitter.fitOption_poly   = anaObj->GetFitOption(3);
	fitter.fitOption_omega  = anaObj->GetFitOption(4);
	fitter.fitOption_etap   = anaObj->GetFitOption(5);
	anaObj->GetFitRange(fitter.minFitRange, fitter.maxFitRange);
	fitter.excludeRegions.clear();
	return;
}

void EtaAnalyzer::DumpSettings() 
{
	printf("\n\n");
	printf("=================================================\n");
	printf("Extracting eta->gg angular yield.\n\n");
	printf("  Beam Energy Range: %.2f GeV - %.2f GeV\n", minBeamEnergy, maxBeamEnergy);
	printf("  Angular binning: %.2f degrees\n", m_thetaBinSize);
	printf("  Fit functions:\n");
	printf("    Signal: %s\n", GetFitOptionStr(0).Data());
	printf("    Background: %s\n", GetFitOptionStr(1).Data());
	printf("    Omega: %s\n", GetFitOptionStr(2).Data());
	printf("    Eta prime: %s\n", GetFitOptionStr(3).Data());
	printf("    Fitting range: %.3f GeV - %.3f GeV\n", m_minFitRange, m_minFitRange);
	printf("\n=================================================\n");
	printf("\n\n");
	
	return;
}
