#include "CrossSection.h"
#include "YieldFitter.h"

int InitializeYieldFitter(YieldFitter &fitter, EtaAnalyzer &anaObj)
{
	double locBinSize, locMin, locMax;
	
	anaObj.GetBeamEnergyBinning(locBinSize, locMin, locMax);
	fitter.SetBeamEnergyBinning(locBinSize, locMin, locMax);
	
	anaObj.GetReconAngleBinning(locBinSize, locMin, locMax);
	fitter.SetReconAngleBinning(locBinSize, locMin, locMax);
	
	anaObj.GetThrownAngleBinning(locBinSize, locMin, locMax);
	fitter.SetThrownAngleBinning(locBinSize, locMin, locMax);
	
	fitter.SetLuminosity(anaObj.GetLuminosity());
	
	// check that the angular matrices have been loaded. If not, load them:
	if(!anaObj.IsMatrixLoaded()) {
		if(anaObj.LoadAngularMatrix()) {
			printf("\n\nProblem loading angular matrices.\n\n");
			return 1;
		}
	}
	
	fitter.SetAngularMatrix((TH3F*)anaObj.GetAngularMatrix());
	fitter.SetAngularMatrixFine((TH3F*)anaObj.GetAngularMatrixFine());
	fitter.SetFluxWeights((TH1F*)anaObj.GetFluxWeights());
	
	fitter.SetYield((TH1F*)anaObj.GetAngularYield(1));
	
	return 0;
}

void YieldFitter::FitAngularYield(double minFitRange, double maxFitRange)
{
	if(LoadTheoryHists()) return;
	
	printf("Fitting from: %.3f GeV - %.3f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
	printf("Flux Weights Integral: %f\n", h_fluxWeights->Integral());
	
	InitializeFitFunction(&f_yield, "f_yield");
	f_yield->SetParameters(0.515, 1.0, 60.0, 1.0, 1.0);
	f_yield->SetParLimits(0,    0.1,   2.5);
	f_yield->SetParLimits(1,    0.1,   2.0);
	f_yield->SetParLimits(2,  -90.0, 360.0);
	f_yield->SetParLimits(3,    0.0,   1.5);
	if(m_components.size()>3) {
		f_yield->SetParLimits(4, 0.1, 5.0);
	}
	
	if(0) {
		// Fix incoherent normalization in fit:
		f_yield->FixParameter(3, 0.552842);
		if(m_components.size()>3) {
			f_yield->FixParameter(4, 0.6);
		}
	}
	f_yield->SetParameter(2, 90.0);
	
	f_yield->SetRange(minFitRange, maxFitRange);
	
	h_yield->Fit(f_yield, "R0");
	
	TF1 *fitresult = h_yield->GetFunction("f_yield");
	printf("Fit Result:\n");
	printf("  Chi-squared: %f\n", fitresult->GetChisquare());
	printf("  NDF:         %d\n", fitresult->GetNDF());
	
	ofstream outf("results.txt");
	char buf[256];
	sprintf(buf, "%f  %f\n", f_yield->GetParameter(0), f_yield->GetParError(0));
	outf << buf;
	sprintf(buf, "%f  %f\n", f_yield->GetParameter(2), f_yield->GetParError(2));
	outf << buf;
	sprintf(buf, "%f  %f\n", f_yield->GetParameter(1), f_yield->GetParError(1));
	outf << buf;
	sprintf(buf, "%f  %f\n", f_yield->GetParameter(3), f_yield->GetParError(3));
	outf << buf;
	sprintf(buf, "%f  %d\n", fitresult->GetChisquare(), fitresult->GetNDF());
	outf << buf;
	outf.close();
	
	DrawFitResult(0.0, maxFitRange);
	
	return;
}

void YieldFitter::CombinedFit(double minFitRange, double maxFitRange)
{
	if(LoadTheoryHists()) return;
	
	if(h_yield_8GeV==nullptr) return;
	if(h_yield_9GeV==nullptr) return;
	if(h_yield_10GeV==nullptr) return;
	
	//=======================================================================================================//
	// Define the struct needed to compute global chi2 from combined fit:
	
	struct GlobalChi2 {
		const ROOT::Math::IMultiGenFunction *chiSq_8GeV;
		const ROOT::Math::IMultiGenFunction *chiSq_9GeV;
		const ROOT::Math::IMultiGenFunction *chiSq_10GeV;
		
		vector<int> index_8GeV, index_9GeV, index_10GeV;
		
		GlobalChi2(
			const ROOT::Math::IMultiGenFunction &c1,
			const ROOT::Math::IMultiGenFunction &c2,
			const ROOT::Math::IMultiGenFunction &c3,
			const vector<int> &index1, 
			const vector<int> &index2,
			const vector<int> &index3
		) : chiSq_8GeV(&c1), chiSq_9GeV(&c2), chiSq_10GeV(&c3), index_8GeV(index1), index_9GeV(index2), index_10GeV(index3) {}
		GlobalChi2() {};
		
		double operator()(const double* p) const {
			vector<double> p1(index_8GeV.size()), p2(index_9GeV.size()), p3(index_10GeV.size());
			
			for(size_t ipar=0; ipar<index_8GeV.size(); ipar++) p1[ipar] = p[index_8GeV[ipar]];
			for(size_t ipar=0; ipar<index_9GeV.size(); ipar++) p2[ipar] = p[index_9GeV[ipar]];
			for(size_t ipar=0; ipar<index_10GeV.size(); ipar++) p3[ipar] = p[index_10GeV[ipar]];
			return (*chiSq_8GeV)(p1.data()) + (*chiSq_9GeV)(p2.data()) + (*chiSq_10GeV)(p3.data());
		}
	};
	//=======================================================================================================//
	
	return;
}

void YieldFitter::DrawFitResult(double min, double max)
{
	c_fit = new TCanvas("yield_fit", "Yield Fit", 950, 700);
	styleCanvas(c_fit);
	
	TF1 *fDraw, *fPrim, *fCoh, *fQFP, *fQFN, *fInt;
	
	InitializeDrawFunction(&fDraw, "f_Draw",         kBlack,   2, 3);
	InitializeDrawFunction(&fPrim, "f_Primakoff",    kRed,     2, 2);
	InitializeDrawFunction(&fCoh,  "f_Coherent",     kBlue,    2, 2);
	InitializeInterFunction(&fInt, "f_Interference", kMagenta, 2, 2);
	if(m_components.size()==4) {
		InitializeDrawFunction(&fQFP, "f_QuasifreeProton",  kGreen,   2, 2);
		InitializeDrawFunction(&fQFN, "f_QuasifreeNeutron", kGreen-7, 2, 2);
	}
	else {
		InitializeDrawFunction(&fQFP, "f_Incoherent", kGreen, 2, 2);
	}
	
	fDraw->SetParameters(f_yield->GetParameters());
	fPrim->SetParameters(f_yield->GetParameter(0), 0.0, 0.0, 0.0, 0.0);
	fCoh->SetParameters(0.0, f_yield->GetParameter(1), 0.0, 0.0, 0.0);
	fInt->SetParameters(f_yield->GetParameter(0), f_yield->GetParameter(1), f_yield->GetParameter(2));
	
	fQFP->SetParameters(0.0, 0.0, 0.0, f_yield->GetParameter(3), 0.0);
	if(m_components.size()>3) fQFN->SetParameters(0.0, 0.0, 0.0, 0.0, f_yield->GetParameter(4));
	
	h_yield->GetXaxis()->SetRangeUser(min, max);
	h_yield->SetMinimum(0.0);
	
	//h_yield->SetMarkerStyle(1);
	
	c_fit->cd();
	h_yield->Draw("PE1X0");
	if(m_components.size()>3) fQFN->Draw("same");
	fQFP->Draw("same");
	fInt->Draw("same");
	fCoh->Draw("same");
	fPrim->Draw("same");
	fDraw->Draw("same");
	
	TLegend *locLeg = new TLegend(0.154, 0.657, 0.454, 0.907);
	locLeg->AddEntry(fDraw, "Total Fit", "l");
	locLeg->AddEntry(fPrim, "Primakoff", "l");
	locLeg->AddEntry(fCoh, "Nuclear Coherent", "l");
	locLeg->AddEntry(fInt, "Interference", "l");
	locLeg->AddEntry(fQFP, "Nuclear Incoherent", "l");
	//locLeg->Draw();
	
	
	
	
	// Draw the dN/dOmega distribution too:
	
	h_yield_dOmega = (TH1F*)h_yield->Clone("h_yield_dOmega");
	double locMax_dOmega = 0.0;
	
	for(int ibin=1; ibin<=h_yield->GetXaxis()->GetNbins(); ibin++) {
		double locYield    = h_yield->GetBinContent(ibin);
		double locYieldErr = h_yield->GetBinError(ibin);
		double locAngle    = h_yield->GetXaxis()->GetBinCenter(ibin);
		double locSinTh    = 1.0/(2.0*TMath::Pi()*sin(locAngle*TMath::DegToRad()));
		
		h_yield_dOmega->SetBinContent(ibin, locSinTh*locYield);
		h_yield_dOmega->SetBinError(ibin, locSinTh*locYieldErr);
		
		if((locSinTh*locYield)>locMax_dOmega) locMax_dOmega = locSinTh*locYield;
	}
	locMax_dOmega += h_yield_dOmega->GetBinError(1);
	
	h_yield_dOmega->GetYaxis()->SetTitle("dN(#eta#rightarrow#gamma#gamma) / d#Omega");
	h_yield_dOmega->GetYaxis()->SetRangeUser(0.0, 1.2*locMax_dOmega);
	h_yield_dOmega->GetYaxis()->SetMaxDigits(3);
	
	c_fit_dOmega = new TCanvas("yield_fit_dO", "dN/dOmega fit", 950, 700);
	styleCanvas(c_fit_dOmega);
	
	TF1 *fDraw_dO, *fPrim_dO, *fCoh_dO, *fQFP_dO, *fQFN_dO, *fInt_dO;
	
	InitializeDrawFunction_dOmega(&fDraw_dO, "f_Draw_dO",         kBlack,   2, 3);
	InitializeDrawFunction_dOmega(&fPrim_dO, "f_Primakoff_dO",    kRed,     2, 2);
	InitializeDrawFunction_dOmega(&fCoh_dO,  "f_Coherent_dO",     kBlue,    2, 2);
	InitializeInterFunction_dOmega(&fInt_dO, "f_Interference_dO", kMagenta, 2, 2);
	if(m_components.size()==4) {
		InitializeDrawFunction_dOmega(&fQFP_dO, "f_QuasifreeProton_dO",  kGreen,   2, 2);
		InitializeDrawFunction_dOmega(&fQFN_dO, "f_QuasifreeNeutron_dO", kGreen-7, 2, 2);
	}
	else {
		InitializeDrawFunction_dOmega(&fQFP_dO, "f_Incoherent_dO", kGreen, 2, 2);
	}
	
	fDraw_dO->SetParameters(f_yield->GetParameters());
	fPrim_dO->SetParameters(f_yield->GetParameter(0), 0.0, 0.0, 0.0, 0.0);
	fCoh_dO->SetParameters(0.0, f_yield->GetParameter(1), 0.0, 0.0, 0.0);
	fInt_dO->SetParameters(f_yield->GetParameter(0), f_yield->GetParameter(1), f_yield->GetParameter(2));
	
	fQFP_dO->SetParameters(0.0, 0.0, 0.0, f_yield->GetParameter(3), 0.0);
	if(m_components.size()>3) fQFN_dO->SetParameters(0.0, 0.0, 0.0, 0.0, f_yield->GetParameter(4));
	
	h_yield_dOmega->GetXaxis()->SetRangeUser(min, max);
	h_yield_dOmega->SetMinimum(0.0);
	
	
	//h_yield_dOmega->SetMarkerStyle(1);
	
	c_fit_dOmega->cd();
	h_yield_dOmega->Draw("PE1X0");
	if(m_components.size()>3) fQFN_dO->Draw("same");
	fQFP_dO->Draw("same");
	fInt_dO->Draw("same");
	fCoh_dO->Draw("same");
	fPrim_dO->Draw("same");
	fDraw_dO->Draw("same");
	
	TLegend *locLeg_dO = new TLegend(0.597, 0.648, 0.898, 0.898);
	locLeg_dO->AddEntry(fDraw_dO, "Total Fit", "l");
	locLeg_dO->AddEntry(fPrim_dO, "Primakoff", "l");
	locLeg_dO->AddEntry(fCoh_dO, "Nuclear Coherent", "l");
	locLeg_dO->AddEntry(fInt_dO, "Interference", "l");
	locLeg_dO->AddEntry(fQFP_dO, "Nuclear Incoherent", "l");
	locLeg_dO->Draw();
	
	return;
}

void YieldFitter::InitializeFitFunction(TF1 **f1, TString funcName, int lineColor)
{
	// initialize fit function for each angular bin:
	
	int nParameters = m_components.size() + 1;;
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldFitFunction, m_minReconAngle, m_maxReconAngle, nParameters);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "#phi[#circ]");
	if(m_components.size()==4) {
		(*f1)->SetParName(3, "A_{QFP}");
		(*f1)->SetParName(4, "A_{QFN}");
	}
	else {
		(*f1)->SetParName(3, "A_{Inc}");
	}
	
	(*f1)->SetLineColor(lineColor);
	
	return;
}

void YieldFitter::InitializeDrawFunction(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	int nParameters = m_components.size() + 1;
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldDrawFunction, m_minReconAngle, m_maxReconAngle, nParameters);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "#phi[#circ]");
	if(m_components.size()==4) {
		(*f1)->SetParName(3, "A_{QFP}");
		(*f1)->SetParName(4, "A_{QFN}");
	}
	else {
		(*f1)->SetParName(3, "A_{Inc}");
	}
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

void YieldFitter::InitializeDrawFunction_dOmega(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	int nParameters = m_components.size() + 1;
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldDrawFunction_dOmega, m_minReconAngle, m_maxReconAngle, nParameters);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "#phi[#circ]");
	if(m_components.size()==4) {
		(*f1)->SetParName(3, "A_{QFP}");
		(*f1)->SetParName(4, "A_{QFN}");
	}
	else {
		(*f1)->SetParName(3, "A_{Inc}");
	}
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

void YieldFitter::InitializeInterFunction(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldDrawFunctionInterference, m_minReconAngle, m_maxReconAngle, 3);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "#phi[#circ]");
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

void YieldFitter::InitializeInterFunction_dOmega(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldDrawFunctionInterference_dOmega, m_minReconAngle, m_maxReconAngle, 3);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "#phi[#circ]");
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

int YieldFitter::LoadTheoryHists()
{
	printf("\nREADING THEORY CALCULATIONS...\n");
	
	int nComponents = InitializeModelComponents();
	
	vector<TString> theoryFileNames;
	for(int i=0; i<nComponents; i++) {
		theoryFileNames.push_back("");
	}
	
	switch(m_model) {
		case AFIX:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			break;
		}
		case SGEVORKYAN:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			break;
		}
		case MIXED_V1:
		{
			for(int i=1; i<2; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			for(int i=2; i<4; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			break;
		}
		case MIXED_V2:
		{
			theoryFileNames[0] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			for(int i=1; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			break;
		}
		case SGEVORKYAN_FERMI:
		{
			for(int i=0; i<(nComponents-1); i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			theoryFileNames[nComponents-1] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/test-sergey-gevorkyanb.root";
			break;
		}
		case SGEVORKYAN_UPD_V0:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v0.root";
			break;
		}
		case SGEVORKYAN_UPD_V1:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v1.root";
			break;
		}
		case SGEVORKYAN_UPD_V2:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root";
			break;
		}
		case SGEVORKYAN_UPD_V3:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v3.root";
			break;
		}
		case SGEVORKYAN_UPD_FERMI:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-sgevorkyan-v2-folded.root";
			break;
		}
	}
	
	for(int i=0; i<nComponents; i++) {
		
		TFile *fTheory  = new TFile(theoryFileNames[i].Data(), "READ");
		TH2F *hTheory2D = (TH2F*)fTheory->Get(Form("%s",m_components[i].Data()));
		
		double locBeamEnergy = m_minBeamEnergy + 0.5*m_beamEnergyBinSize;
		while(locBeamEnergy < m_maxBeamEnergy) {
			int locMinEnergyBin = hTheory2D->GetXaxis()->FindBin(locBeamEnergy-0.5*m_beamEnergyBinSize);
			int locMaxEnergyBin = hTheory2D->GetXaxis()->FindBin(locBeamEnergy+0.5*m_beamEnergyBinSize);
			
			double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin);
			
			TH1F *loch1 = (TH1F*)hTheory2D->ProjectionY(Form("%s_%f_%f", m_components[i].Data(), 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
			loch1->Scale(1.0/nEnergyBins);
			
			// rebin histograms to match acceptance matrices:
			
			double theoryAngleBinSize = loch1->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				loch1->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				loch1->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			loch1->SetDirectory(0);
			h_Theory[i].push_back(loch1);
			locBeamEnergy += m_beamEnergyBinSize;
		}
		
		fTheory->Close();
	}
	
	if(m_model>=5) {
		
		TFile *fTheory = new TFile(theoryFileNames[0].Data(), "READ");
		
		TH2F *hPrimReal2D = (TH2F*)fTheory->Get("amp_prim_real_vs_egam");
		TH2F *hPrimImag2D = (TH2F*)fTheory->Get("amp_prim_imag_vs_egam");
		TH2F *hStrongReal2D = (TH2F*)fTheory->Get("amp_coh_real_vs_egam");
		TH2F *hStrongImag2D = (TH2F*)fTheory->Get("amp_coh_imag_vs_egam");
		
		double locBeamEnergy = m_minBeamEnergy + 0.5*m_beamEnergyBinSize;
		while(locBeamEnergy < m_maxBeamEnergy) {
			int locMinEnergyBin = hPrimReal2D->GetXaxis()->FindBin(locBeamEnergy-0.5*m_beamEnergyBinSize);
			int locMaxEnergyBin = hPrimReal2D->GetXaxis()->FindBin(locBeamEnergy+0.5*m_beamEnergyBinSize);
			
			double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin);
			
			TH1F *locPrimReal = (TH1F*)hPrimReal2D->ProjectionY(Form("prim_real_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
			locPrimReal->Scale(1.0/nEnergyBins);
			
			TH1F *locPrimImag = (TH1F*)hPrimImag2D->ProjectionY(Form("prim_imag_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
			locPrimImag->Scale(1.0/nEnergyBins);
			
			TH1F *locStrongReal = (TH1F*)hStrongReal2D->ProjectionY(Form("coh_real_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
			locStrongReal->Scale(1.0/nEnergyBins);
			
			TH1F *locStrongImag = (TH1F*)hStrongImag2D->ProjectionY(Form("coh_imag_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
			locStrongImag->Scale(1.0/nEnergyBins);
			
			// rebin histograms to match acceptance matrices:
			
			double theoryAngleBinSize;
			
			theoryAngleBinSize = locPrimReal->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locPrimReal->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locPrimReal->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			locPrimReal->SetDirectory(0);
			h_PrimReal.push_back(locPrimReal);
			
			theoryAngleBinSize = locPrimImag->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locPrimImag->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locPrimImag->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			locPrimImag->SetDirectory(0);
			h_PrimImag.push_back(locPrimImag);
			
			theoryAngleBinSize = locStrongReal->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locStrongReal->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locStrongReal->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			locStrongReal->SetDirectory(0);
			h_StrongReal.push_back(locStrongReal);
			
			theoryAngleBinSize = locStrongImag->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locStrongImag->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locStrongImag->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			locStrongImag->SetDirectory(0);
			h_StrongImag.push_back(locStrongImag);
			
			locBeamEnergy += m_beamEnergyBinSize;
		}
		
		fTheory->Close();
	}
	
	h_TheoryTulio = new TH1F("TulioInc", "; #theta_{#eta} [#circ]; d#sigma/d#theta_{#eta} [#mub/rad]", 51, 0.0, 5.1);
	
	ifstream inf("/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/NI_ETA_4He_MCMC_Preliminary.dat");
	double aa, bb, cc;
	for(int iline=0; iline<51; iline++) {
		inf >> aa >> bb >> cc;
		int locBin = h_TheoryTulio->FindBin(aa);
		h_TheoryTulio->SetBinContent(locBin, cc);
	}
	inf.close();
	
	
	return 0;
}

int YieldFitter::InitializeModelComponents() 
{
	int nComponents = 0;
	switch(m_model) {
		case AFIX:
			nComponents = 4;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_qfp_vs_egam");
			m_components.push_back("xs_qfn_vs_egam");
			break;
		case SGEVORKYAN:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case MIXED_V1:
			nComponents = 4;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_qfp_vs_egam");
			m_components.push_back("xs_qfn_vs_egam");
			break;
		case MIXED_V2:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_FERMI:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("qf_txs_lab");
			break;
		case SGEVORKYAN_UPD_V0:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V1:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V2:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V3:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_FERMI:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("qf_txs_lab");
			break;
	}
	
	for(int i=0; i<nComponents; i++) {
		vector<TH1F*> locVec;
		h_Theory.push_back(locVec);
	}
	
	return nComponents;
}

double YieldFitter::YieldFitFunction(double *x, double *par)
{
	if(h_matrix==NULL) return 0.0;
	
	double reconAngle = x[0];
	
	if(reconAngle>0.08 && reconAngle<0.12) {
		//TF1::RejectPoint();
		//return 0.;
	}
	
	int reconAngleBin = h_matrix->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double Phi   = par[2];
	double AincP = par[3];
	double AincN = (m_components.size()>3) ? par[4] : 0.0;
	
	int locMinEnergyBin = h_matrix->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrix->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrix->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrix->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrix->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSection(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, AincP, AincN, Phi, locEnergy);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * EtaAnalyzer::m_branchingRatio * angleBinSize;
	return yield;
}

double YieldFitter::YieldDrawFunction(double *x, double *par)
{
	if(h_matrixFine==NULL) return 0.0;
	
	double reconAngle = x[0];
	int reconAngleBin = h_matrixFine->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	
	double scaleFactor = h_matrix->GetYaxis()->GetBinWidth(1) / h_matrixFine->GetYaxis()->GetBinWidth(1);
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double Phi   = par[2];
	double AincP = par[3];
	double AincN = (m_components.size()>3) ? par[4] : 0.0;
	
	int locMinEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrixFine->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSection(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, AincP, AincN, Phi, locEnergy);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * EtaAnalyzer::m_branchingRatio * angleBinSize * scaleFactor;
	return yield;
}

double YieldFitter::YieldDrawFunctionInterference(double *x, double *par)
{
	if(h_matrixFine==NULL) return 0.0;
	
	double reconAngle = x[0];
	int reconAngleBin = h_matrixFine->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	double scaleFactor  = h_matrix->GetYaxis()->GetBinWidth(1) / h_matrixFine->GetYaxis()->GetBinWidth(1);
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double Phi   = par[2];
	
	int locMinEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrixFine->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSectionInterference(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, Phi, locEnergy);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * EtaAnalyzer::m_branchingRatio * angleBinSize * scaleFactor;
	return yield;
}

double YieldFitter::YieldDrawFunction_dOmega(double *x, double *par)
{
	if(h_matrixFine==NULL) return 0.0;
	
	double reconAngle = x[0];
	int reconAngleBin = h_matrixFine->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	
	double scaleFactor = h_matrix->GetYaxis()->GetBinWidth(1) / h_matrixFine->GetYaxis()->GetBinWidth(1);
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double Phi   = par[2];
	double AincP = par[3];
	double AincN = (m_components.size()>3) ? par[4] : 0.0;
	
	int locMinEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrixFine->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSection(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, AincP, AincN, Phi, locEnergy);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * EtaAnalyzer::m_branchingRatio * angleBinSize * scaleFactor;
	double sinth = 2.0*TMath::Pi()*sin(reconAngle*TMath::DegToRad());
	yield /= sinth;
	return yield;
}

double YieldFitter::YieldDrawFunctionInterference_dOmega(double *x, double *par)
{
	if(h_matrixFine==NULL) return 0.0;
	
	double reconAngle = x[0];
	int reconAngleBin = h_matrixFine->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	double scaleFactor  = h_matrix->GetYaxis()->GetBinWidth(1) / h_matrixFine->GetYaxis()->GetBinWidth(1);
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double Phi   = par[2];
	
	int locMinEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrixFine->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSectionInterference(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, Phi, locEnergy);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * EtaAnalyzer::m_branchingRatio * angleBinSize * scaleFactor;
	double sinth = 2.0*TMath::Pi()*sin(reconAngle*TMath::DegToRad());
	yield /= sinth;
	return yield;
}

double YieldFitter::GetCrossSection(int energyBin, int thetaBin, 
	double Gamma, double Acoh, double AincP, double AincN, double Phi, double beamEnergy) 
{
	double gammaGen = (m_model>=5) ? 0.515 : 0.510;
	
	double locPrim = (Gamma/gammaGen) * h_Theory[0][energyBin]->GetBinContent(thetaBin);
	double locCoh  = Acoh * h_Theory[1][energyBin]->GetBinContent(thetaBin);
	double locInc  = AincP * h_Theory[2][energyBin]->GetBinContent(thetaBin);
	if(m_components.size()>3) {
		locInc += (AincN * h_Theory[3][energyBin]->GetBinContent(thetaBin));
	}
	
	// use tulio's incoherent:
	//int locTulioBin = h_TheoryTulio->FindBin(h_Theory[0][energyBin]->GetXaxis()->GetBinCenter(thetaBin));
	//locInc = AincP * h_TheoryTulio->GetBinContent(locTulioBin);
	
	// TEST:
	//locCoh /= pow(beamEnergy,0.2);
	
	double locInt  = GetCrossSectionInterference(energyBin, thetaBin, Gamma, Acoh, Phi, beamEnergy);
	return (locPrim+locCoh+locInc+locInt);
}

double YieldFitter::GetCrossSectionInterference(int energyBin, int thetaBin, 
	double Gamma, double Acoh, double Phi, double beamEnergy) 
{
	double gammaGen = (m_model>=5) ? 0.515 : 0.510;
	
	double locInt = 0.0;
	if(m_model>=5) {
		double locPrimReal = h_PrimReal[energyBin]->GetBinContent(thetaBin);
		double locPrimImag = h_PrimImag[energyBin]->GetBinContent(thetaBin);
		double locStrongReal = h_StrongReal[energyBin]->GetBinContent(thetaBin);
		double locStrongImag = h_StrongImag[energyBin]->GetBinContent(thetaBin);
		
		// TEST:
		//locStrongReal /= pow(beamEnergy,0.2);
		//locStrongImag /= pow(beamEnergy,0.2);
		
		double int1 = (locPrimReal*locStrongReal + locPrimImag*locStrongImag) * cos(Phi*TMath::DegToRad());
		double int2 = (locPrimImag*locStrongReal - locPrimReal*locStrongImag) * sin(Phi*TMath::DegToRad());
		locInt = 2.0*sqrt(Gamma/gammaGen)*sqrt(Acoh)*(int1+int2);
	} else {
		double locPrim = (Gamma/gammaGen) * h_Theory[0][energyBin]->GetBinContent(thetaBin);
		double locCoh  = Acoh * h_Theory[1][energyBin]->GetBinContent(thetaBin);
		
		// TEST:
		//locCoh /= pow(beamEnergy,0.4);
		
		locInt  = 2.0*sqrt(locPrim*locCoh)*cos(Phi*TMath::DegToRad());
	}
	return locInt;
}
