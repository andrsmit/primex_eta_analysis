#include "CrossSection.h"
#include "YieldFitter.h"

int InitializeYieldFitter(YieldFitter &fitter, EtaAnalyzer anaObj)
{
	double locBinSize, locMin, locMax;
	
	anaObj.GetBeamEnergyBinning(locBinSize, locMin, locMax);
	fitter.SetBeamEnergyBinning(locBinSize, locMin, locMax);
	
	anaObj.GetReconAngleBinning(locBinSize, locMin, locMax);
	fitter.SetReconAngleBinning(locBinSize, locMin, locMax);
	
	anaObj.GetThrownAngleBinning(locBinSize, locMin, locMax);
	fitter.SetThrownAngleBinning(locBinSize, locMin, locMax);
	
	fitter.SetLuminosity(anaObj.GetLuminosity());
	fitter.SetAngularMatrix((TH3F*)anaObj.GetAngularMatrix());
	fitter.SetAngularMatrixFine((TH3F*)anaObj.GetAngularMatrixFine());
	fitter.SetFluxWeights((TH1F*)anaObj.GetFluxWeights());
	
	fitter.SetYield((TH1F*)anaObj.GetAngularYield(1));
	
	printf("Loading theory calculations...");
	if(fitter.LoadTheoryHists()) return 1;
	printf("Finished.\n");
	
	return 0;
}

void YieldFitter::FitAngularYield(int drawOption)
{
	InitializeFitFunction(&f_yield, "f_yield");
	f_yield->SetParameters(0.45, 0.5, 0.6, 1.5, 90.0);
	f_yield->SetParLimits(0, 0.1, 2.5);
	f_yield->SetParLimits(1, 0.1, 2.0);
	f_yield->SetParLimits(2, 0.1, 5.0);
	f_yield->SetParLimits(3, 0.1, 5.0);
	f_yield->SetParLimits(4, -180.0, 180.0);
	/*
	f_yield->FixParameter(1,  6.75607e-01);
	f_yield->FixParameter(2,  4.82658e-01);
	f_yield->FixParameter(3,  1.22440e+00);
	f_yield->FixParameter(4, -2.15931e+01);
	*/
	f_yield->SetRange(0.0, 4.0);
	h_yield->Fit(f_yield, "R0");
	
	if(drawOption) 
	{
		DrawFitResult();
	}
	
	return;
}

void YieldFitter::DrawFitResult()
{
	c_fit = new TCanvas("yield_fit", "Yield Fit", 950, 700);
	styleCanvas(c_fit);
	
	TF1 *fDraw, *fPrim, *fCoh, *fQFP, *fQFN, *fInt;
	
	InitializeDrawFunction(&fDraw, "f_Draw",             kBlack,   2, 3);
	InitializeDrawFunction(&fPrim, "f_Primakoff",        kRed,     2, 2);
	InitializeDrawFunction(&fCoh,  "f_Coherent",         kBlue,    2, 2);
	InitializeDrawFunction(&fQFP,  "f_QuasifreeProton",  kGreen+2, 2, 2);
	InitializeDrawFunction(&fQFN,  "f_QuasifreeNeutron", kGreen-7, 2, 2);
	InitializeInterFunction(&fInt, "f_Interference",     kMagenta, 2, 2);
	
	fDraw->SetParameters(f_yield->GetParameters());
	fPrim->SetParameters(f_yield->GetParameter(0), 0.0, 0.0, 0.0, 0.0);
	fCoh->SetParameters(0.0, f_yield->GetParameter(1), 0.0, 0.0, 0.0);
	fQFP->SetParameters(0.0, 0.0, f_yield->GetParameter(2), 0.0, 0.0);
	fQFN->SetParameters(0.0, 0.0, 0.0, f_yield->GetParameter(3), 0.0);
	fInt->SetParameters(f_yield->GetParameter(0), f_yield->GetParameter(1), f_yield->GetParameter(4));
	
	c_fit->cd();
	h_yield->Draw("PE1X0");
	fQFN->Draw("same");
	fQFP->Draw("same");
	fInt->Draw("same");
	fCoh->Draw("same");
	fPrim->Draw("same");
	fDraw->Draw("same");
	
	return;
}

void YieldFitter::InitializeFitFunction(TF1 **f1, TString funcName, int lineColor)
{
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldFitFunction, m_minReconAngle, m_maxReconAngle, 5);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "A_{QFP}");
	(*f1)->SetParName(3, "A_{QFN}");
	(*f1)->SetParName(4, "#phi[#circ]");
	
	(*f1)->SetLineColor(lineColor);
	
	return;
}

void YieldFitter::InitializeDrawFunction(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldDrawFunction, m_minReconAngle, m_maxReconAngle, 5);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(1, "A_{Coh}");
	(*f1)->SetParName(2, "A_{QFP}");
	(*f1)->SetParName(3, "A_{QFN}");
	(*f1)->SetParName(4, "#phi[#circ]");
	
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

int YieldFitter::LoadTheoryHists()
{
	TString theoryFileName = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
	if(gSystem->AccessPathName(theoryFileName.Data())) {
		printf("\nProblem accessing theory file name.\n");
		return 1;
	}
	
	TH2F *hTheory[4];
	TFile *fTheory = new TFile(theoryFileName.Data(), "READ");
	hTheory[0] = (TH2F*)fTheory->Get("xs_prim_vs_egam");
	hTheory[1] = (TH2F*)fTheory->Get("xs_coh_vs_egam");
	hTheory[2] = (TH2F*)fTheory->Get("xs_qfp_vs_egam");
	hTheory[3] = (TH2F*)fTheory->Get("xs_qfn_vs_egam");
	for(int i=0; i<4; i++) hTheory[i]->SetDirectory(0);
	fTheory->Close();
	
	//-------------------------------------------------------------------//
	// Get 1-dimensional cross sections as average over each energy bin:
	
	h_TheoryPrim.clear();
	h_TheoryCoh.clear();
	h_TheoryQFP.clear();
	h_TheoryQFN.clear();
	
	double locBeamEnergy = m_minBeamEnergy + 0.5*m_beamEnergyBinSize;
	while(locBeamEnergy < m_maxBeamEnergy) {
		int locMinEnergyBin = hTheory[0]->GetXaxis()->FindBin(locBeamEnergy-0.5*m_beamEnergyBinSize);
		int locMaxEnergyBin = hTheory[0]->GetXaxis()->FindBin(locBeamEnergy+0.5*m_beamEnergyBinSize);
		
		double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin);
		
		TH1F *h1Prim = (TH1F*)hTheory[0]->ProjectionY(Form("h1Prim_%f_%f",
			locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
		TH1F *h1Coh  = (TH1F*)hTheory[1]->ProjectionY(Form("h1Coh_%f_%f", 
			locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
		TH1F *h1QFP = (TH1F*)hTheory[2]->ProjectionY(Form("h1QFP_%f_%f", 
			locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
		TH1F *h1QFN = (TH1F*)hTheory[3]->ProjectionY(Form("h1QFN_%f_%f", 
			locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin-1);
		
		h1Prim->Scale(1.0/nEnergyBins);
		h1Coh->Scale(1.0/nEnergyBins);
		h1QFP->Scale(1.0/nEnergyBins);
		h1QFN->Scale(1.0/nEnergyBins);
		
		h_TheoryPrim.push_back(h1Prim);
		h_TheoryCoh.push_back(h1Coh);
		h_TheoryQFP.push_back(h1QFP);
		h_TheoryQFN.push_back(h1QFN);
		
		locBeamEnergy += m_beamEnergyBinSize;
	}
	
	delete fTheory;
	return 0;
}

double YieldFitter::YieldFitFunction(double *x, double *par)
{
	if(h_matrix==NULL) return 0.0;
	
	double reconAngle = x[0];
	/*
	if(reconAngle>0.55 && reconAngle<1.00) {
		TF1::RejectPoint();
		return 0.;
	}
	*/
	int reconAngleBin = h_matrix->GetYaxis()->FindBin(reconAngle);
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	
	double Gamma = par[0];
	double Acoh  = par[1];
	double AincP = par[2];
	double AincN = par[3];
	double Phi   = par[4];
	
	int locMinEnergyBin = h_matrix->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrix->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrix->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrix->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrix->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSection(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, AincP, AincN, Phi);
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
	double AincP = par[2];
	double AincN = par[3];
	double Phi   = par[4];
	
	int locMinEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFine->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locMatrix = h_matrixFine->GetBinContent(iThetaBin, reconAngleBin, iEnergyBin);
			double locCS     = GetCrossSection(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, AincP, AincN, Phi);
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
			double locCS     = GetCrossSectionInterference(iEnergyBin-locMinEnergyBin, iThetaBin, Gamma, Acoh, Phi);
			dNdTheta += (locMatrix * locCS * h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(locEnergy)));
		}
	}
	
	double yield = dNdTheta * m_luminosity * 0.985 * 0.981 * EtaAnalyzer::m_branchingRatio * angleBinSize * scaleFactor;
	return yield;
}

double YieldFitter::GetCrossSection(int energyBin, int thetaBin, 
	double Gamma, double Acoh, double AincP, double AincN, double Phi) 
{
	double locPrim = (Gamma/0.515) * h_TheoryPrim[energyBin]->GetBinContent(thetaBin);
	double locCoh  = Acoh * h_TheoryCoh[energyBin]->GetBinContent(thetaBin);
	double locIncP = AincP * h_TheoryQFP[energyBin]->GetBinContent(thetaBin);
	double locIncN = AincN * h_TheoryQFN[energyBin]->GetBinContent(thetaBin);
	double locInt  = 2.0*sqrt(locPrim*locCoh)*cos(Phi*TMath::DegToRad());
	return (locPrim+locCoh+locIncP+locIncN+locInt);
}

double YieldFitter::GetCrossSectionInterference(int energyBin, int thetaBin, 
	double Gamma, double Acoh, double Phi) 
{
	double locPrim = (Gamma/0.515) * h_TheoryPrim[energyBin]->GetBinContent(thetaBin);
	double locCoh  = Acoh * h_TheoryCoh[energyBin]->GetBinContent(thetaBin);
	double locInt  = 2.0*sqrt(locPrim*locCoh)*cos(Phi*TMath::DegToRad());
	return locInt;
}
