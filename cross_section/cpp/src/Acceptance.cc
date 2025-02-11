#include "EtaAnalyzer.h"
#include "CrossSection.h"

int EtaAnalyzer::LoadAngularMatrix()
{
	int locPhase = m_phase;
	
	//------------------------------------------------------------//
	// ROOT file where angular and acceptance matrix is stored:
	
	if(m_analysisOption==1 && m_matrixHistName!="AngularMatrix") {
		m_matrixHistName = "AngularMatrix";
	}
	
	TString matrixFileName = Form("phase%d_matrix.root", locPhase);
	if(m_matrixHistName.Contains("FCAL"))
	{
		matrixFileName = Form("phase%d_FCAL.root", locPhase);
	}
	else if(m_matrixHistName.Contains("TOF"))
	{
		matrixFileName = Form("phase%d_TOF.root", locPhase);
	}
	
	TString matrixFileNameFull = Form(
		"/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles/phase%d/%s", 
		locPhase, matrixFileName.Data());
	
	// return if filename is not accessible:
	if(gSystem->AccessPathName(matrixFileNameFull.Data())) return 1;
	
	//------------------------------------------------------------//
	
	TFile *fMatrix = new TFile(matrixFileNameFull.Data(), "READ");
	
	h_matrix = (TH3F*)fMatrix->Get(m_matrixHistName.Data());
	h_matrix->SetDirectory(0);
	
	// Make a copy of 'h_matrix' before rebinning:
	
	h_matrixFine = (TH3F*)h_matrix->Clone("matrixFine");
	h_matrixFine->SetDirectory(0);
	
	//-----------------------------------------------------------//
	// Figure out how to re-bin matrix according to specified bin sizes:
	
	int nRebinsRecon, nRebinsThrown, nRebinsEnergy;
	double testBinSize;
	
	// Thrown Angle:
	nRebinsThrown = (int)(m_thrownAngleBinSize / h_matrix->GetXaxis()->GetBinWidth(1));
	testBinSize   = (double)(nRebinsThrown * h_matrix->GetXaxis()->GetBinWidth(1));
	if(testBinSize != m_thrownAngleBinSize) {
		printf("\n\nWarning: Requested bin size for angular matrix is off.\n");
		printf("  Requested thrown angle bin size: %.3f\n", m_thrownAngleBinSize);
		printf("  Actual thrown angle bin size: %.3f\n", testBinSize);
	}
	
	// Reconstructed Angle:
	nRebinsRecon = (int)(m_reconAngleBinSize / h_matrix->GetYaxis()->GetBinWidth(1));
	testBinSize  = (double)(nRebinsRecon * h_matrix->GetYaxis()->GetBinWidth(1));
	if(testBinSize != m_reconAngleBinSize) {
		printf("\n\nWarning: Requested bin size for angular matrix is off.\n");
		printf("  Requested reconstructed angle bin size: %.3f\n", m_reconAngleBinSize);
		printf("  Actual reconstructed angle bin size: %.3f\n", testBinSize);
	}
	
	// Beam Energy:
	nRebinsEnergy = (int)(m_beamEnergyBinSize / h_matrix->GetZaxis()->GetBinWidth(1));
	testBinSize   = (double)(nRebinsEnergy * h_matrix->GetZaxis()->GetBinWidth(1));
	if(testBinSize != m_beamEnergyBinSize) {
		printf("\n\nWarning: Requested bin size for angular matrix is off.\n");
		printf("  Requested beam energy bin size: %.3f\n", m_beamEnergyBinSize);
		printf("  Actual beam energy bin size: %.3f\n", testBinSize);
	}
	
	printf("RebinsThrown = %d\n", nRebinsThrown);
	printf("RebinsRecon  = %d\n", nRebinsRecon);
	printf("RebinsEnergy = %d\n", nRebinsEnergy);
	
	
	h_matrix->RebinX(nRebinsThrown);
	h_matrix->RebinY(nRebinsRecon);
	h_matrix->RebinZ(nRebinsEnergy);
	h_matrix->SetDirectory(0);
	
	/*
	m_beamEnergyBinSize  = h_matrix->GetZaxis()->GetBinWidth(1);
	m_thrownAngleBinSize = h_matrix->GetXaxis()->GetBinWidth(1);
	m_reconAngleBinSize  = h_matrix->GetYaxis()->GetBinWidth(1);
	*/
	
	//-----------------------------------------------------------//
	
	TH2F *hThrown = (TH2F*)fMatrix->Get("thrown");
	hThrown->SetDirectory(0);
	
	for(int iEnergyBin=1; iEnergyBin<=h_matrix->GetZaxis()->GetNbins(); iEnergyBin++) {
		double locEnergy = h_matrix->GetZaxis()->GetBinCenter(iEnergyBin);
		
		for(int iThetaBin=1; iThetaBin<=h_matrix->GetXaxis()->GetNbins(); iThetaBin++) {
			double locTheta = h_matrix->GetXaxis()->GetBinCenter(iThetaBin);
			
			double locThrown = hThrown->GetBinContent(hThrown->GetXaxis()->FindBin(locEnergy), 
				hThrown->GetYaxis()->FindBin(locTheta));
			
			for(int iRecBin=1; iRecBin<=h_matrix->GetYaxis()->GetNbins(); iRecBin++) {
				double locAcc = 0.0, locAccErr = 0.0;
				if(locThrown > 0.0) {
					locAcc    = h_matrix->GetBinContent(iThetaBin, iRecBin, iEnergyBin) / locThrown;
					locAccErr = sqrt(locThrown * locAcc * (1.0-locAcc)) / locThrown;
				}
				h_matrix->SetBinContent(iThetaBin, iRecBin, iEnergyBin, locAcc);
				h_matrix->SetBinError(iThetaBin, iRecBin, iEnergyBin, locAccErr);
			}
		}
	}
	
	for(int iEnergyBin=1; iEnergyBin<=h_matrixFine->GetZaxis()->GetNbins(); iEnergyBin++) {
		double locEnergy = h_matrixFine->GetZaxis()->GetBinCenter(iEnergyBin);
		
		for(int iThetaBin=1; iThetaBin<=h_matrixFine->GetXaxis()->GetNbins(); iThetaBin++) {
			double locTheta = h_matrixFine->GetXaxis()->GetBinCenter(iThetaBin);
			
			double locThrown = hThrown->GetBinContent(hThrown->GetXaxis()->FindBin(locEnergy), 
				hThrown->GetYaxis()->FindBin(locTheta));
			
			for(int iRecBin=1; iRecBin<=h_matrixFine->GetYaxis()->GetNbins(); iRecBin++) {
				double locAcc = 0.0, locAccErr = 0.0;
				if(locThrown > 0.0) {
					locAcc    = h_matrixFine->GetBinContent(iThetaBin, iRecBin, iEnergyBin) / locThrown;
					locAccErr = sqrt(locThrown * locAcc * (1.0-locAcc)) / locThrown;
				}
				h_matrixFine->SetBinContent(iThetaBin, iRecBin, iEnergyBin, locAcc);
				h_matrixFine->SetBinError(iThetaBin, iRecBin, iEnergyBin, locAccErr);
			}
		}
	}
	
	fMatrix->Close();
	
	delete fMatrix;
	delete hThrown;
	
	return 0;
}

int EtaAnalyzer::CalcAcceptance()
{
	if(m_binningSet==false) InitializeBinning();
	
	h_Acceptance = new TH1F("h_Acceptance", ";#theta(rec) [#circ]", 
		m_angularBin.size(), 0.0, m_reconAngleBinSize*(double)(m_angularBin.size()));
	h_Acceptance->SetDirectory(0);
	h_Acceptance->SetLineColor(kRed);
	h_Acceptance->SetLineWidth(2);
	h_Acceptance->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_Acceptance->GetXaxis()->SetTitleSize(0.05);
	h_Acceptance->GetXaxis()->CenterTitle(true);
	h_Acceptance->GetYaxis()->SetTitle("Acceptance");
	h_Acceptance->GetYaxis()->SetTitleSize(0.05);
	h_Acceptance->GetYaxis()->CenterTitle(true);
	
	//if(LoadAngularMatrix()) return 1;
	
	//------------------------------------------------------------//
	// ROOT file where angular and acceptance matrix is stored:
	
	if(m_analysisOption==1 && m_matrixHistName!="AngularMatrix") {
		m_matrixHistName = "AngularMatrix";
	}
	
	TString matrixFileName = Form("phase%d_matrix.root", m_phase);
	if(m_matrixHistName.Contains("FCAL"))
	{
		matrixFileName = Form("phase%d_FCAL.root", m_phase);
	}
	else if(m_matrixHistName.Contains("TOF"))
	{
		matrixFileName = Form("phase%d_TOF.root", m_phase);
	}
	
	TString matrixFileNameFull = Form(
		"/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles/phase%d/%s", 
		m_phase, matrixFileName.Data());
	
	// return if filename is not accessible:
	if(gSystem->AccessPathName(matrixFileNameFull.Data())) return 1;
	
	//------------------------------------------------------------//
	
	TFile *fMatrix = new TFile(matrixFileNameFull.Data(), "READ");
	
	int nRebinsThrown = 1;
	
	TH3F *hRecMatrix = (TH3F*)fMatrix->Get(m_matrixHistName.Data());
	hRecMatrix->RebinX(nRebinsThrown);
	hRecMatrix->RebinY(m_rebinsTheta);
	hRecMatrix->SetDirectory(0);
	
	TH2F *hThrownAngle_vs_ThrownEnergy = (TH2F*)fMatrix->Get("thrown");
	hThrownAngle_vs_ThrownEnergy->RebinY(nRebinsThrown);
	hThrownAngle_vs_ThrownEnergy->SetDirectory(0);
	/*
	hThrownAngle_vs_ThrownEnergy is a 2-dimensional histogram where the 
	beam energy axis ranges from 5.0 GeV to 12.0 GeV in 0.1 GeV step sizes, and the 
	thrown angle axis ranges from 0.0 deg to 5.0 deg in 0.01 deg step sizes.
	*/
	
	double locThrownAngleBinSize = hRecMatrix->GetXaxis()->GetBinWidth(1);
	double  locReconAngleBinSize = hRecMatrix->GetYaxis()->GetBinWidth(1);
	double  locBeamEnergyBinSize = hRecMatrix->GetZaxis()->GetBinWidth(1);
	
	//---------------------------------------------------------------------//
	
	int locMinEnergyBin = hRecMatrix->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBeamEnergyBinSize);
	int locMaxEnergyBin = hRecMatrix->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBeamEnergyBinSize);
	double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin+1);
	
	//---------------------------------------------------------------------//
	// To get the acceptance as a function of the reconstructed angle, we sum the matrix over reconstructed angle bins for a particular beam energy
	// and divide by the sum of thrown etas in that thrown angular range and beam energy bin.
	// Then we take the flux-weighted average over all energy bins.
	// In doing so, we need to make sure that the bin edges of the reconstructed angle axis of hMatrix align with edges of the thrown
	// angle axis bins. This is because we need to integrate the number of thrown etas in the same angular range as defined by the 
	// reconstructed angle bin size. 
	// Hence, we require that m_rebins_theta must be divisible by n_RebinsThrown.
	
	for(int iRecAngleBin=1; iRecAngleBin<=h_Acceptance->GetXaxis()->GetNbins(); iRecAngleBin++) {
		
		double minRecAngle = hRecMatrix->GetYaxis()->GetBinCenter(iRecAngleBin) - 0.5*locReconAngleBinSize;
		double maxRecAngle = hRecMatrix->GetYaxis()->GetBinCenter(iRecAngleBin) + 0.5*locReconAngleBinSize;
		
		int minThrownAngleBin = hRecMatrix->GetXaxis()->FindBin(minRecAngle + 0.5*locThrownAngleBinSize);
		int maxThrownAngleBin = hRecMatrix->GetXaxis()->FindBin(maxRecAngle - 0.5*locThrownAngleBinSize);
		
		// check that the bins line up correctly:
		
		double locMinThrownAngle = hRecMatrix->GetXaxis()->GetBinCenter(minThrownAngleBin) - 0.5*locThrownAngleBinSize;
		double locMaxThrownAngle = hRecMatrix->GetXaxis()->GetBinCenter(maxThrownAngleBin) + 0.5*locThrownAngleBinSize;
		if((fabs(locMinThrownAngle-minRecAngle)>1.e-6) || (fabs(locMaxThrownAngle-maxRecAngle)>1.e-6)) {
			printf("\nWarning: Bin edges of thrown and reconstructed axes do not coincide.\n");
			printf("  rec bin edges: %f - %f\n", minRecAngle, maxRecAngle);
			printf("  thrown bin edges: %f - %f\n", locMinThrownAngle, locMaxThrownAngle);
		}
		
		double locAccSum        = 0.0;
		double locAccESum       = 0.0;
		double locFluxWeightSum = 0.0;
		
		for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
			
			double beamEnergy    = hRecMatrix->GetZaxis()->GetBinCenter(iEnergyBin);
			double locFluxWeight = h_fluxWeights->GetBinContent(h_fluxWeights->FindBin(beamEnergy));
			
			double locRec    = hRecMatrix->Integral(1,hRecMatrix->GetXaxis()->GetNbins(), iRecAngleBin,iRecAngleBin, iEnergyBin,iEnergyBin);
			double locThrown = hThrownAngle_vs_ThrownEnergy->Integral(
				hThrownAngle_vs_ThrownEnergy->GetXaxis()->FindBin(beamEnergy), hThrownAngle_vs_ThrownEnergy->GetXaxis()->FindBin(beamEnergy), 
				minThrownAngleBin, maxThrownAngleBin);
			if(locThrown<=0.0) continue;
			
			double locAcc  = locRec/locThrown;
			double locAccE = sqrt(locRec*(1.0 - locAcc))/locThrown;
			
			locAccSum  += locAcc*locFluxWeight;
			locAccESum += pow(locAccE*locFluxWeight,2.0);
			locFluxWeightSum += locFluxWeight;
		}
		
		locAccSum = locFluxWeightSum > 0.0 ? (locAccSum/locFluxWeightSum) : 0.0;
		
		//printf("Angle: %f-%f, Acceptance: %f\n", minRecAngle, maxRecAngle, locAccSum);
		
		if(locAccSum < 0.05) {
			h_Acceptance->SetBinContent(iRecAngleBin, 0.0);
			h_Acceptance->SetBinError(iRecAngleBin, 0.0);
		} else {
			h_Acceptance->SetBinContent(iRecAngleBin, locAccSum);
			h_Acceptance->SetBinError(iRecAngleBin, sqrt(locAccESum)/locFluxWeightSum);
		}
	}
	
	fMatrix->Close();
	
	delete fMatrix;
	delete hRecMatrix;
	delete hThrownAngle_vs_ThrownEnergy;
	
	h_Acceptance->GetYaxis()->SetRangeUser(0.0,1.0);
	
	if(cAcceptance==NULL) {
		cAcceptance = new TCanvas("cAcceptance", "Acceptance", 950, 700);
		styleCanvas(cAcceptance);
	}
	h_Acceptance->Draw("PE");
	
	return 0;
}
