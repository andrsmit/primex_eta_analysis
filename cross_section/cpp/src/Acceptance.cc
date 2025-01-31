#include "EtaAnalyzer.h"
#include "CrossSection.h"

int EtaAnalyzer::GetAcceptance(TString matrixName)
{
	if(m_binningSet==false) InitializeBinning();
	
	h_Acceptance = new TH1F("h_Acceptance", ";#theta(rec) [#circ]", m_angularBin.size(), 0.0, m_thetaBinSize*(double)(m_angularBin.size()));
	h_Acceptance->SetDirectory(0);
	h_Acceptance->SetLineColor(kRed);
	h_Acceptance->SetLineWidth(2);
	
	h_Acceptance->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_Acceptance->GetXaxis()->SetTitleSize(0.05);
	h_Acceptance->GetXaxis()->CenterTitle(true);
	h_Acceptance->GetYaxis()->SetTitle("Acceptance");
	h_Acceptance->GetYaxis()->SetTitleSize(0.05);
	h_Acceptance->GetYaxis()->CenterTitle(true);
	
	//------------------------------------------------------------//
	// ROOT file where angular and acceptance matrix is stored:
	
	TString matrixFileName = "phase1_matrix.root";
	if(matrixName.Contains("FCAL")) {
		matrixFileName = "phase1_FCAL.root";
	} else if(matrixName.Contains("TOF")) {
		matrixFileName = "phase1_TOF.root";
	}
	
	TString matrixFileNameFull = Form(
		"/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles/phase1/%s", 
		matrixFileName.Data());
	
	// return if filename is not accessible:
	if(gSystem->AccessPathName(matrixFileNameFull.Data())) return 1;
	
	TFile *fMatrix = new TFile(matrixFileNameFull.Data(), "READ");
	
	int nRebinsThrown = 2;
	
	TH3F *hRecMatrix = (TH3F*)fMatrix->Get(matrixName.Data());
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
	double    locRecAngleBinSize = hRecMatrix->GetYaxis()->GetBinWidth(1);
	double  locBeamEnergyBinSize = hRecMatrix->GetZaxis()->GetBinWidth(1);
	
	//---------------------------------------------------------------------//
	
	int locMinEnergyBin = hRecMatrix->GetZaxis()->FindBin(minBeamEnergy + 0.5*locBeamEnergyBinSize);
	int locMaxEnergyBin = hRecMatrix->GetZaxis()->FindBin(maxBeamEnergy - 0.5*locBeamEnergyBinSize);
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
		
		double minRecAngle = hRecMatrix->GetYaxis()->GetBinCenter(iRecAngleBin) - 0.5*locRecAngleBinSize;
		double maxRecAngle = hRecMatrix->GetYaxis()->GetBinCenter(iRecAngleBin) + 0.5*locRecAngleBinSize;
		
		int minThrownAngleBin = hRecMatrix->GetXaxis()->FindBin(minRecAngle + 0.5*locThrownAngleBinSize);
		int maxThrownAngleBin = hRecMatrix->GetXaxis()->FindBin(maxRecAngle - 0.5*locThrownAngleBinSize);
		
		double locAccSum        = 0.0;
		double locAccESum       = 0.0;
		double locFluxWeightSum = 0.0;
		
		for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
			
			double beamEnergy = hRecMatrix->GetZaxis()->GetBinCenter(iEnergyBin);
			
			double locRec    = hRecMatrix->Integral(1,hRecMatrix->GetXaxis()->GetNbins(), iRecAngleBin,iRecAngleBin, iEnergyBin,iEnergyBin);
			double locThrown = hThrownAngle_vs_ThrownEnergy->Integral(
				hThrownAngle_vs_ThrownEnergy->GetXaxis()->FindBin(beamEnergy), hThrownAngle_vs_ThrownEnergy->GetXaxis()->FindBin(beamEnergy), 
				minThrownAngleBin, maxThrownAngleBin);
			if(locThrown<=0.0) continue;
			
			double locAcc  = locRec/locThrown;
			double locAccE = sqrt(locRec*(1.0 - locAcc))/locThrown;
			
			locAccSum  += locAcc*m_fluxWeights[iEnergyBin-locMinEnergyBin];
			locAccESum += pow(locAccE*m_fluxWeights[iEnergyBin-locMinEnergyBin],2.0);
			locFluxWeightSum += m_fluxWeights[iEnergyBin-locMinEnergyBin];
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
	
	h_Acceptance->GetYaxis()->SetRangeUser(0.0,1.0);
	
	if(cAcceptance==NULL) {
		cAcceptance = new TCanvas("cAcceptance", "Acceptance", 950, 700);
		styleCanvas(cAcceptance);
	}
	h_Acceptance->Draw("PE");
	
	return 0;
}
