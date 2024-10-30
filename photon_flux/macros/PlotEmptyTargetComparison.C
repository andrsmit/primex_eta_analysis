/*
One of the questions I asked is whether or not the 2-photon invariant mass spectrum for the 
empty target data has strong energy-dependence.

In this ROOT macro, I plot a comparison of the invariant mass
distribution with different beam energies for comparison.

Motivation:
Because of the variation in the endpoint energy over the full run period,
the ratio of photon flux collected with the target full to the photon flux 
collected with the target empty is energy-dependent. Therefore, to properly subtract
the empty target background, the background distribution should be scaled by 
an energy-dependent factor. This is accomplished by binning the invariant mass spectrum
on beam energy, then scaling the distributions in each energy bin by the ratio of 
integrated flux on the full target to the integrated flux on the empty target
within that beam energy range.
The scaled distributions in each beam energy range are then summed and combined into one 
overall distribution.
This is expected to be a small improvement over the normal approach of scaling the integrated
empty target distribution by the ratio of flux with full target to the flux on empty target, integrated
over the full energy range.
But just to check, a comparison is drawn between these two cases in this macro.

Conclusion:
It really makes almost no difference whether the scaling of the empty target
background is applied bin-by-bin or integrated over the full energy range.

This is done by using the TH3F object created by the "FillInvmassMatrix" function in EtaAna.cc.
The 3-D histogram plots invmass vs. production angle and beam energy:
	X-axis: production angle (600 bins from 0.0 - 6.0 degrees)
	Y-axis: 2-photon invariant mass with energy-constraint (1200 bins from 0.0 - 1.2 GeV/c^2)
	Z-axis: Beam Energy (50 bins from 7.0 - 12.0 GeV) 
*/

// 

void PlotEmptyTargetComparison(int phase=1, double minAngle=0.0, double maxAngle=6.0) {
	
	TString fieldString = phase==1 ? "nobfield" : "bfield";
	
	TString primexDir = getenv("PRIMEXDIR");
	TString emptyTargetFilename = primexDir + Form("/analyze_trees/rootFiles/phase%d/empty_target_%s.root", phase, fieldString.Data());
	
	if(gSystem->AccessPathName(emptyTargetFilename)) {
		std::cout << "Problem accessing empty target ROOT file." << std::endl;
		return;
	}
	
	//--------------------------------------------------//
	// Get invariant mass matrix:
	
	TFile *fEmpty = new TFile(emptyTargetFilename.Data(), "READ");
	
	TH3F *h3Invmass = (TH3F*)fEmpty->Get("invmassMatrix");
	if(h3Invmass==NULL) {
		std::cout << "Unable to get invmassMatrix from ROOT file." << std::endl;
		fEmpty->Close();
		return;
	}
	h3Invmass->SetDirectory(0);
	fEmpty->Close();
	
	//--------------------------------------------------//
	
	double angleBinSize = h3Invmass->GetXaxis()->GetBinCenter(2) - h3Invmass->GetXaxis()->GetBinCenter(1);
	int angleBinLo = h3Invmass->GetXaxis()->FindBin(minAngle+0.5*angleBinSize);
	int angleBinHi = h3Invmass->GetXaxis()->FindBin(maxAngle-0.5*angleBinSize);
	
	double beamEnergyBinSize = h3Invmass->GetZaxis()->GetBinCenter(2) - h3Invmass->GetZaxis()->GetBinCenter(1);
	int energyBin1A = h3Invmass->GetZaxis()->FindBin( 8.0 + 0.5*beamEnergyBinSize);
	int energyBin1B = h3Invmass->GetZaxis()->FindBin( 9.0 - 0.5*beamEnergyBinSize);
	int energyBin2A = h3Invmass->GetZaxis()->FindBin( 9.0 + 0.5*beamEnergyBinSize);
	int energyBin2B = h3Invmass->GetZaxis()->FindBin(10.0 - 0.5*beamEnergyBinSize);
	int energyBin3A = h3Invmass->GetZaxis()->FindBin(10.0 + 0.5*beamEnergyBinSize);
	int energyBin3B = h3Invmass->GetZaxis()->FindBin(11.0 - 0.5*beamEnergyBinSize);
	
	TH1F *hInvmass1 = (TH1F*)h3Invmass->ProjectionY("Invmass_low",  angleBinLo, angleBinHi, energyBin1A, energyBin1B);
	TH1F *hInvmass2 = (TH1F*)h3Invmass->ProjectionY("Invmass_mid",  angleBinLo, angleBinHi, energyBin2A, energyBin2B);
	TH1F *hInvmass3 = (TH1F*)h3Invmass->ProjectionY("Invmass_high", angleBinLo, angleBinHi, energyBin3A, energyBin3B);
	
	double locMax = hInvmass1->GetMaximum();
	if(hInvmass2->GetMaximum() > locMax) locMax = hInvmass2->GetMaximum();
	if(hInvmass3->GetMaximum() > locMax) locMax = hInvmass3->GetMaximum();
	
	hInvmass1->GetXaxis()->SetTitleSize(0.05);
	hInvmass1->GetXaxis()->SetTitleOffset(1.0);
	hInvmass1->GetXaxis()->CenterTitle(true);
	hInvmass1->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
	hInvmass1->Rebin(5);
	
	hInvmass2->GetXaxis()->SetTitleSize(0.05);
	hInvmass2->GetXaxis()->SetTitleOffset(1.0);
	hInvmass2->GetXaxis()->CenterTitle(true);
	hInvmass2->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
	hInvmass2->Rebin(5);
	
	hInvmass3->GetXaxis()->SetTitleSize(0.05);
	hInvmass3->GetXaxis()->SetTitleOffset(1.0);
	hInvmass3->GetXaxis()->CenterTitle(true);
	hInvmass3->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
	hInvmass3->Rebin(5);
	
	hInvmass1->SetMarkerColor(kBlue);
	hInvmass1->SetLineColor(kBlue);
	
	hInvmass2->SetMarkerColor(kRed);
	hInvmass2->SetLineColor(kRed);
	
	hInvmass3->SetMarkerColor(kBlack);
	hInvmass3->SetLineColor(kBlack);
	
	TCanvas *cInvmass = new TCanvas("cInvmass", "Invmass", 950, 700);
	cInvmass->SetTickx(); cInvmass->SetTicky();
	cInvmass->SetLeftMargin(0.13); cInvmass->SetRightMargin(0.07);
	cInvmass->SetBottomMargin(0.13); cInvmass->SetTopMargin(0.07);
	
	hInvmass1->Draw("hist");
	hInvmass2->Draw("PE same");
	hInvmass3->Draw("PE same");
	
	//--------------------------------------------------//
	// Loop over each beam energy bin and scale each 1-D projection of invmass by the ratio 
	// of photon flux on full to flux on empty in that energy bin:
	
	TString  fullTargetFluxFilename = primexDir + Form("/photon_flux/rootFiles/phase%d/full.root",  phase);
	TString emptyTargetFluxFilename = primexDir + Form("/photon_flux/rootFiles/phase%d/empty.root", phase);
	
	if(gSystem->AccessPathName(fullTargetFluxFilename.Data()) || gSystem->AccessPathName(emptyTargetFluxFilename.Data())) {
		std::cout << "Problem accessing flux ROOT files." << std::endl;
		return;
	}
	
	TFile *fFluxFull  = new TFile( fullTargetFluxFilename.Data(), "READ");
	TH1F  *hFluxFull  = (TH1F*)fFluxFull->Get("flux_vs_egamma")->Clone("hFluxFull");
	hFluxFull->SetDirectory(0);
	fFluxFull->Close();
	
	TFile *fFluxEmpty = new TFile(emptyTargetFluxFilename.Data(), "READ");
	TH1F  *hFluxEmpty = (TH1F*)fFluxEmpty->Get("flux_vs_egamma")->Clone("hFluxEmpty");
	hFluxEmpty->SetDirectory(0);
	fFluxEmpty->Close();
	
	// match bin size to invariant mass matrix:
	hFluxFull->Rebin(100);
	hFluxEmpty->Rebin(100);
	
	int minEnergyBin = h3Invmass->GetZaxis()->FindBin( 8.0 + 0.5*beamEnergyBinSize);
	int maxEnergyBin = h3Invmass->GetZaxis()->FindBin(10.9 - 0.5*beamEnergyBinSize);
	
	double integratedFluxFull  = hFluxFull->Integral(
		hFluxFull->GetXaxis()->FindBin(h3Invmass->GetZaxis()->GetBinCenter(minEnergyBin)),
		hFluxFull->GetXaxis()->FindBin(h3Invmass->GetZaxis()->GetBinCenter(maxEnergyBin)));
	double integratedFluxEmpty = hFluxEmpty->Integral(
		hFluxEmpty->GetXaxis()->FindBin(h3Invmass->GetZaxis()->GetBinCenter(minEnergyBin)),
		hFluxEmpty->GetXaxis()->FindBin(h3Invmass->GetZaxis()->GetBinCenter(maxEnergyBin)));
	double integratedFluxRatio = integratedFluxFull / integratedFluxEmpty;
	//printf("Integrated Flux Ratio = %f\n", integratedFluxRatio);
	
	TH1F *hInvmass = (TH1F*)h3Invmass->ProjectionY("hInvmass", angleBinLo, angleBinHi, minEnergyBin, maxEnergyBin);
	hInvmass->Scale(integratedFluxFull/integratedFluxEmpty);
	
	TH1F *hInvmassCorr = (TH1F*)hInvmass->Clone("hInvmassCorr");
	hInvmassCorr->Reset();
	
	for(int iEnergyBin = minEnergyBin; iEnergyBin <= maxEnergyBin; iEnergyBin++) {
		int locFluxBin = hFluxFull->GetXaxis()->FindBin(h3Invmass->GetZaxis()->GetBinCenter(iEnergyBin));
		TH1F *locProj = (TH1F*)h3Invmass->ProjectionY(Form("locProj_%d",iEnergyBin), angleBinLo, angleBinHi, iEnergyBin, iEnergyBin);
		double locFluxRatio = hFluxFull->GetBinContent(locFluxBin) / hFluxEmpty->GetBinContent(locFluxBin);
		locProj->Scale(locFluxRatio);
		hInvmassCorr->Add(locProj);
		//printf("Energy bin: %d (Egamma = %.2f)\n", iEnergyBin, h3Invmass->GetZaxis()->GetBinCenter(iEnergyBin));
		//printf("  Flux ratio = %f\n", locFluxRatio);
	}
	
	cInvmass->cd();
	
	hInvmass->GetXaxis()->SetTitleSize(0.05);
	hInvmass->GetXaxis()->SetTitleOffset(1.0);
	hInvmass->GetXaxis()->CenterTitle(true);
	hInvmass->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
	hInvmass->Rebin(1);
	
	hInvmassCorr->GetXaxis()->SetTitleSize(0.05);
	hInvmassCorr->GetXaxis()->SetTitleOffset(1.0);
	hInvmassCorr->GetXaxis()->CenterTitle(true);
	hInvmassCorr->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
	hInvmassCorr->Rebin(1);
		
	hInvmass->SetMarkerColor(kBlue);
	hInvmass->SetLineColor(kBlue);
	hInvmass->SetLineWidth(2);
	
	hInvmassCorr->SetMarkerColor(kRed);
	hInvmassCorr->SetLineColor(kRed);
	hInvmassCorr->SetLineWidth(2);
	
	hInvmass->Draw("hist");
	hInvmassCorr->Draw("same hist");
	
	return;
}
