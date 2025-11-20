void AddYield(TH1F *hy, TH1F *hb);

void getBGGENBkgd(int vetoOption=2) {
	
	gStyle->SetOptStat(0);
	
	TFile *fData = new TFile(Form("yield_phase3_VetoOption%d.root",vetoOption), "READ");
	TH1F *hYield = (TH1F*)fData->Get("AngularYield");
	TH1F *hEtaPi = (TH1F*)fData->Get("EtaPionYield");
	TH1F *hHadBG = (TH1F*)fData->Get("HadronicBkgdYield");
	
	TH1F *hYieldFit = (TH1F*)fData->Get("AngularYieldFit");
	/*
	AddYield(hYield, hEtaPi);
	AddYield(hYield, hHadBG);
	*/
	TFile *fIn = new TFile(Form(
		"/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana/analyze_trees/rootFiles/phase3/noBCALThreshold/Helium_VetoOption%d.root", 
		vetoOption), "READ");
	
	TH1F *hthrown = (TH1F*)fIn->Get("thrown_reactions_bggen");
	double nthrown = hthrown->Integral();
	double cs      = 1.18673e+08;
	double lumi    = 10.0*nthrown / cs;
	
	TH2F *hhSignal = (TH2F*)fIn->Get("mgg_const_bggen_signal");
	TH2F *hhEtaPi  = (TH2F*)fIn->Get("mgg_const_bggen_etapion");
	TH2F *hhEta2Pi = (TH2F*)fIn->Get("mgg_const_bggen_eta2pion");
	TH2F *hhEta3Pi = (TH2F*)fIn->Get("mgg_const_bggen_eta3pion");
	TH2F *hhOmega  = (TH2F*)fIn->Get("mgg_const_bggen_omega");
	TH2F *hhRho    = (TH2F*)fIn->Get("mgg_const_bggen_rho");
	TH2F *hhBkgd   = (TH2F*)fIn->Get("mgg_const_bggen_bkgd");
	hhBkgd->Add(hhRho);
	hhBkgd->Add(hhOmega);
	hhBkgd->Add(hhEta3Pi);
	hhBkgd->Add(hhEta2Pi);
	hhBkgd->Add(hhEtaPi);
	
	TH1F *hSignal = (TH1F*)hhSignal->ProjectionX("hSignal", 201,300);
	TH1F *hBkgd   = (TH1F*)hhBkgd->ProjectionX("hBkgd", 201,300);
	
	TH1F *hTotal = (TH1F*)hSignal->Clone("hTotal");
	hTotal->Add(hBkgd);
	
	hTotal->SetLineColor(kBlack);
	hSignal->SetLineColor(kBlue);
	hBkgd->SetLineColor(kRed);
	
	hTotal->Rebin(4);
	hSignal->Rebin(4);
	hBkgd->Rebin(4);
	
	double data_lumi = 2.0*17.6;
	
	hTotal->Scale(data_lumi/lumi);
	hSignal->Scale(data_lumi/lumi);
	hBkgd->Scale(data_lumi/lumi);
	
	TCanvas *c1 = new TCanvas("c1","c1",950,700);
	c1->cd();
	hYield->Draw("PE1X0");
	//hTotal->Draw("PE1X0 same");
	//hSignal->Draw("PE1X0 same");
	hBkgd->Draw("PE1X0 same");
	
	double totalIntegral = hTotal->Integral();
	double signalIntegral = hSignal->Integral();
	double bkgdIntegral = hBkgd->Integral();
	
	TH1F *hBkgdSub = (TH1F*)hYield->Clone("hBkgdSub");
	for(int ibin=1; ibin<=hYield->GetXaxis()->GetNbins(); ibin++) {
		double y = hYield->GetBinContent(ibin);
		double b = hBkgd->GetBinContent(ibin);
		hBkgdSub->SetBinContent(ibin, y-b);
	}
	hBkgdSub->SetMarkerColor(kRed);
	hBkgdSub->SetLineColor(kRed);
	hBkgdSub->Draw("PE1X0 same");
	
	hYieldFit->SetMarkerColor(kBlack);
	hYieldFit->SetLineColor(kBlack);
	hYieldFit->Draw("PE1X0 same");
	
	printf("Bkgd Level: %f\n", bkgdIntegral/totalIntegral);
	
	return;
}
void AddYield(TH1F *hy, TH1F *hb) {
	for(int ibin=1; ibin<=hy->GetXaxis()->GetNbins(); ibin++) {
		double y = hy->GetBinContent(ibin);
		double b = hb->GetBinContent(ibin);
		hy->SetBinContent(ibin, y+b);
	}
}
