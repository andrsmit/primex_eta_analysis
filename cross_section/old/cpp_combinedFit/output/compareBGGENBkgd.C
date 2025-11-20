void StyleHist(TH1F*, int lineColor=kBlack, int markerStyle=8);
void StyleCanvas(TCanvas*);

void compareBGGENBkgd(int vetoOption=6) {
	
	gStyle->SetOptStat(0);
	
	//------------------------------------------------------------------------------//
	
	TFile *fData = new TFile(Form("yield_phase3_VetoOption%d_8.0GeV_11.3GeV.root",vetoOption), "READ");
	TH1F *hYield = (TH1F*)fData->Get("AngularYield");
	TH1F *hEtaPi = (TH1F*)fData->Get("EtaPionYield");
	
	TH1F *hFracData = (TH1F*)hEtaPi->Clone("FracData");
	for(int ibin=1; ibin<hFracData->GetXaxis()->GetNbins(); ibin++) {
		double locFrac    = 0.0;
		double locFracErr = 0.0;
		
		double num    = hEtaPi->GetBinContent(ibin);
		double numErr = hEtaPi->GetBinError(ibin);
		double den    = hYield->GetBinContent(ibin);
		double denErr = hYield->GetBinError(ibin);
		if((den>0.0) && (denErr>0.0)) {
			locFrac    = num/den;
			locFracErr = sqrt(pow(numErr/den,2.0) + pow(denErr*num/(den*den),2.0));
		}
		hFracData->SetBinContent(ibin,locFrac);
		hFracData->SetBinError(ibin,locFracErr);
	}
	
	hFracData->SetDirectory(0);
	fData->Close();
	
	//------------------------------------------------------------------------------//
	
	
	TFile *fBGGEN = new TFile(Form(	"/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana/analyze_trees/rootFiles/phase3/Helium_VetoOption%d.root", 
		vetoOption), "READ");
	
	TH2F *hhSignal = (TH2F*)fBGGEN->Get("mgg_const_bggen_signal_cut");
	TH2F *hhEtaPi  = (TH2F*)fBGGEN->Get("mgg_const_bggen_etapion_cut");
	
	int minMggCutBin = hhSignal->GetYaxis()->FindBin(0.5);
	int maxMggCutBin = hhSignal->GetYaxis()->FindBin(0.6);
	
	TH1F *hSignal = (TH1F*)hhSignal->ProjectionX("hSignal", minMggCutBin, maxMggCutBin);
	TH1F *hBkgd   = (TH1F*)hhEtaPi->ProjectionX("hBkgd", minMggCutBin, maxMggCutBin);
	
	hSignal->Rebin(5);
	hBkgd->Rebin(5);
	
	TH1F *hFracBGGEN = (TH1F*)hBkgd->Clone("FracBGGEN");
	for(int ibin=1; ibin<hFracBGGEN->GetXaxis()->GetNbins(); ibin++) {
		double locFrac    = 0.0;
		double locFracErr = 0.0;
		
		double num    = hBkgd->GetBinContent(ibin);
		double numErr = hBkgd->GetBinError(ibin);
		double den    = hSignal->GetBinContent(ibin);
		double denErr = hSignal->GetBinError(ibin);
		if((den>0.0) && (denErr>0.0)) {
			locFrac    = num/den;
			locFracErr = sqrt(pow(numErr/den,2.0) + pow(denErr*num/(den*den),2.0));
		}
		hFracBGGEN->SetBinContent(ibin,locFrac);
		hFracBGGEN->SetBinError(ibin,locFracErr);
	}
	
	hFracBGGEN->SetDirectory(0);
	fBGGEN->Close();
	
	//------------------------------------------------------------------------------//
	
	StyleHist(hFracData, kBlack);
	StyleHist(hFracBGGEN, kRed);
	
	hFracData->GetYaxis()->SetRangeUser(0.0,1.0);
	
	TCanvas *c1 = new TCanvas("c1","c1",950,700);
	StyleCanvas(c1);
	hFracData->Draw("PE1X0");
	hFracBGGEN->Draw("PE1X0 same");
	
	TLegend *leg = new TLegend(0.174, 0.735, 0.487, 0.899);
	leg->AddEntry(hFracData, "Data (fit results)", "PE1X0");
	leg->AddEntry(hFracBGGEN, "BGGEN MC", "PE");
	leg->Draw();
	
	return;
}

void StyleHist(TH1F *h1, int lineColor, int markerStyle)
{
	h1->SetLineWidth(2);
	h1->SetLineColor(lineColor);
	h1->SetMarkerColor(lineColor);
	h1->SetMarkerStyle(markerStyle);
	
	h1->SetTitle("");
	
	h1->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetXaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetTitle("Fraction of #eta+#pi background");
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetMaxDigits(2);
}

void StyleCanvas(TCanvas *c1)
{
	c1->SetLeftMargin(0.15); c1->SetRightMargin(0.05);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	c1->SetTickx(); c1->SetTicky();
	c1->cd();
}
