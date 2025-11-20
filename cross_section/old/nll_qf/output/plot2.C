void StyleHist(TH1F*, int lineColor=kBlack, int markerStyle=4);
void StyleCanvas(TCanvas*);

void GetCrossSection(TH1F*,TH1F*);

void plot2()
{
	int veto1 = 4;
	int veto2 = 6;
	
	//---------------------------------------------------//
	
	TFile *fIn1 = new TFile(Form("yield_phase3_VetoOption%d_new.root",veto1), "READ");
	
	TH1F *hYield1 = (TH1F*)fIn1->Get("AngularYieldFit")->Clone("yield1");
	TH1F *hb11 = (TH1F*)fIn1->Get("EtaPionYield")->Clone("hb11");
	hYield1->Add(hb11);
	
	if(veto1<6) {
		TH1F *hb12 = (TH1F*)fIn1->Get("HadronicBkgdYield")->Clone("hb12");
		hYield1->Add(hb12);
	}
	
	TH1F *hEtaPi1 = (TH1F*)fIn1->Get("EtaPionYield")->Clone("hEtaPi1");
	TF1  *fEtaPi1 = new TF1("fEtaPi1", "pol5", 0.0, 4.5);
	fEtaPi1->SetParLimits(0, 0.0, 1.e3);
	fEtaPi1->FixParameter(1, 0.0);
	hEtaPi1->Fit(fEtaPi1, "R0QL");
	
	TH1F *hEtaPiPi1 = (TH1F*)fIn1->Get("HadronicBkgdYield")->Clone("hEtaPiPi1");
	TF1  *fEtaPiPi1 = new TF1("fEtaPiPi1", "pol5", 0.0, 4.5);
	fEtaPiPi1->SetParLimits(0, 0.0, 1.e3);
	fEtaPiPi1->FixParameter(1, 0.0);
	if(veto1) hEtaPiPi1->Fit(fEtaPiPi1, "R0QL");
	
	for(int ibin=1; ibin<=hYield1->GetXaxis()->GetNbins(); ibin++) {
		double locYield = hYield1->GetBinContent(ibin);
		double locBG1   = fEtaPi1->Eval(hYield1->GetBinCenter(ibin));
		double locBG2   = veto1<6 ? fEtaPiPi1->Eval(hYield1->GetBinCenter(ibin)) : 0.0;
		hYield1->SetBinContent(ibin, locYield - locBG1 - locBG2);
	}
	
	TH1F *hAcc1 = (TH1F*)fIn1->Get("h_Acceptance")->Clone("hAcc1");
	
	hYield1->SetDirectory(0);
	hAcc1->SetDirectory(0);
	fIn1->Close();
	
	//---------------------------------------------------//
	
	TFile *fIn2 = new TFile(Form("yield_phase3_VetoOption%d_new.root",veto2), "READ");
	
	TH1F *hYield2 = (TH1F*)fIn2->Get("AngularYieldFit")->Clone("yield2");
	TH1F *hb21 = (TH1F*)fIn2->Get("EtaPionYield")->Clone("hb21");
	hYield2->Add(hb21);
	
	if(veto2<6) {
		TH1F *hb22 = (TH1F*)fIn2->Get("HadronicBkgdYield")->Clone("hb22");
		hYield2->Add(hb22);
	}
	
	TH1F *hEtaPi2 = (TH1F*)fIn2->Get("EtaPionYield")->Clone("hEtaPi2");
	TF1  *fEtaPi2 = new TF1("fEtaPi2", "pol5", 0.0, 4.5);
	fEtaPi2->SetParLimits(0, 0.0, 1.e3);
	fEtaPi2->FixParameter(1, 0.0);
	hEtaPi2->Fit(fEtaPi2, "R0QL");
	
	TH1F *hEtaPiPi2 = (TH1F*)fIn2->Get("HadronicBkgdYield")->Clone("hEtaPiPi2");
	TF1  *fEtaPiPi2 = new TF1("fEtaPiPi2", "pol5", 0.0, 4.5);
	fEtaPiPi2->SetParLimits(0, 0.0, 1.e3);
	fEtaPiPi2->FixParameter(1, 0.0);
	if(veto2) hEtaPiPi2->Fit(fEtaPiPi2, "R0QL");
	
	for(int ibin=1; ibin<=hYield2->GetXaxis()->GetNbins(); ibin++) {
		double locYield = hYield2->GetBinContent(ibin);
		double locBG1   = fEtaPi2->Eval(hYield2->GetBinCenter(ibin));
		double locBG2   = veto2<6 ? fEtaPiPi2->Eval(hYield2->GetBinCenter(ibin)) : 0.0;
		hYield2->SetBinContent(ibin, locYield - locBG1 - locBG2);
	}
	
	TH1F *hAcc2 = (TH1F*)fIn2->Get("h_Acceptance")->Clone("hAcc2");
	
	hYield2->SetDirectory(0);
	hAcc2->SetDirectory(0);
	fIn2->Close();
	
	//---------------------------------------------------//
	
	StyleHist(hYield1, kRed);
	StyleHist(hYield2, kBlue);
	
	GetCrossSection(hYield1, hAcc1);
	GetCrossSection(hYield2, hAcc2);
	
	TCanvas *c1 = new TCanvas("c1", "c1", 950, 700);
	StyleCanvas(c1);
	hYield2->Draw("PE1X0");
	//hYield2->Draw("PE1X0 same");
	
	TLegend *leg1 = new TLegend(0.155, 0.713, 0.455, 0.904);
	leg1->AddEntry(hYield1, "BCAL Veto", "PE1X0");
	leg1->AddEntry(hYield2, "BCAL+SC Veto", "PE1X0");
	//leg1->Draw();
	
	return;
}

void GetCrossSection(TH1F *hY, TH1F *hA) {
	
	double BR = 0.3936;
	double lumi = 17.265406 / (1.e-6);
	double binSize = hY->GetBinWidth(1) * TMath::DegToRad();
	
	for(int ibin=1; ibin<=hY->GetXaxis()->GetNbins(); ibin++) {
		double locY = hY->GetBinContent(ibin);
		double locE = hY->GetBinError(ibin);
		double locA = hA->GetBinContent(ibin);
		double locC = locY / (locA * lumi * BR * binSize);
		double locCErr = locE / (locA * lumi * BR * binSize);
		hY->SetBinContent(ibin, locC);
		hY->SetBinError(ibin, locCErr);
	}
	hY->GetYaxis()->SetTitle("d#sigma/d#theta [#mub/rad]");
	return;
}

void StyleHist(TH1F *h1, int lineColor, int markerStyle) {
	
	for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
		double binC = h1->GetBinContent(ibin);
		double binE = h1->GetBinError(ibin);
		
		// check for nans:
		if(binE!=binE) {
			h1->SetBinError(ibin, 1.5*sqrt(binC));
			continue;
		}
		
		// make sure the error bars are not smaller than sqrt(yield):
		if(sqrt(binC) > binE) {
			h1->SetBinError(ibin, sqrt(binC));
			continue;
		}
		
		// make sure the error bars are not anomalously large:
		if(binE > 2.5*sqrt(binC)) h1->SetBinError(ibin, 2.5*sqrt(binC));
	}
	
	h1->SetDirectory(0);
	
	h1->SetLineColor(lineColor);
	h1->SetMarkerColor(lineColor);
	h1->SetMarkerStyle(markerStyle);
	
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetXaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetMaxDigits(2);
}

void StyleCanvas(TCanvas *c1)
{
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	c1->SetTickx(); c1->SetTicky();
	c1->cd();
}
