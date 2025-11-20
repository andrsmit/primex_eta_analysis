void StyleHist(TH1F*, int lineColor=kBlack, int markerStyle=4);
void StyleCanvas(TCanvas*);

void plot1(int vetoOption=6)
{
	TFile *fIn = new TFile(Form("yield_phase3_VetoOption%d_new.root",vetoOption), "READ");
	
	TH1F *hyi  = (TH1F*)fIn->Get("AngularYieldInclusive");
	TH1F *hys  = (TH1F*)fIn->Get("AngularYieldFit");
	TH1F *hyb1 = (TH1F*)fIn->Get("EtaPionYield");
	TH1F *hyb  = (TH1F*)hyb1->Clone("hyb");
	if(vetoOption<6) {
		TH1F *hyb2 = (TH1F*)fIn->Get("HadronicBkgdYield");
		hyb->Add(hyb2);
	}
	
	StyleHist(hyi, kBlack);
	StyleHist(hys, kBlue);
	StyleHist(hyb, kRed);
	
	TCanvas *c1 = new TCanvas("c1", "c1", 950, 700);
	StyleCanvas(c1);
	hyi->Draw("PE1X0");
	hys->Draw("PE1X0 same");
	hyb->Draw("PE1X0 same");
	
	
	TH1F *hsum = (TH1F*)hys->Clone("hsum");
	hsum->Add(hyb);
	hsum->Draw("same PE1X0");
	
	
	
	/*
	TF1 *f1 = new TF1("fb", "pol5", 0.0, 4.5);
	f1->SetParLimits(0, 0.0, 1.e3);
	f1->FixParameter(1, 0.0);
	hyb->Fit(f1, "R0QL");
	//f1->Draw("same");
	
	TLegend *leg = new TLegend(0.154, 0.701, 0.441, 0.895);
	leg->AddEntry(hys, "Single-#eta Yield", "PE1X0");
	leg->AddEntry(hyb, "#eta#pi Yield", "PE1X0");
	leg->AddEntry(hyi, "Sum", "PE1X0");
	leg->SetBorderSize(0);
	leg->Draw();
	
	// construct exclusive yield as inclusive - background fit:
	
	TH1F *hy = (TH1F*)hyi->Clone("hy");
	for(int ibin=1; ibin<=hy->GetXaxis()->GetNbins(); ibin++) {
		double locX  = hy->GetBinCenter(ibin);
		double locY1 = hyi->GetBinContent(ibin);
		double locE1 = hyi->GetBinError(ibin);
		double locY2 = f1->Eval(locX);
		double locE2 = hyb->GetBinError(ibin);
		hy->SetBinContent(ibin, locY1-locY2);
		hy->SetBinError(ibin, sqrt(pow(locE1,2.0)+pow(locE2,2.0)));
	}
	hy->Draw("PE1X0");
	//hys->Draw("PE1X0 same");
	
	TH1F *hr = (TH1F*)hy->Clone("hr");
	hr->Divide(hys);
	//hr->Draw("PE1X0");
	*/
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
