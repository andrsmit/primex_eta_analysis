void AddYield(TH1F*, TH1F*);
void AddCrossSection(TH1F*, TH1F*, TH1F*);

void compareEnergies()
{
	int vetoOption = 6;
	int phase = 3;
	
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	
	TFile *f0 = new TFile(Form("yield_phase%d_VetoOption%d_8.0GeV_11.3GeV.root", phase, vetoOption), "READ");
	TH1F  *hc0 = (TH1F*)f0->Get("CrossSectionFit")->Clone("cs0");
	hc0->SetDirectory(0);
	TH1F  *hy0 = (TH1F*)f0->Get("AngularYieldFit")->Clone("yield0");
	hy0->SetDirectory(0);
	TH1F  *hb0 = (TH1F*)f0->Get("EtaPionYield")->Clone("bkgd0");
	hb0->SetDirectory(0);
	f0->Close();
	
	TFile *f1 = new TFile(Form("yield_phase%d_VetoOption%d_10.0GeV_11.3GeV.root", phase, vetoOption), "READ");
	TH1F  *hc1 = (TH1F*)f1->Get("CrossSectionFit")->Clone("cs1");
	hc1->SetDirectory(0);
	TH1F  *hy1 = (TH1F*)f1->Get("AngularYieldFit")->Clone("yield1");
	hy1->SetDirectory(0);
	TH1F  *hb1 = (TH1F*)f1->Get("EtaPionYield")->Clone("bkgd1");
	hb1->SetDirectory(0);
	f1->Close();
	
	TFile *f2 = new TFile(Form("yield_phase%d_VetoOption%d_9.0GeV_10.0GeV.root", phase, vetoOption), "READ");
	TH1F  *hc2 = (TH1F*)f2->Get("CrossSectionFit")->Clone("cs2");
	hc2->SetDirectory(0);
	TH1F  *hy2 = (TH1F*)f2->Get("AngularYieldFit")->Clone("yield2");
	hy2->SetDirectory(0);
	TH1F  *hb2 = (TH1F*)f2->Get("EtaPionYield")->Clone("bkgd2");
	hb2->SetDirectory(0);
	f2->Close();
	
	TFile *f3 = new TFile(Form("yield_phase%d_VetoOption%d_8.0GeV_9.0GeV.root", phase, vetoOption), "READ");
	TH1F  *hc3 = (TH1F*)f3->Get("CrossSectionFit")->Clone("cs3");
	hc3->SetDirectory(0);
	TH1F  *hy3 = (TH1F*)f3->Get("AngularYieldFit")->Clone("yield3");
	hy3->SetDirectory(0);
	TH1F  *hb3 = (TH1F*)f3->Get("EtaPionYield")->Clone("bkgd3");
	hb3->SetDirectory(0);
	f3->Close();
	
	
	if(0) {
		AddCrossSection(hc0, hy0, hb0);
		AddCrossSection(hc1, hy1, hb1);
		AddCrossSection(hc2, hy2, hb2);
		AddCrossSection(hc3, hy3, hb3);
	}
	
	hc0->SetLineWidth(1);
	hc0->SetLineColor(kBlack);
	hc0->SetMarkerColor(kBlack);
	hc0->SetMarkerStyle(8);
	hc0->SetMarkerSize(0.7);
	
	hc1->SetLineWidth(1);
	hc1->SetLineColor(kBlue);
	hc1->SetMarkerColor(kBlue);
	hc1->SetMarkerStyle(8);
	hc1->SetMarkerSize(0.7);
	
	hc2->SetLineWidth(1);
	hc2->SetLineColor(kRed);
	hc2->SetMarkerColor(kRed);
	hc2->SetMarkerStyle(26);
	hc2->SetMarkerSize(0.7);
	
	hc3->SetLineWidth(1);
	hc3->SetLineColor(kBlack);
	hc3->SetMarkerColor(kBlack);
	hc3->SetMarkerStyle(4);
	hc3->SetMarkerSize(0.7);
	
	hc0->GetYaxis()->SetRangeUser(0.,2.2);
	
	c->cd();
	//hc0->Draw("PE1X0");
	hc1->Draw("PE1X0");
	hc2->Draw("PE1X0 same");
	hc3->Draw("PE1X0 same");
	
	TLegend *leg1 = new TLegend(0.13, 0.73, 0.43, 0.93);
	//leg1->AddEntry(hc0, " 8.0 - 11.3 GeV", "PE1X0");
	leg1->AddEntry(hc3, " 8.0 -  9.0 GeV", "PE1X0");
	leg1->AddEntry(hc2, " 9.0 - 10.0 GeV", "PE1X0");
	leg1->AddEntry(hc1, "10.0 - 11.3 GeV", "PE1X0");
	leg1->Draw();
	
	return;
}

void AddYield(TH1F *hy, TH1F *hb) {
	for(int ibin=1; ibin<=hy->GetXaxis()->GetNbins(); ibin++) {
		double y = hy->GetBinContent(ibin);
		double b = hb->GetBinContent(ibin);
		hy->SetBinContent(ibin, y+b);
	}
}
void AddCrossSection(TH1F *hc, TH1F *hy, TH1F *hb) {
	for(int ibin=1; ibin<=hc->GetXaxis()->GetNbins(); ibin++) {
		double c = hc->GetBinContent(ibin);
		double y = hy->GetBinContent(ibin);
		double b = hb->GetBinContent(ibin);
		hc->SetBinContent(ibin, c*(y+b)/y);
	}
}
