void AddYield(TH1F*, TH1F*);
void AddCrossSection(TH1F*, TH1F*, TH1F*);

void compareShifts(int vetoOption=6)
{
	gStyle->SetOptStat(0);
	
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	
	TLegend *leg1 = new TLegend(0.146, 0.703, 0.514, 0.963);
	
	vector<TString> ShiftStrings = {"1.5MeV","2.5MeV","3.5MeV",""};
	vector<int>  colorOptions = {kBlack, kBlue, kRed, kMagenta, kGreen+1};
	vector<int> markerOptions = {     8,    25,   22,       23,        8};
	
	
	TH1F *hmax, *hmin;
	
	for(int i=0; i<ShiftStrings.size(); i++) {
		
		TFile *fIn = new TFile(Form("yield_phase3_VetoOption%d_%s.root", vetoOption, ShiftStrings[i].Data()), "READ");
		
		TH1F *hc = (TH1F*)fIn->Get("CrossSectionFit")->Clone(Form("CrossSection_VetoOption%d",vetoOption));
		TH1F *hy = (TH1F*)fIn->Get("AngularYieldFit")->Clone(Form("AngularYield_VetoOption%d",vetoOption));
		TH1F *hb = (TH1F*)fIn->Get("EtaPionYield")->Clone(Form("EtaPionYield_VetoOption%d",vetoOption));
		
		hc->SetDirectory(0);
		hy->SetDirectory(0);
		hb->SetDirectory(0);
		fIn->Close();
		
		hc->SetLineWidth(1);
		hc->SetLineColor(colorOptions[i]);
		hc->SetMarkerColor(colorOptions[i]);
		hc->SetMarkerStyle(markerOptions[i]);
		
		hy->SetLineWidth(1);
		hy->SetLineColor(colorOptions[i]);
		hy->SetMarkerColor(colorOptions[i]);
		hy->SetMarkerStyle(markerOptions[i]);
		
		hb->SetLineWidth(1);
		hb->SetLineColor(colorOptions[i]);
		hb->SetMarkerColor(colorOptions[i]);
		hb->SetMarkerStyle(markerOptions[i]);
		
		if(0) {
			AddYield(hy,hb);
			AddCrossSection(hc,hy,hb);
		}
		
		if(i==0) hmin = (TH1F*)hc->Clone("hmin");
		if(i==(ShiftStrings.size()-1)) hmax = (TH1F*)hc->Clone("hmax");
		
		hc->GetYaxis()->SetRangeUser(0.0,2.4);
		
		if(i==0) {
			hc->Draw("PE1X0");
			//hy->Draw("PE1X0");
			//hb->Draw("PE1X0");
		}
		else {
			hc->Draw("PE1X0 same");
			//hy->Draw("PE1X0 same");
			//hb->Draw("PE1X0 same");
		}
		TString locString = ShiftStrings[i];
		if(locString=="") {
			locString = "Variable shift";
		}
		else {
			locString = Form("%s shift", ShiftStrings[i].Data());
		}
		leg1->AddEntry(hc, Form("%s",locString.Data()), "PE1X0");
	}
	leg1->Draw();
	
	hmax->Add(hmin,-1.0);
	hmax->Draw("same hist");
	
	TLatex lat;
	lat.DrawLatexNDC(0.55, 0.825, "#scale[1.0]{Modified BCAL Veto}");
	
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
