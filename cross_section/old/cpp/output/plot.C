TString GetVetoString(int option);

void AddYield(TH1F *hy, TH1F *hb0, TH1F *hb1=nullptr, TH1F *hb2=nullptr, TH1F *hb3=nullptr);
void AddCrossSection(TH1F *hc, TH1F *hy, TH1F *hb0, TH1F *hb1=nullptr, TH1F *hb2=nullptr, TH1F *hb3=nullptr);

void plot()
{
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	c->SetGrid();
	
	TLegend *leg1 = new TLegend(0.146, 0.703, 0.514, 0.963);
	/*
	vector<int>   vetoOptions = {     2,     1,    4,        6,        7};
	vector<int>  colorOptions = {kBlack, kBlue, kRed, kMagenta, kGreen+1};
	vector<int> markerOptions = {     4,    25,   22,       23,        8};
	*/
	vector<int>   vetoOptions = {     4,        6};
	vector<int>  colorOptions = {kBlack, kBlue, kRed, kMagenta, kGreen+1};
	vector<int> markerOptions = {     4,    25,   22,       23,        8};
	
	for(int i=0; i<vetoOptions.size(); i++) {
		int locVetoOption = vetoOptions[i];
		
		TFile *fIn = new TFile(
			Form("yield_phase3_VetoOption%d.root", locVetoOption), "READ");
		
		TH1F *hc  = (TH1F*)fIn->Get("CrossSectionFit")->Clone(Form("CrossSection_VetoOption%d",locVetoOption));
		TH1F *hy  = (TH1F*)fIn->Get("AngularYieldFit")->Clone(Form("AngularYield_VetoOption%d",locVetoOption));
		
		TH1F *hb0 = (TH1F*)fIn->Get("BkgdYield")->Clone(Form("BkgdYield_VetoOption%d",locVetoOption));
		TH1F *hb1 = (TH1F*)fIn->Get("OmegaYield")->Clone(Form("OmegaYield_VetoOption%d",locVetoOption));
		TH1F *hb2 = (TH1F*)fIn->Get("HadronicBkgdYield")->Clone(Form("HadronicBkgdYield_VetoOption%d",locVetoOption));
		TH1F *hb3 = (TH1F*)fIn->Get("EtaPionYield")->Clone(Form("EtaPionYield_VetoOption%d",locVetoOption));
		
		hc->SetDirectory(0);
		hy->SetDirectory(0);
		hb0->SetDirectory(0);
		hb1->SetDirectory(0);
		hb2->SetDirectory(0);
		hb3->SetDirectory(0);
		fIn->Close();
		
		hc->SetLineWidth(1);
		hc->SetLineColor(colorOptions[i]);
		hc->SetMarkerColor(colorOptions[i]);
		hc->SetMarkerStyle(markerOptions[i]);
		
		hy->SetLineWidth(1);
		hy->SetLineColor(colorOptions[i]);
		hy->SetMarkerColor(colorOptions[i]);
		hy->SetMarkerStyle(markerOptions[i]);
		
		hb0->SetLineWidth(1);
		hb0->SetLineColor(colorOptions[i]);
		hb0->SetMarkerColor(colorOptions[i]);
		hb0->SetMarkerStyle(markerOptions[i]);
		
		hb1->SetLineWidth(1);
		hb1->SetLineColor(colorOptions[i]);
		hb1->SetMarkerColor(colorOptions[i]);
		hb1->SetMarkerStyle(markerOptions[i]);
		
		hb2->SetLineWidth(1);
		hb2->SetLineColor(colorOptions[i]);
		hb2->SetMarkerColor(colorOptions[i]);
		hb2->SetMarkerStyle(markerOptions[i]);
		
		hb3->SetLineWidth(1);
		hb3->SetLineColor(colorOptions[i]);
		hb3->SetMarkerColor(colorOptions[i]);
		hb3->SetMarkerStyle(markerOptions[i]);
		
		if(0) {
			//AddYield(hy,hb1,hb2);
			AddCrossSection(hc,hy,hb0,hb2,hb3);
		}
		
		hc->GetYaxis()->SetRangeUser(0.0,2.4);
		
		if(i==0) {
			hc->Draw("PE1X0");
			//hy->Draw("PE1X0");
			//hb1->Draw("PE1X0");
			//hb2->Draw("PE1X0");
		}
		else {
			hc->Draw("PE1X0 same");
			//hy->Draw("PE1X0 same");
			//hb1->Draw("PE1X0 same");
			//hb2->Draw("PE1X0 same");
		}
		leg1->AddEntry(hc, Form("%s",GetVetoString(locVetoOption).Data()), "PE1X0");
	}
	leg1->Draw();
	
	return;
}

void AddYield(TH1F *hy, TH1F *hb0, TH1F *hb1, TH1F *hb2, TH1F *hb3) {
	for(int ibin=1; ibin<=hy->GetXaxis()->GetNbins(); ibin++) {
		double y  =  hy->GetBinContent(ibin);
		double b0 = hb0->GetBinContent(ibin);
		double b1 = hb1->GetBinContent(ibin);
		double b2 = hb2->GetBinContent(ibin);
		double b3 = hb3->GetBinContent(ibin);
		hy->SetBinContent(ibin, y+b0+b1+b2+b3);
	}
}
void AddCrossSection(TH1F *hc, TH1F *hy, TH1F *hb0, TH1F *hb1, TH1F *hb2, TH1F *hb3) {
	for(int ibin=1; ibin<=hc->GetXaxis()->GetNbins(); ibin++) {
		double c  =  hc->GetBinContent(ibin);
		double y  =  hy->GetBinContent(ibin);
		double b0 = 0., b1  = 0., b2 = 0., b3 = 0.;
		if(hb0) b0 = hb0->GetBinContent(ibin);
		if(hb1) b1 = hb1->GetBinContent(ibin);
		if(hb2) b2 = hb2->GetBinContent(ibin);
		if(hb3) b3 = hb3->GetBinContent(ibin);
		hc->SetBinContent(ibin, c*(y+b0+b1+b2+b3)/y);
	}
}

TString GetVetoString(int option)
{
	switch(option) {
		case 0:
			return "No Veto";
		case 1:
			return "Strict BCAL Veto";
		case 2:
			return "Loose BCAL Veto";
		case 3:
			return "1 ns BCAL w/ coplanarity";
		case 4:
			return "Loose BCAL Veto";
		case 5:
			return "1 ns BCAL w/ coplanarity + SC Veto";
		case 6:
			return "Loose BCAL+SC Veto";
		case 7:
			return "Strict BCAL+SC Veto";
		default:
			return "";
	}
}
