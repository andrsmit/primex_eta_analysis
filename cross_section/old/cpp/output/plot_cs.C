void plot_cs(int vetoOption=4)
{
	gStyle->SetOptStat(0);
	
	//TFile *f1 = new TFile(Form("testingShifts/yield_phase3_VetoOption%d_floatingShift.root",2), "READ");
	TFile *f1 = new TFile(Form("yield_phase3_VetoOption%d_.root",vetoOption), "READ");
	
	
	TString igString = "";
	if(vetoOption==1) igString = "strict-bcal";
	else if(vetoOption==2) igString = "loose-bcal-and-no-coplanarity-cuts";
	else if(vetoOption==4) igString = "loose-bcal-and-coplanarity-cuts";
	else if(vetoOption==7) igString = "strict-bcal-and-sc";
	else igString = "loose-bcal-and-coplanarity-cuts";
	
	
	TFile *fIgal = new TFile(Form("/work/halld/home/ijaegle/public/ForDrew/preliminary-%s-2022-08.root", igString.Data()), "READ");
	TGraphErrors *gIgal = (TGraphErrors*)fIgal->Get("gr_xs")->Clone("gIgal");
	
	TH1F *h1c  = (TH1F*)f1->Get("CrossSectionFit")->Clone("h1c");
	TH1F *h1y  = (TH1F*)f1->Get("AngularYieldFit")->Clone("h1y");
	TH1F *h1b1 = (TH1F*)f1->Get("EtaPionYield")->Clone("h1b1");
	TH1F *h1b2 = (TH1F*)f1->Get("HadronicBkgdYield")->Clone("h1b2");
	TH1F *h1b3 = (TH1F*)f1->Get("BkgdYield")->Clone("h1b3");
	
	TH1F *h1cb   = (TH1F*)f1->Get("CrossSectionFit")->Clone("h1cb");
	TH1F *h1cSub = (TH1F*)f1->Get("CrossSectionFit")->Clone("h1cSub");
	
	TH1F *h1ci = (TH1F*)f1->Get("CrossSectionInclusive")->Clone("h1ci");
	h1ci->SetLineColor(kBlue);
	h1ci->SetMarkerColor(kBlue);
	
	TH1F *h1co = (TH1F*)f1->Get("Counts")->Clone("h1co");
	TH1F *h1eco = (TH1F*)f1->Get("EmptyCounts")->Clone("h1eco");
	h1co->Add(h1eco,-1.0);
	
	TH1F *h1c_co  = (TH1F*)h1c->Clone("h1c_co");
	h1c_co->SetMarkerColor(kBlack);
	h1c_co->SetLineColor(kBlack);
	/*
	for(int ib=1; ib<h1c->GetXaxis()->GetNbins(); ib++) {
		double y  = h1y->GetBinContent(ib);
		double b1 = h1b1->GetBinContent(ib);
		double b2 = h1b2->GetBinContent(ib);
		double b3 = h1b3->GetBinContent(ib);
		double co = h1co->GetBinContent(ib);
		double c  = h1c->GetBinContent(ib);
		h1c_co->SetBinContent(ib, c*(co/y));
		h1c->SetBinContent(ib,c*(y+b1+b2+b3)/y);
		h1cb->SetBinContent(ib,c*(b1+b2+b3)/y);
	}
	*/
	h1c->SetLineColor(kBlue);
	h1c->SetMarkerColor(kBlue);
	h1c->SetMarkerStyle(8);
	
	h1c_co->GetYaxis()->SetRangeUser(0.0,1.2*h1c_co->GetMaximum());
	h1c->GetYaxis()->SetRangeUser(0.0,1.2*h1c->GetMaximum());
	
	h1c->SetMarkerColor(kBlack);
	h1c->SetLineColor(kBlack);
	
	
	gIgal->SetMarkerColor(kRed);
	gIgal->SetLineColor(kRed);
	gIgal->SetMarkerStyle(8);
	
	
	TCanvas *c1 = new TCanvas("c1","c1",950,700);
	c1->SetLeftMargin(0.13);
	c1->SetRightMargin(0.07);
	c1->SetTopMargin(0.07);
	c1->SetBottomMargin(0.13);
	c1->cd();
	
	//gIgal->Draw("APE");
	h1c->Draw("PE1X0");
	h1ci->Draw("PE1X0 same");
	
	return;
}
