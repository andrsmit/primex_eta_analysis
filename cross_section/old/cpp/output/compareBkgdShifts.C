void compareBkgdShifts(int vetoOption=4)
{
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	
	TFile *f1 = new TFile(Form("bkgdShifted_0.000/yield_phase3_VetoOption%d.root", vetoOption), "READ");
	TH1F  *h1 = (TH1F*)f1->Get("CrossSection")->Clone(Form("CrossSection1_VetoOption%d",vetoOption));
	h1->SetDirectory(0);
	f1->Close();
	
	h1->SetLineWidth(1);
	h1->SetLineColor(kBlue);
	h1->SetMarkerColor(kBlue);
	h1->SetMarkerStyle(4);
	
	TFile *f2 = new TFile(Form("bkgdShifted_0.006/yield_phase3_VetoOption%d.root", vetoOption), "READ");
	TH1F  *h2 = (TH1F*)f2->Get("CrossSection")->Clone(Form("CrossSection2_VetoOption%d",vetoOption));
	h2->SetDirectory(0);
	f2->Close();
	
	h2->SetLineWidth(1);
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->SetMarkerStyle(8);
	
	c->cd();
	h1->Draw("PE1X0");
	h2->Draw("PE1X0 same");
	
	return;
}
