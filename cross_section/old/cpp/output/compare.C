TString fileName1 = "yield_phase3_VetoOption4_1.5MeV.root";
TString fileName2 = "yield_phase3_VetoOption6_1.5MeV.root";

void compare()
{
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	
	TFile *f1 = new TFile(fileName1.Data(), "READ");
	TH1F  *h1 = (TH1F*)f1->Get("CrossSectionFit")->Clone("CrossSection1");
	TH1F  *h1b = (TH1F*)f1->Get("EtaPionYield")->Clone("CrossSection1b");
	//h1->Add(h1b);
	h1->SetDirectory(0);
	f1->Close();
	
	h1->SetLineWidth(1);
	h1->SetLineColor(kBlue);
	h1->SetMarkerColor(kBlue);
	h1->SetMarkerStyle(4);
	
	TFile *f2 = new TFile(fileName2.Data(), "READ");
	TH1F  *h2 = (TH1F*)f2->Get("CrossSectionFit")->Clone("CrossSection2");
	TH1F  *h2b = (TH1F*)f2->Get("EtaPionYield")->Clone("CrossSection2b");
	//h2->Add(h2b);
	h2->SetDirectory(0);
	f2->Close();
	
	h2->SetLineWidth(1);
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->SetMarkerStyle(8);
	
	c->cd();
	h1->Draw("PE1X0");
	h2->Draw("PE1X0 same");
	
	TLegend *leg1 = new TLegend(0.15, 0.7, 0.45, 0.9);
	leg1->AddEntry(h1, "Loose BCAL Veto", "PE1X0");
	leg1->AddEntry(h2, "Loose BCAL+SC Veto", "PE1X0");
	leg1->Draw();
	
	return;
}
