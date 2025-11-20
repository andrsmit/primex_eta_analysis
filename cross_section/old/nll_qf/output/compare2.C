//TString fileName1 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption2.root";
//TString fileName2 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";

TString fileName1 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";
TString fileName2 = "floatingFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";

void compare2()
{
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	
	TFile *f1  = new TFile(fileName1.Data(), "READ");
	TH1F  *h1  = (TH1F*)f1->Get("CrossSectionFit")->Clone("CrossSection1");
	TH1F  *h1y = (TH1F*)f1->Get("AngularYieldFit")->Clone("AngularYield1");
	h1->SetDirectory(0);
	h1y->SetDirectory(0);
	f1->Close();
	
	h1->SetLineWidth(1);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetMarkerStyle(4);
	
	for(int ibin=1; ibin<=h1y->GetXaxis()->GetNbins(); ibin++) {
		double locYield = h1y->GetBinContent(ibin);
		double locErr   = h1y->GetBinError(ibin);
		if(locErr > (2.5*sqrt(locYield))) locErr = 2.5*sqrt(locYield);
		double locRelErr = locErr/locYield;
		h1->SetBinError(ibin, h1->GetBinContent(ibin)*locRelErr);
	}
	
	
	TFile *f2  = new TFile(fileName2.Data(), "READ");
	TH1F  *h2  = (TH1F*)f2->Get("CrossSectionFit")->Clone("CrossSection2");
	TH1F  *h2y = (TH1F*)f2->Get("AngularYieldFit")->Clone("AngularYield2");
	h2->SetDirectory(0);
	h2y->SetDirectory(0);
	f2->Close();
	
	h2->SetLineWidth(1);
	h2->SetLineColor(kBlue);
	h2->SetMarkerColor(kBlue);
	h2->SetMarkerStyle(8);
	
	for(int ibin=1; ibin<=h2y->GetXaxis()->GetNbins(); ibin++) {
		double locYield = h2y->GetBinContent(ibin);
		double locErr   = h2y->GetBinError(ibin);
		if(locErr > (2.5*sqrt(locYield))) locErr = 2.5*sqrt(locYield);
		double locRelErr = locErr/locYield;
		h2->SetBinError(ibin, h2->GetBinContent(ibin)*locRelErr);
	}
	
	c->cd();
	h1->Draw("PE1X0");
	h2->Draw("PE1X0 same");
	
	TLegend *leg1 = new TLegend(0.15, 0.7, 0.45, 0.9);
	leg1->AddEntry(h1, "z_{qf} Floating",     "PE1X0");
	leg1->AddEntry(h2, "z_{qf} Fixed to 1.0", "PE1X0");
	leg1->Draw();
	
	TH1F *hDiff1 = (TH1F*)h2->Clone("hDiff1");
	hDiff1->GetYaxis()->SetTitle("Relative Difference");
	hDiff1->GetYaxis()->SetRangeUser(0.0, 1.0);
	for(int ibin=1; ibin<=hDiff1->GetXaxis()->GetNbins(); ibin++) {
		double c1 = h1->GetBinContent(ibin);
		double e1 = h1->GetBinError(ibin);
		double c2 = h2->GetBinContent(ibin);
		double e2 = h2->GetBinError(ibin);
		
		double locDiff, locDiffErr;
		if((c1+c2)<=0.0) {
			locDiff    = 0.;
			locDiffErr = 0.0;
		} else {
			locDiff    = (c1-c2)/(0.5*(c1+c2));
			locDiffErr = e1/(0.5*(c1+c2));
		}
		hDiff1->SetBinContent(ibin, locDiff);
		hDiff1->SetBinError(ibin, locDiffErr);
	}
	
	TCanvas *cDiff = new TCanvas("cDiff","cDiff",950,700);
	cDiff->SetLeftMargin(0.13); cDiff->SetRightMargin(0.07);
	cDiff->SetBottomMargin(0.13); cDiff->SetTopMargin(0.07);
	cDiff->cd();
	hDiff1->Draw("PE1X0");
	
	return;
}
