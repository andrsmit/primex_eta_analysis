//TString fileName1 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption2.root";
//TString fileName2 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";

//TString fileName1 = "fixedFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";
//TString fileName2 = "floatingFraction_bkgdShifted_2MeV/yield_phase3_VetoOption6.root";

vector<pair<double,double>> energies = {{8.0,9.0}, {9.0,10.0}, {10.0,11.3}};

vector<int>   lineColors = {kBlack, kBlue, kRed, kMagenta};
vector<int> markerStyles = {4, 20, 23, 25};

void StyleHist(TH1F *h1, int lineColor=kBlack, int markerStyle=4);
void StyleCanvas(TCanvas*);

void compare_energies()
{
	gStyle->SetOptStat(0);
	
	int n_energies = (int)energies.size();
	
	TCanvas *c = new TCanvas("c","c",950,700);
	StyleCanvas(c);
	
	TLegend *leg = new TLegend(0.150, 0.70, 0.445, 0.90);
	
	bool first = true;
	for(int ien=0; ien<n_energies; ien++) {
		
		double locMinE = energies[ien].first;
		double locMaxE = energies[ien].second;
		
		TString fileName = Form("yield_phase3_VetoOption6_%.1fGeV_%.1fGeV.root", locMinE, locMaxE);
		if(gSystem->AccessPathName(fileName.Data())) continue;
		
		TFile *fIn = new TFile(fileName.Data(), "READ");
		
		TH1F *h1 = (TH1F*)fIn->Get("CrossSectionFit")->Clone(Form("cs_%.1fGeV_%.1fGeV",locMinE,locMaxE));
		h1->SetDirectory(0);
		fIn->Close();
		
		int locLineColor;
		if(ien>lineColors.size()) locLineColor = lineColors[lineColors.size()-1] + 1;
		else locLineColor = lineColors[ien];
		
		int locMarkerStyle;
		if(ien>markerStyles.size()) locMarkerStyle = markerStyles[markerStyles.size()-1] + 1;
		else locMarkerStyle = markerStyles[ien];
		
		StyleHist(h1, locLineColor, locMarkerStyle);
		
		TString drawOption = first ? "PE1X0" : "PE1X0 same";
		c->cd();
		h1->Draw(drawOption.Data());
		
		leg->AddEntry(h1, Form("%.1f GeV < E_{#gamma} < %.1f GeV", locMinE, locMaxE), "PE1X0");
		
		first = false;
	}
	leg->Draw();
	
	return;
}

void StyleHist(TH1F *h1, int lineColor, int markerStyle)
{
	h1->GetYaxis()->SetRangeUser(0.0,2.0);
	
	h1->SetLineColor(lineColor);
	h1->SetMarkerColor(lineColor);
	h1->SetMarkerStyle(markerStyle);
	h1->SetLineWidth(2);
	
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetXaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetMaxDigits(2);
	
	for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
		double locC = h1->GetBinContent(ibin);
		double locE = h1->GetBinError(ibin);
		if(locE > 0.25*locC) h1->SetBinError(ibin, 0.25*locC);
	}
}

void StyleCanvas(TCanvas *c1)
{
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	c1->SetTickx(); c1->SetTicky();
	c1->cd();
}
