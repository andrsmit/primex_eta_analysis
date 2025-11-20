TString GetVetoString(int option);

void compareTOF()
{
	gStyle->SetOptStat(0);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	c->SetGrid();
	
	TLegend *leg1 = new TLegend(0.146, 0.703, 0.514, 0.963);
	
	vector<TString> TOFOptions = {"noTOF", "TOF", "singleTOF"};
	vector<TString>   TOFNames = {"No TOF Veto", "Loose TOF Veto", "Strict TOF Veto"};
	vector<int>   colorOptions = {kMagenta,  kBlue,        kRed};
	vector<int>  markerOptions = {     22,     25,          8};
	
	TH1F *h0;
	bool first = true;
	
	for(int i=0; i<TOFOptions.size(); i++) {
		TString locTOFOption = TOFOptions[i];
		
		TFile *fIn = new TFile(Form("yield_phase3_VetoOption6_%s.root", locTOFOption.Data()), "READ");
		
		TH1F *hc  = (TH1F*)fIn->Get("CrossSectionFit")->Clone(Form("CrossSection_%s",locTOFOption.Data()));
		hc->SetDirectory(0);
		
		hc->SetLineWidth(1);
		hc->SetLineColor(colorOptions[i]);
		hc->SetMarkerColor(colorOptions[i]);
		hc->SetMarkerStyle(markerOptions[i]);
		hc->SetMarkerSize(0.9);
		
		if(first) {
			h0 = (TH1F*)hc->Clone("h0");
			first = false;
		}
		
		TH1F *hr = (TH1F*)hc->Clone(Form("Ratio_%s",locTOFOption.Data()));
		hr->Divide(h0);
		
		hc->GetYaxis()->SetRangeUser(0.0,2.1);
		
		if(i==0) {
			hc->Draw("PE1X0");
		}
		else {
			hc->Draw("PE1X0 same");
		}
		leg1->AddEntry(hc, Form("%s",TOFNames[i].Data()), "PE1X0");
	}
	leg1->Draw();
	
	return;
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
			return "Modified BCAL Veto";
		case 5:
			return "1 ns BCAL w/ coplanarity + SC Veto";
		case 6:
			return "Modified BCAL+SC Veto";
		case 7:
			return "Strict BCAL+SC Veto";
		default:
			return "";
	}
}
