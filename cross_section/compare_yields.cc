#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/style.cc"

void compare_yields(TString yield1_fname="output1.root", TString yield2_fname="output2.root") {
	
	//-------------------------------------------------------------------------//
	// edit these to give a description of the two data sets for the plot legend:
	
	TString leg_title_1 = "Yield1";
	TString leg_title_2 = "Yield2";
	
	leg_title_1 = "Strict BCAL Veto";
	leg_title_2 = "Loose BCAL Veto";
	
	//-------------------------------------------------------------------------//
	
	if(gSystem->AccessPathName(yield1_fname.Data()) || gSystem->AccessPathName(yield2_fname.Data())) {
		cout << "Input files not found." << endl;
		return;
	}
	
	TFile *fIn1 = new TFile(yield1_fname.Data(), "READ");
	TFile *fIn2 = new TFile(yield2_fname.Data(), "READ");
	
	TGraphErrors *g1 = (TGraphErrors*)fIn1->Get("eta_gg_yield")->Clone("yield1");
	TGraphErrors *g2 = (TGraphErrors*)fIn2->Get("eta_gg_yield")->Clone("yield2");
	
	TCanvas *c = new TCanvas("c","c",700,500);
	styleCanvas(c);
	
	// find the maximum value amongst the two graphs:
	int max_graph    = 1;
	double max_yield = 0.;
	for(int ibin=0; ibin<g1->GetN(); ibin++) {
		if(g1->GetY()[ibin] > max_yield) {
			max_yield = g1->GetY()[ibin];
		}
	}
	for(int ibin=0; ibin<g2->GetN(); ibin++) {
		if(g2->GetY()[ibin] > max_yield) {
			max_yield = g2->GetY()[ibin];
			max_graph = 2;
		}
	}
	
	g1->GetYaxis()->SetRangeUser(0.0, 1.1*max_yield);
	g2->GetYaxis()->SetRangeUser(0.0, 1.1*max_yield);
	
	g2->SetMarkerColor(kTeal-6);
	g2->SetLineColor(kTeal-6);
	
	c->cd();
	
	g1->Draw("AP");
	g2->Draw("P same");
	
	if(g1->GetN() != g2->GetN()) return;
	
	double *g1_angle     = g1->GetX();
	double *g1_angle_err = g1->GetEX();
	
	double *g1_yield     = g1->GetY();
	double *g1_yield_err = g1->GetEY();
	double *g2_yield     = g2->GetY();
	double *g2_yield_err = g2->GetEY();
	
	TGraphErrors *g3;
	if(max_graph==1) {
		g3 = (TGraphErrors*)g1->Clone("yield_difference");
		for(int ibin=0; ibin<g3->GetN(); ibin++) {
			double loc_yield_err = sqrt(pow(g1_yield_err[ibin],2.0) + pow(g2_yield_err[ibin],2.0));
			g3->SetPointY(ibin, g1_yield[ibin]-g2_yield[ibin]);
			g3->SetPointError(ibin, g1_angle_err[ibin], loc_yield_err);
		}
	} else {
		g3 = (TGraphErrors*)g1->Clone("yield_difference");
		for(int ibin=0; ibin<g3->GetN(); ibin++) {
			double loc_yield_err = sqrt(pow(g1_yield_err[ibin],2.0) + pow(g2_yield_err[ibin],2.0));
			g3->SetPointY(ibin, g2_yield[ibin]-g1_yield[ibin]);
			g3->SetPointError(ibin, g1_angle_err[ibin], loc_yield_err);
		}
	}
	g3->SetLineColor(kRed);
	g3->SetMarkerColor(kRed);
	
	c->cd();
	g3->Draw("P same");
	
	TLegend *leg = new TLegend(0.6, 0.6, 0.8, 0.8);
	leg->AddEntry(g1, leg_title_1, "PE");
	leg->AddEntry(g2, leg_title_2, "PE");
	leg->AddEntry(g3, "Difference", "PE");
	leg->Draw();
	
	return;
}
