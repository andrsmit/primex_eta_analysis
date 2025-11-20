void plot_energy_dependence() {
	
	double energy[3] = {8.5, 9.5, 10.65};
	double energyE[3] = {0.5, 0.5, 0.65};
	
	double avg_gamma = 0.438655;
	double gamma[3] = {0.411684, 0.467358, 0.418881};
	double gammaE[3] = {0.0218696, 0.029355, 0.0176076};
	
	double avg_coh = 0.681597;
	double coh[3] = {0.566996, 0.665491, 0.769842};
	double cohE[3] = {0.0687378, 0.0683878, 0.0581345};
	
	double avg_inc = 0.433723;
	double inc[3] = {0.397115, 0.45854, 0.463012};
	double incE[3] = {0.0445649, 0.0358378, 0.0247029};
	
	for(int i=0; i<3; i++) {
		gamma[i]  = gamma[i]/avg_gamma;
		gammaE[i] = gammaE[i]/avg_gamma;
		
		coh[i]  = coh[i]/avg_coh;
		cohE[i] = cohE[i]/avg_coh;
		
		inc[i]  = inc[i]/avg_inc;
		incE[i] = incE[i]/avg_inc;
	}
	
	TGraphErrors *gGamma = new TGraphErrors(3, energy, gamma, energyE, gammaE);
	gGamma->GetXaxis()->SetTitle("Photon Beam Energy Range [GeV]");
	gGamma->GetXaxis()->SetTitleSize(0.05);
	gGamma->GetXaxis()->SetTitleOffset(1.0);
	gGamma->GetXaxis()->CenterTitle(true);
	gGamma->GetYaxis()->SetTitle("Fitted #Gamma(#eta#rightarrow#gamma#gamma) / Fit to Full Energy Range");
	gGamma->GetYaxis()->SetTitleSize(0.05);
	gGamma->GetYaxis()->SetTitleOffset(1.0);
	gGamma->GetYaxis()->CenterTitle(true);
	gGamma->SetTitle("");
	gGamma->SetMarkerStyle(8);
	gGamma->SetMarkerSize(1.0);
	gGamma->SetMarkerColor(kBlue);
	gGamma->SetLineColor(kBlue);
	gGamma->SetLineWidth(2);
	gGamma->GetYaxis()->SetRangeUser(0.76,1.24);
	
	TCanvas *cGamma = new TCanvas("cGamma", "cGamma", 1000, 500);
	cGamma->SetLeftMargin(0.13); cGamma->SetRightMargin(0.07);
	cGamma->SetBottomMargin(0.13); cGamma->SetTopMargin(0.07);
	gGamma->Draw("AP");
	
	//--------------------------------------//
	
	TGraphErrors *gCoh = new TGraphErrors(3, energy, coh, energyE, cohE);
	gCoh->GetXaxis()->SetTitle("Photon Beam Energy Range [GeV]");
	gCoh->GetXaxis()->SetTitleSize(0.05);
	gCoh->GetXaxis()->SetTitleOffset(1.0);
	gCoh->GetXaxis()->CenterTitle(true);
	gCoh->GetYaxis()->SetTitle("Fitted A_{NC} / Fit to Full Energy Range");
	gCoh->GetYaxis()->SetTitleSize(0.05);
	gCoh->GetYaxis()->SetTitleOffset(1.0);
	gCoh->GetYaxis()->CenterTitle(true);
	gCoh->SetTitle("");
	gCoh->SetMarkerStyle(8);
	gCoh->SetMarkerSize(1.0);
	gCoh->SetMarkerColor(kBlue);
	gCoh->SetLineColor(kBlue);
	gCoh->SetLineWidth(2);
	gCoh->GetYaxis()->SetRangeUser(0.76,1.24);
	
	TCanvas *cCoh = new TCanvas("cCoh", "cCoh", 1000, 500);
	cCoh->SetLeftMargin(0.13); cCoh->SetRightMargin(0.07);
	cCoh->SetBottomMargin(0.13); cCoh->SetTopMargin(0.07);
	gCoh->Draw("AP");
	
	//--------------------------------------//
	
	TGraphErrors *gInc = new TGraphErrors(3, energy, inc, energyE, incE);
	gInc->GetXaxis()->SetTitle("Photon Beam Energy Range [GeV]");
	gInc->GetXaxis()->SetTitleSize(0.05);
	gInc->GetXaxis()->SetTitleOffset(1.0);
	gInc->GetXaxis()->CenterTitle(true);
	gInc->GetYaxis()->SetTitle("Fitted A_{NI} / Fit to Full Energy Range");
	gInc->GetYaxis()->SetTitleSize(0.05);
	gInc->GetYaxis()->SetTitleOffset(1.0);
	gInc->GetYaxis()->CenterTitle(true);
	gInc->SetTitle("");
	gInc->SetMarkerStyle(8);
	gInc->SetMarkerSize(1.0);
	gInc->SetMarkerColor(kBlue);
	gInc->SetLineColor(kBlue);
	gInc->SetLineWidth(2);
	gInc->GetYaxis()->SetRangeUser(0.76,1.24);
	
	TCanvas *cInc = new TCanvas("cInc", "cInc", 1000, 500);
	cInc->SetLeftMargin(0.13); cInc->SetRightMargin(0.07);
	cInc->SetBottomMargin(0.13); cInc->SetTopMargin(0.07);
	gInc->Draw("AP");
	
	return;
}
