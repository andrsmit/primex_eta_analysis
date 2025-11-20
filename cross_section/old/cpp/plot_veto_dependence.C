void plot_veto_dependence() {
	
	double option[4] = {1, 2, 3, 4};
	double optionE[4] = {0, 0, 0, 0};
	
	/*
	double avg_gamma_sub = 0.0;
	double gamma_sub[4] = {0.418410, 0.423691, 0.454536, 0.431751};
	double gammaE_sub[4] = {0.0209677, 0.0208966, 0.0181942, 0.0165401};
	
	double avg_gamma_noSub = 0.0;
	double gamma_noSub[4] = {0.441248, 0.443424, 0.461281, 0.438486};
	double gammaE_noSub[4] = {0.0210331, 0.0132748, 0.0181968, 0.0165890};
	*/
	double avg_gamma_sub = 0.0;
	double gamma_sub[4] = {0.4215, 0.433, 0.4272, 0.435};
	double gammaE_sub[4] = {0.02, 0.02, 0.02, 0.02};
	
	double avg_gamma_noSub = 0.0;
	double gamma_noSub[4] = {0.4416, 0.4455, 0.4441, 0.439};
	double gammaE_noSub[4] = {0.02, 0.02, 0.02, 0.02};
	
	
	for(int i=0; i<4; i++) {
		avg_gamma_sub += gamma_sub[i];
		avg_gamma_noSub += gamma_noSub[i];
	}
	avg_gamma_sub /= 4.0;
	avg_gamma_noSub /= 4.0;
	
	avg_gamma_sub = gamma_sub[0];
	
	for(int i=0; i<4; i++) {
		gamma_sub[i]  = gamma_sub[i]/avg_gamma_sub;
		gammaE_sub[i] = gammaE_sub[i]/avg_gamma_sub;
		
		gamma_noSub[i]  = gamma_noSub[i]/avg_gamma_sub;
		gammaE_noSub[i] = gammaE_noSub[i]/avg_gamma_sub;
	}
	
	TGraphErrors *gGammaSub = new TGraphErrors(4, option, gamma_sub, optionE, gammaE_sub);
	gGammaSub->GetXaxis()->SetTitle("");
	gGammaSub->GetXaxis()->SetTitleSize(0.05);
	gGammaSub->GetXaxis()->SetLabelSize(0.045);
	gGammaSub->GetXaxis()->SetTitleOffset(1.0);
	gGammaSub->GetXaxis()->CenterTitle(true);
	gGammaSub->GetYaxis()->SetTitle("Fitted #Gamma(#eta#rightarrow#gamma#gamma) / First Point");
	gGammaSub->GetYaxis()->SetTitleSize(0.05);
	gGammaSub->GetYaxis()->SetTitleOffset(1.0);
	gGammaSub->GetYaxis()->CenterTitle(true);
	gGammaSub->SetTitle("");
	gGammaSub->SetMarkerStyle(8);
	gGammaSub->SetMarkerSize(1.0);
	gGammaSub->SetMarkerColor(kBlue);
	gGammaSub->SetLineColor(kBlue);
	gGammaSub->SetLineWidth(2);
	gGammaSub->GetYaxis()->SetRangeUser(0.76,1.24);
	
	gGammaSub->GetXaxis()->SetBinLabel(gGammaSub->GetXaxis()->FindBin(1), "Strict BCAL Veto");
	gGammaSub->GetXaxis()->SetBinLabel(gGammaSub->GetXaxis()->FindBin(2), "Loose BCAL Veto");
	gGammaSub->GetXaxis()->SetBinLabel(gGammaSub->GetXaxis()->FindBin(3), "Loose BCAL+SC Veto");
	gGammaSub->GetXaxis()->SetBinLabel(gGammaSub->GetXaxis()->FindBin(4), "Strict BCAL+SC Veto");
	
	TGraphErrors *gGammaNoSub = new TGraphErrors(4, option, gamma_noSub, optionE, gammaE_noSub);
	gGammaNoSub->GetXaxis()->SetTitle("Photon Beam Energy Range [GeV]");
	gGammaNoSub->GetXaxis()->SetTitleSize(0.05);
	gGammaNoSub->GetXaxis()->SetTitleOffset(1.0);
	gGammaNoSub->GetXaxis()->CenterTitle(true);
	gGammaNoSub->GetYaxis()->SetTitle("Fitted #Gamma(#eta#rightarrow#gamma#gamma) / Fit to Full Energy Range");
	gGammaNoSub->GetYaxis()->SetTitleSize(0.05);
	gGammaNoSub->GetYaxis()->SetTitleOffset(1.0);
	gGammaNoSub->GetYaxis()->CenterTitle(true);
	gGammaNoSub->SetTitle("");
	gGammaNoSub->SetMarkerStyle(8);
	gGammaNoSub->SetMarkerSize(1.0);
	gGammaNoSub->SetMarkerColor(kRed);
	gGammaNoSub->SetLineColor(kRed);
	gGammaNoSub->SetLineWidth(2);
	gGammaNoSub->GetYaxis()->SetRangeUser(0.76,1.24);
	
	TCanvas *cGammaSub = new TCanvas("cGammaSub", "cGammaSub", 1000, 500);
	cGammaSub->SetLeftMargin(0.13); cGammaSub->SetRightMargin(0.07);
	cGammaSub->SetBottomMargin(0.13); cGammaSub->SetTopMargin(0.07);
	gGammaSub->Draw("APE");
	gGammaNoSub->Draw("PE same");
	
	TLegend *leg = new TLegend(0.1, 0.7, 0.4, 0.9);
	leg->AddEntry(gGammaSub, "Had. Bkgd Subtracted Before Fit", "PE");
	leg->AddEntry(gGammaNoSub, "Had. Bkgd Not Subtracted", "PE");
	leg->Draw();
	
	return;
}
