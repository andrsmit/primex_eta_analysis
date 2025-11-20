TString folderName = "VetoOption7_AngularYieldFit";

void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

void plotResults_NI() {
	
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>> phiVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	double NI=0.0;
	while(NI<2.1) {
		TString filename = Form("%s/results_NI_%.1f.txt", folderName.Data(), NI);
		if(gSystem->AccessPathName(filename.Data())) {
			NI += 0.1;
			continue;
		}
		ifstream inf(filename.Data());
		
		double locGamma, locGammaErr;
		double locPhi,   locPhiErr;
		double locCoh,   locCohErr;
		double locInc,   locIncErr;
		double locChi2;
		int    locNDF;
		
		inf >> locGamma >> locGammaErr;
		inf >> locPhi   >> locPhiErr;
		inf >> locCoh   >> locCohErr;
		inf >> locInc   >> locIncErr;
		inf >> locChi2  >> locNDF;
		inf.close();
		
		gammaVec.push_back({locGamma, locGammaErr});
		phiVec.push_back({locPhi, locPhiErr});
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
		
		NI += 0.1;
	}
	
	int nPoints = gammaVec.size();
	
	double *gamma    = new double[nPoints];
	double *gammaErr = new double[nPoints];
	double *phi      = new double[nPoints];
	double *phiErr   = new double[nPoints];
	double *coh      = new double[nPoints];
	double *cohErr   = new double[nPoints];
	double *inc      = new double[nPoints];
	double *incErr   = new double[nPoints];
	double *chi2     = new double[nPoints];
	double *chi2Err  = new double[nPoints];
	
	for(int i=0; i<nPoints; i++) {
		gamma[i]    = gammaVec[i].first;
		gammaErr[i] = gammaVec[i].second;
		phi[i]      = phiVec[i].first;
		phiErr[i]   = phiVec[i].second;
		coh[i]      = cohVec[i].first;
		cohErr[i]   = cohVec[i].second;
		inc[i]      = incVec[i].first;
		incErr[i]   = incVec[i].second;
		chi2[i]     = chi2Vec[i].first / ((double)chi2Vec[i].second);
		chi2Err[i]  = 0.0;
	}
	
	TGraphErrors *gGamma = new TGraphErrors(nPoints, inc, gamma, incErr, gammaErr);
	TGraphErrors *gPhi   = new TGraphErrors(nPoints, inc, phi,   incErr,   phiErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, inc, coh,   incErr,   cohErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, inc, chi2,  incErr,  chi2Err);
	
	StyleGraph(gGamma, kBlue, "#Gamma(#eta#rightarrow#gamma#gamma) [keV]");
	StyleGraph(gPhi,   kBlue, "Interference Angle [#circ]");
	StyleGraph(gCoh,   kBlue, "N.C. Normalization Factor");
	StyleGraph(gChi2,  kRed,  "#chi^{2} / n.d.f.");
	
	Double_t xmin, xmax, dx;
	Double_t ymin, ymax, dy;
	
	//----------------------------------------------------------//
	
	TCanvas *cGamma = new TCanvas("cGamma","Gamma",950,600);
	TPad *p1Gamma = new TPad("p1Gamma","", 0, 0, 1, 1);
	TPad *p2Gamma = new TPad("p2Gamma","", 0, 0, 1, 1);
	p2Gamma->SetFillStyle(4000);
	
	p1Gamma->Draw();
	p1Gamma->cd();
	gGamma->Draw("APE");
	gPad->Update();
	
	xmin = p1Gamma->GetUxmin();
	xmax = p1Gamma->GetUxmax();
	dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	ymin = gChi2->GetHistogram()->GetMinimum();
	ymax = gChi2->GetHistogram()->GetMaximum();
	dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2Gamma->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2Gamma->Draw();
	p2Gamma->cd();
	gChi2->Draw("PE");
	gPad->Update();
	
	TGaxis *axisGamma = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axisGamma->SetLineColor(kRed);
	axisGamma->SetLabelColor(kRed);
	axisGamma->SetTitleColor(kRed);
	axisGamma->SetTitle("#chi^{2}/n.d.f.");
	axisGamma->Draw();
	gPad->Update();
	
	//----------------------------------------------------------//
	
	TCanvas *cCoh = new TCanvas("cCoh","Coherent",950,600);
	TPad *p1Coh = new TPad("p1Coh","", 0, 0, 1, 1);
	TPad *p2Coh = new TPad("p2Coh","", 0, 0, 1, 1);
	p2Coh->SetFillStyle(4000);
	
	p1Coh->Draw();
	p1Coh->cd();
	gCoh->Draw("APE");
	gPad->Update();
	
	xmin = p1Coh->GetUxmin();
	xmax = p1Coh->GetUxmax();
	dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	ymin = gChi2->GetHistogram()->GetMinimum();
	ymax = gChi2->GetHistogram()->GetMaximum();
	dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2Coh->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2Coh->Draw();
	p2Coh->cd();
	gChi2->Draw("PE");
	gPad->Update();
	
	TGaxis *axisCoh = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axisCoh->SetLineColor(kRed);
	axisCoh->SetLabelColor(kRed);
	axisCoh->SetTitleColor(kRed);
	axisCoh->SetTitle("#chi^{2}/n.d.f.");
	axisCoh->Draw();
	gPad->Update();
	
	//----------------------------------------------------------//
	
	TCanvas *cPhi = new TCanvas("cPhi","Interference",950,600);
	TPad *p1Phi = new TPad("p1Phi","", 0, 0, 1, 1);
	TPad *p2Phi = new TPad("p2Phi","", 0, 0, 1, 1);
	p2Phi->SetFillStyle(4000);
	
	p1Phi->Draw();
	p1Phi->cd();
	gPhi->Draw("APE");
	gPad->Update();
	
	xmin = p1Phi->GetUxmin();
	xmax = p1Phi->GetUxmax();
	dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	ymin = gChi2->GetHistogram()->GetMinimum();
	ymax = gChi2->GetHistogram()->GetMaximum();
	dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2Phi->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2Phi->Draw();
	p2Phi->cd();
	gChi2->Draw("PE");
	gPad->Update();
	
	TGaxis *axisPhi = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axisPhi->SetLineColor(kRed);
	axisPhi->SetLabelColor(kRed);
	axisPhi->SetTitleColor(kRed);
	axisPhi->SetTitle("#chi^{2}/n.d.f.");
	axisPhi->Draw();
	gPad->Update();
	
	return;
}

void StyleGraph(TGraphErrors *g, int color, TString axisName) 
{
	g->SetTitle("");
	
	g->GetXaxis()->SetTitle("N.I. Normarlization Factor (fixed)");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(1.0);
	g->GetXaxis()->CenterTitle(true);
	
	g->GetYaxis()->SetTitle(axisName.Data());
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->CenterTitle(true);
	g->GetYaxis()->SetLabelColor(color);
	g->GetYaxis()->SetTitleColor(color);
	
	g->SetMarkerStyle(8);
	g->SetMarkerColor(color);
	g->SetLineColor(color);
	
	g->GetXaxis()->SetRangeUser(-0.5, 1.05);
}

void StyleCanvas(TCanvas *c) 
{
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
}
