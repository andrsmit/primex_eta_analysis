TString folderName = "VetoOption7_CohLineshape_SGevorkyanFermi";

void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);
void StylePad(TPad *p);

void plotResults_Int() {
	
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>> phiVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	double Phi=0.0;
	while(Phi<=60.0) {
		TString filename = Form("%s/results_Int_%.1f.txt", folderName.Data(), Phi);
		if(gSystem->AccessPathName(filename.Data())) {
			Phi += 5.0;
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
		
		Phi += 5.0;
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
	
	TGraphErrors *gGamma = new TGraphErrors(nPoints, phi, gamma, phiErr, gammaErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, phi, coh,   phiErr,   cohErr);
	TGraphErrors *gInc   = new TGraphErrors(nPoints, phi, inc,   phiErr,   incErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, phi, chi2,  phiErr,  chi2Err);
	
	StyleGraph(gGamma, kBlue, "#Gamma(#eta#rightarrow#gamma#gamma) [keV]");
	StyleGraph(gCoh,   kBlue, "N.C. Normalization Factor");
	StyleGraph(gInc,   kBlue, "N.I. Normalization Factor");
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
	
	TCanvas *cInc = new TCanvas("cInc","Incoherent",950,600);
	TPad *p1Inc = new TPad("p1Inc","", 0, 0, 1, 1);
	TPad *p2Inc = new TPad("p2Inc","", 0, 0, 1, 1);
	p2Inc->SetFillStyle(4000);
	
	p1Inc->Draw();
	p1Inc->cd();
	gInc->Draw("APE");
	gPad->Update();
	
	xmin = p1Inc->GetUxmin();
	xmax = p1Inc->GetUxmax();
	dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	ymin = gChi2->GetHistogram()->GetMinimum();
	ymax = gChi2->GetHistogram()->GetMaximum();
	dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2Inc->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2Inc->Draw();
	p2Inc->cd();
	gChi2->Draw("PE");
	gPad->Update();
	
	TGaxis *axisInc = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axisInc->SetLineColor(kRed);
	axisInc->SetLabelColor(kRed);
	axisInc->SetTitleColor(kRed);
	axisInc->SetTitle("#chi^{2}/n.d.f.");
	axisInc->Draw();
	gPad->Update();
	
	//----------------------------------------------------------//
	
	//TCanvas *cChi2 = new TCanvas("cChi2","Chi2",950,600);
	//StyleCanvas(cChi2);
	//gChi2->Draw("APE");
	
	return;
}

void StyleGraph(TGraphErrors *g, int color, TString axisName) 
{
	g->SetTitle("");
	
	g->GetXaxis()->SetTitle("Interfernce Angle (fixed) [#circ]");
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
}

void StyleCanvas(TCanvas *c) 
{
	c->SetLeftMargin(0.13); c->SetRightMargin(0.13);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
}

void StylePad(TPad *c) 
{
	c->SetLeftMargin(0.13); c->SetRightMargin(0.13);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
}
