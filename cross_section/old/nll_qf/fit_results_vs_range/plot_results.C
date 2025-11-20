void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

void plot_results()
{
	
	vector<double> angleVec;
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>> phiVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	double locAngle=0.5;
	while(locAngle<4.5) {
		TString filename = Form("results_%.2f.txt", locAngle);
		if(gSystem->AccessPathName(filename.Data())) {
			locAngle += 0.25;
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
		
		angleVec.push_back(locAngle);
		gammaVec.push_back({locGamma, locGammaErr});
		phiVec.push_back({locPhi, locPhiErr});
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
		
		locAngle += 0.25;
	}
	
	int nPoints = angleVec.size();
	
	double *angle    = new double[nPoints];
	double *angleErr = new double[nPoints];
	
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
		angle[i]    = angleVec[i];
		angleErr[i] = 0.0;
		gamma[i]    = gammaVec[i].first;
		gammaErr[i] = gammaVec[i].second;
		phi[i]      = phiVec[i].first;
		phiErr[i]   = phiVec[i].second;
		coh[i]      = cohVec[i].first;
		cohErr[i]   = cohVec[i].second;
		inc[i]      = incVec[i].first;
		incErr[i]   = incVec[i].second;
		chi2[i]     = chi2Vec[i].first   / ((double)chi2Vec[i].second);
		chi2Err[i]  = 0.0;
	}
	/*
	double gammaSave = gamma[0];
	for(int i=0; i<nPoints; i++) {
		gamma[i]    /= gammaSave;
		gammaErr[i] /= gammaSave;
	}
	*/
	TGraphErrors *gGamma = new TGraphErrors(nPoints, angle, gamma, angleErr, gammaErr);
	TGraphErrors *gPhi   = new TGraphErrors(nPoints, angle, phi,   angleErr,   phiErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, angle, coh,   angleErr,   cohErr);
	TGraphErrors *gInc   = new TGraphErrors(nPoints, angle, inc,   angleErr,   incErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, angle, chi2,  angleErr,  chi2Err);
	
	StyleGraph(gGamma, kRed, "#Gamma(#eta#rightarrow#gamma#gamma) [keV]");
	StyleGraph(gPhi,   kBlue, "Interference Angle [#circ]");
	StyleGraph(gCoh,   kBlack, "N.C. Normalization Factor");
	StyleGraph(gInc,   kGreen, "N.I. Normalization Factor");
	StyleGraph(gChi2,  kCyan,  "#chi^{2} / n.d.f.");
	
	TCanvas *cResults = new TCanvas("cResults","Results",950,600);
	StyleCanvas(cResults);
	gGamma->Draw("AP");
	//gPhi->Draw("P same");
	
	TCanvas *cResults2 = new TCanvas("cResults2","Coh",950,600);
	StyleCanvas(cResults2);
	gCoh->Draw("AP");
	//gInc->Draw("P same");
	//gChi2->Draw("P same");
	
	return;
}

void StyleGraph(TGraphErrors *g, int color, TString axisName) 
{
	g->SetTitle("");
	
	g->GetXaxis()->SetTitle("Angular Fit Range [#circ]");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(1.0);
	g->GetXaxis()->CenterTitle(true);
	
	g->GetYaxis()->SetTitle(axisName.Data());
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->CenterTitle(true);
	//g->GetYaxis()->SetLabelColor(color);
	//g->GetYaxis()->SetTitleColor(color);
	
	g->SetMarkerStyle(8);
	g->SetMarkerColor(color);
	g->SetLineColor(color);
}

void StyleCanvas(TCanvas *c) 
{
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
}
