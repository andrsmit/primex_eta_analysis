void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

void plot_results()
{
	
	vector<pair<double,double>> enVec_pair = {{8.0,9.0}, {9.0,10.0}, {10.0,11.3}};
	
	vector<double> enVec;
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>> phiVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	double locAngle=0.5;
	for(int ie=0; ie<enVec_pair.size(); ie++) {
		
		double locMinEn = enVec_pair[ie].first;
		double locMaxEn = enVec_pair[ie].second;
		double locEn    = 0.5*(locMinEn+locMaxEn);
		
		TString filename = Form("results_%.1fGeV_%.1fGeV.txt", locMinEn, locMaxEn);
		if(gSystem->AccessPathName(filename.Data())) {
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
		
		enVec.push_back(locEn);
		gammaVec.push_back({locGamma, locGammaErr});
		phiVec.push_back({locPhi, locPhiErr});
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
	}
	
	int nPoints = enVec.size();
	
	double *energy    = new double[nPoints];
	double *energyErr = new double[nPoints];
	
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
		energy[i]   = enVec[i];
		energyErr[i]= 0.5*(enVec_pair[i].second - enVec_pair[i].first);
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
	TGraphErrors *gGamma = new TGraphErrors(nPoints, energy, gamma, energyErr, gammaErr);
	TGraphErrors *gPhi   = new TGraphErrors(nPoints, energy, phi,   energyErr,   phiErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, energy, coh,   energyErr,   cohErr);
	TGraphErrors *gInc   = new TGraphErrors(nPoints, energy, inc,   energyErr,   incErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, energy, chi2,  energyErr,  chi2Err);
	
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
	
	g->GetXaxis()->SetTitle("Photon Beam Energy Range [GeV]");
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
