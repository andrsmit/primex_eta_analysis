void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

void plot_results()
{
	
	vector<double> smearVec;
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>> phiVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	int locSmear=0;
	while(locSmear<11) {
		TString filename = Form("smear_%02d.txt", locSmear);
		if(gSystem->AccessPathName(filename.Data())) {
			locSmear += 1;
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
		
		smearVec.push_back(0.5*(double)(locSmear));
		gammaVec.push_back({locGamma, locGammaErr});
		phiVec.push_back({locPhi, locPhiErr});
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
		
		locSmear += 1;
	}
	
	int nPoints = smearVec.size();
	
	double *smear    = new double[nPoints];
	double *smearErr = new double[nPoints];
	
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
		smear[i]    = smearVec[i];
		smearErr[i] = 0.0;
		gamma[i]    = gammaVec[i].first;
		gammaErr[i] = gammaVec[i].second;
		phi[i]      = phiVec[i].first    / phiVec[nPoints-1].first;
		phiErr[i]   = phiVec[i].second   / phiVec[nPoints-1].first;
		coh[i]      = cohVec[i].first    / cohVec[nPoints-1].first;
		cohErr[i]   = cohVec[i].second   / cohVec[nPoints-1].first;
		inc[i]      = incVec[i].first    / incVec[nPoints-1].first;
		incErr[i]   = incVec[i].second   / incVec[nPoints-1].first;
		chi2[i]     = chi2Vec[i].first   / ((double)chi2Vec[i].second);
		chi2Err[i]  = 0.0;
	}
	
	double gammaSave = gamma[0];
	for(int i=0; i<nPoints; i++) {
		gamma[i]    /= gammaSave;
		gammaErr[i] /= gammaSave;
	}
	
	TGraphErrors *gGamma = new TGraphErrors(nPoints, smear, gamma, smearErr, gammaErr);
	TGraphErrors *gPhi   = new TGraphErrors(nPoints, smear, phi,   smearErr,   phiErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, smear, coh,   smearErr,   cohErr);
	TGraphErrors *gInc   = new TGraphErrors(nPoints, smear, inc,   smearErr,   incErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, smear, chi2,  smearErr,  chi2Err);
	
	StyleGraph(gGamma, kRed, "#Gamma(#eta#rightarrow#gamma#gamma) [keV]");
	StyleGraph(gPhi,   kBlue, "Interference Angle [#circ]");
	StyleGraph(gCoh,   kBlack, "N.C. Normalization Factor");
	StyleGraph(gInc,   kGreen, "N.I. Normalization Factor");
	StyleGraph(gChi2,  kCyan,  "#chi^{2} / n.d.f.");
	
	TCanvas *cResults = new TCanvas("cResults","Results",950,600);
	StyleCanvas(cResults);
	gGamma->Draw("AP");
	//gPhi->Draw("P same");
	//gCoh->Draw("P same");
	//gInc->Draw("P same");
	//gChi2->Draw("P same");
	
	return;
}

void StyleGraph(TGraphErrors *g, int color, TString axisName) 
{
	g->SetTitle("");
	
	g->GetXaxis()->SetTitle("Additional Smearing on Reconstructed Photon Angle [%]");
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
