void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

double CohFunc(double *x, double *par) {
	double f = par[0]*pow(x[0]-par[1],par[2]);
	return f;
}

void plotResults() {
	
	vector<double> shiftVec;
	vector<pair<double,double>> gammaVec;
	vector<pair<double,double>>   cohVec;
	vector<pair<double,double>>   incVec;
	vector<pair<double,int   >>  chi2Vec;
	
	double dt=-0.25;
	while(dt<0.3) {
		int dt_int = (int)((dt+0.2501)/0.05);
		printf("dt, dt_int: %f %d\n",dt,dt_int);
		
		TString filename = Form("results_%02d.txt", dt_int);
		if(gSystem->AccessPathName(filename.Data())) {
			dt += 0.05;
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
		
		shiftVec.push_back(dt);
		gammaVec.push_back({locGamma, locGammaErr});
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
		
		dt += 0.05;
	}
	
	int nPoints = shiftVec.size();
	
	double *shift    = new double[nPoints];
	double *shiftErr = new double[nPoints];
	double *gamma    = new double[nPoints];
	double *gammaErr = new double[nPoints];
	double *coh      = new double[nPoints];
	double *cohErr   = new double[nPoints];
	double *inc      = new double[nPoints];
	double *incErr   = new double[nPoints];
	double *chi2     = new double[nPoints];
	double *chi2Err  = new double[nPoints];
	
	for(int i=0; i<nPoints; i++) {
		shift[i]    = shiftVec[i];
		shiftErr[i] = 0.0;
		gamma[i]    = gammaVec[i].first;
		gammaErr[i] = gammaVec[i].second;
		coh[i]      = cohVec[i].first;
		cohErr[i]   = cohVec[i].second;
		inc[i]      = incVec[i].first;
		incErr[i]   = incVec[i].second;
		chi2[i]     = chi2Vec[i].first / ((double)chi2Vec[i].second);
		chi2Err[i]  = 0.0;
	}
	
	TGraphErrors *gGamma = new TGraphErrors(nPoints, shift, gamma, shiftErr, gammaErr);
	TGraphErrors *gCoh   = new TGraphErrors(nPoints, shift,   coh, shiftErr,   cohErr);
	TGraphErrors *gInc   = new TGraphErrors(nPoints, shift,   inc, shiftErr,   incErr);
	TGraphErrors *gChi2  = new TGraphErrors(nPoints, shift,  chi2, shiftErr,  chi2Err);
	
	StyleGraph(gGamma, kBlue, "#Gamma(#eta#rightarrow#gamma#gamma) [keV]");
	StyleGraph(gCoh,   kBlue, "N.C. Normalization Factor");
	StyleGraph(gInc,   kBlue, "N.I. Normalization Factor");
	StyleGraph(gChi2,  kRed,  "#chi^{2} / n.d.f.");
	
	Double_t xmin, xmax, dx;
	Double_t ymin, ymax, dy;
	
	//----------------------------------------------------------//
	
	TCanvas *cGamma = new TCanvas("cGamma","Gamma",950,600);
	cGamma->cd();
	gGamma->GetXaxis()->SetTitle("#Delta#theta_{1,2} [#circ]");
	gGamma->Draw("AP");
	/*
	TF1 *fCoh = new TF1("fCoh", CohFunc, 7.0, 11.0, 3);
	fCoh->SetParameters(0.2, 7.0, 0.2);
	fCoh->SetParLimits(0, 0.0, 1.0);
	fCoh->SetParLimits(1, 0.0, 10.0);
	fCoh->SetParLimits(2, 0.01, 0.9);
	fCoh->FixParameter(1, 7.0);
	gCoh->Fit("fCoh", "R");
	fCoh->SetLineColor(kRed);
	fCoh->SetLineStyle(2);
	fCoh->SetLineWidth(2);
	fCoh->Draw("same");*/
	
	/*TCanvas *cInc = new TCanvas("cInc","Inc",950,600);
	cInc->cd();
	gInc->GetXaxis()->SetTitle("N.I. Normarlization Factor (fixed)");
	gInc->Draw("AP");*/
	
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
