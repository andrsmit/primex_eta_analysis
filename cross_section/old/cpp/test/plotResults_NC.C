void StyleGraph(TGraphErrors *g, int color=kBlack, TString axisName="");
void StyleCanvas(TCanvas *c);

double CohFunc(double *x, double *par) {
	double f = par[0]*pow(x[0]-par[1],par[2]);
	return f;
}

void plotResults_NC() {
	
	vector<double> egammaVec;
	vector<pair<double,double>> cohVec;
	vector<pair<double,double>> incVec;
	vector<pair<double,int>> chi2Vec;
	
	double eb=7.0;
	while(eb<12.0) {
		TString filename = Form("fit_0.0_3.5_%.2fGeV.txt", eb);
		if(gSystem->AccessPathName(filename.Data())) {
			eb += 0.05;
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
		
		egammaVec.push_back(eb);
		cohVec.push_back({locCoh, locCohErr});
		incVec.push_back({locInc, locIncErr});
		chi2Vec.push_back({locChi2, locNDF});
		
		eb += 0.05;
	}
	
	int nPoints = egammaVec.size();
	
	double *egamma    = new double[nPoints];
	double *egammaErr = new double[nPoints];
	double *coh       = new double[nPoints];
	double *cohErr    = new double[nPoints];
	double *inc       = new double[nPoints];
	double *incErr    = new double[nPoints];
	double *chi2      = new double[nPoints];
	double *chi2Err   = new double[nPoints];
	
	for(int i=0; i<nPoints; i++) {
		egamma[i]    = egammaVec[i];
		egammaErr[i] = 0.0;
		coh[i]       = cohVec[i].first;
		cohErr[i]    = cohVec[i].second;
		inc[i]       = incVec[i].first;
		incErr[i]    = incVec[i].second;
		chi2[i]      = chi2Vec[i].first / ((double)chi2Vec[i].second);
		chi2Err[i]   = 0.0;
	}
	
	TGraphErrors *gCoh  = new TGraphErrors(nPoints, egamma,  coh, egammaErr,  cohErr);
	TGraphErrors *gInc  = new TGraphErrors(nPoints, egamma,  inc, egammaErr,  incErr);
	TGraphErrors *gChi2 = new TGraphErrors(nPoints, egamma, chi2, egammaErr, chi2Err);
	
	StyleGraph(gCoh,  kBlue, "N.C. Normalization Factor");
	StyleGraph(gInc,  kBlue, "N.I. Normalization Factor");
	StyleGraph(gChi2, kRed,  "#chi^{2} / n.d.f.");
	
	Double_t xmin, xmax, dx;
	Double_t ymin, ymax, dy;
	
	//----------------------------------------------------------//
	
	TCanvas *cCoh = new TCanvas("cCoh","Coh",950,600);
	cCoh->cd();
	gCoh->GetXaxis()->SetTitle("N.C. Normarlization Factor (fixed)");
	gCoh->Draw("AP");
	
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
	fCoh->Draw("same");
	
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
