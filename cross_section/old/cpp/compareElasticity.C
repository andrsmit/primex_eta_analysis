
void compareElasticity()
{
	
	TString primexDir = "/work/halld/home/andrsmit/primex_eta_analysis";
	
	TFile *fCoh = new TFile(Form("%s/eta_gg_matrix/analyze_trees/rootFiles/phase3/phase3.root", primexDir.Data()) ,"READ");
	TFile *fInc = new TFile(Form("%s/bggen_ana/analyze_trees/rootFiles/phase3/Helium.root", primexDir.Data()) ,"READ");
	
	TH2F *h2Coh = (TH2F*)fCoh->Get("VetoOption6/elasticity_veto_6")->Clone("h2Coh");
	TH2F *h2Inc = (TH2F*)fInc->Get("elas_bggen_signal")->Clone("h2Inc");
	TH2F *h2Bkg = (TH2F*)fInc->Get("elas_bggen_etapion")->Clone("h2Bkg");
	
	TFile *fData  = new TFile(Form("%s/analyze_trees/rootFiles/phase3/full_target_bfield.root", primexDir.Data()), "READ");
	TFile *fEmpty = new TFile(Form("%s/analyze_trees/rootFiles/phase3/empty_target_bfield.root", primexDir.Data()), "READ");
	
	TH2F *h2Data  = (TH2F*)fData->Get("VetoOption6/elasticity_veto_6")->Clone("h2Data");
	TH2F *h2Empty = (TH2F*)fEmpty->Get("VetoOption6/elasticity_veto_6")->Clone("h2Empty");
	h2Empty->Scale(2.4);
	h2Data->Add(h2Empty,-1.0);
	
	TCanvas *cMgg = new TCanvas("cElas", "Elasticity", 1600, 900);
	cMgg->Divide(4,3);
	
	for(int i=0; i<12; i++) {
		
		TPad *locPad = (TPad*)cMgg->cd(i+1);
		
		int minBin = i*35 + 1;
		int maxBin = (i+1)*35;
		
		TH1F *h1Coh = (TH1F*)h2Coh->ProjectionY(Form("h1Coh_%d",i), minBin, maxBin);
		TH1F *h1Inc = (TH1F*)h2Inc->ProjectionY(Form("h1Inc_%d",i), minBin, maxBin);
		TH1F *h1Bkg = (TH1F*)h2Bkg->ProjectionY(Form("h1Bkg_%d",i), minBin, maxBin);
		
		TH1F *h1Data = (TH1F*)h2Data->ProjectionY(Form("h1Data_%d",i), minBin, maxBin);
		
		h1Coh->SetLineColor(kBlue);
		h1Inc->SetLineColor(kRed);
		h1Bkg->SetLineColor(kGreen);
		h1Data->SetLineColor(kBlack);
		h1Inc->Rebin(2);
		h1Coh->Rebin(2);
		h1Bkg->Rebin(2);
		h1Data->Rebin(2);
		
		double dataIntegral = h1Data->Integral(h1Data->FindBin(0.9), h1Data->FindBin(1.1));
		
		h1Data->Scale(0.9*h1Inc->Integral()/dataIntegral);
		h1Coh->Scale(h1Inc->Integral()/h1Coh->Integral());
		
		double locBinWidth = h2Inc->GetXaxis()->GetBinWidth(1);
		double locMinAngle = h2Inc->GetXaxis()->GetBinCenter(minBin) - 0.5*locBinWidth;
		double locMaxAngle = h2Inc->GetXaxis()->GetBinCenter(maxBin) + 0.5*locBinWidth;
		
		h1Data->SetTitle(Form("%.2f#circ < #theta_{rec} < %.2f#circ", locMinAngle, locMaxAngle));
		h1Data->GetXaxis()->SetRangeUser(0.80, 1.15);
		h1Data->GetXaxis()->SetTitleSize(0.05);
		h1Data->GetXaxis()->SetTitleOffset(1.0);
		h1Data->GetXaxis()->CenterTitle(true);
		
		h1Inc->SetTitle(Form("%.2f#circ < #theta_{rec} < %.2f#circ", locMinAngle, locMaxAngle));
		h1Inc->GetXaxis()->SetRangeUser(0.80, 1.15);
		h1Inc->GetXaxis()->SetTitleSize(0.05);
		h1Inc->GetXaxis()->SetTitleOffset(1.0);
		h1Inc->GetXaxis()->CenterTitle(true);
		
		h1Data->Draw("PE");
		h1Inc->Draw("PE same");
		h1Coh->Draw("hist same");
		h1Bkg->Draw("PE same");
		
		locPad->Update();
		
		TLine *l1 = new TLine(1.0, gPad->GetUymin(), 1.0, gPad->GetUymax());
		l1->SetLineStyle(4);
		l1->SetLineWidth(2);
		l1->SetLineColor(2);
		l1->Draw("same");
	}
	
	return;
}
