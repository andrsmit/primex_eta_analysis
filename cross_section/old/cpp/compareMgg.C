
void compareMgg()
{
	gStyle->SetOptStat(0);
	
	TString primexDir = "/work/halld/home/andrsmit/primex_eta_analysis";
	
	TFile *fCoh = new TFile(Form("%s/eta_gg_matrix/analyze_trees/rootFiles/phase3/phase3.root", primexDir.Data()) ,"READ");
	TFile *fInc = new TFile(Form("%s/bggen_ana/analyze_trees/rootFiles/phase3/Helium_VetoOption4.root", primexDir.Data()) ,"READ");
	
	TH2F *h2Coh = (TH2F*)fCoh->Get("VetoOption4/mgg_const_veto_4")->Clone("h2Coh");
	TH2F *h2Inc = (TH2F*)fInc->Get("mgg_const_bggen_signal")->Clone("h2Inc");
	TH2F *h2Bkg = (TH2F*)fInc->Get("mgg_const_bggen_etapion")->Clone("h2Bkg");
	TH2F *h2Bkg2 = (TH2F*)fInc->Get("mgg_const_bggen_eta2pion")->Clone("h2Bkg2");
	h2Bkg->Add(h2Bkg2);
	
	TH2F *h2CohCorr = (TH2F*)fCoh->Get("VetoOption4/mgg_const_coh_veto_4")->Clone("h2CohCorr");
	
	TFile *fData  = new TFile(Form("%s/analyze_trees/rootFiles/phase3/default/full_target_bfield.root", primexDir.Data()), "READ");
	TFile *fEmpty = new TFile(Form("%s/analyze_trees/rootFiles/phase3/default/empty_target_bfield.root", primexDir.Data()), "READ");
	
	TH2F *h2Data  = (TH2F*)fData->Get("VetoOption4/mgg_const_coh_veto_4")->Clone("h2Data");
	TH2F *h2Empty = (TH2F*)fEmpty->Get("VetoOption4/mgg_const_coh_veto_4")->Clone("h2Empty");
	h2Empty->Scale(2.4);
	h2Data->Add(h2Empty,-1.0);
	
	TCanvas *cMgg = new TCanvas("cMgg", "cMgg", 2000, 900);
	cMgg->Divide(4,2);
	
	TCanvas *c1 = new TCanvas("c1","c1",950,700);
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	
	for(int i=0; i<8; i++) {
		
		TPad *locPad = (TPad*)cMgg->cd(i+1);
		
		int minBin = i*50 + 1;
		int maxBin = (i+1)*50;
		
		double locBinWidth = h2Inc->GetXaxis()->GetBinWidth(1);
		double locMinAngle = h2Inc->GetXaxis()->GetBinCenter(minBin) - 0.5*locBinWidth;
		double locMaxAngle = h2Inc->GetXaxis()->GetBinCenter(maxBin) + 0.5*locBinWidth;
		
		TH1F *h1Coh = (TH1F*)h2Coh->ProjectionY(Form("h1Coh_%d",i), minBin, maxBin);
		TH1F *h1Inc = (TH1F*)h2Inc->ProjectionY(Form("h1Inc_%d",i), minBin, maxBin);
		TH1F *h1Bkg = (TH1F*)h2Bkg->ProjectionY(Form("h1Bkg_%d",i), minBin, maxBin);
		
		TH1F *h1Data = (TH1F*)h2Data->ProjectionY(Form("h1Data_%d",i), minBin, maxBin);
		
		TH1F *h1CohCorr = (TH1F*)h2CohCorr->ProjectionY(Form("h1CohCorr_%d",i), minBin, maxBin);
		
		h1Coh->SetLineColor(kBlue);
		h1Coh->SetMarkerColor(kBlue);
		
		h1CohCorr->SetLineColor(kRed);
		h1CohCorr->SetMarkerColor(kRed);
		
		/*
		h1Inc->SetLineColor(kRed);
		h1Inc->SetMarkerColor(kRed);
		h1Inc->SetMarkerStyle(24);
		h1Inc->SetMarkerSize(0.8);
		
		h1Bkg->SetLineColor(kMagenta);
		h1Bkg->SetMarkerColor(kMagenta);
		h1Bkg->SetMarkerStyle(25);
		h1Bkg->SetMarkerSize(0.8);
		
		h1Data->SetLineColor(kBlack);
		h1Data->SetMarkerColor(kBlack);
		h1Data->SetMarkerStyle(20);
		h1Data->SetMarkerSize(0.8);
		
		h1Inc->Rebin(2);
		h1Coh->Rebin(2);
		h1Bkg->Rebin(2);
		h1Data->Rebin(2);
		
		double dataIntegral = h1Data->Integral(h1Data->FindBin(0.5), h1Data->FindBin(0.6));
		double mcIntegral   = h1Inc->Integral(h1Inc->FindBin(0.5), h1Inc->FindBin(0.6)) 
			+ h1Bkg->Integral(h1Bkg->FindBin(0.5), h1Bkg->FindBin(0.6));
		
		h1Inc->Scale(dataIntegral/mcIntegral);
		h1Bkg->Scale(dataIntegral/mcIntegral);
		
		double incIntegral = h1Inc->Integral(h1Inc->FindBin(0.5), h1Inc->FindBin(0.6));
		double cohIntegral = h1Coh->Integral(h1Coh->FindBin(0.5), h1Coh->FindBin(0.6));
		
		h1Coh->Scale(incIntegral/cohIntegral);
		
		h1Data->SetTitle(Form("%.2f#circ < #theta < %.2f#circ", locMinAngle, locMaxAngle));
		h1Data->GetXaxis()->SetRangeUser(0.45, 0.65);
		h1Data->GetXaxis()->SetTitleSize(0.05);
		h1Data->GetXaxis()->SetTitleOffset(1.0);
		h1Data->GetXaxis()->CenterTitle(true);
		
		h1Inc->SetTitle(Form("%.2f#circ < #theta < %.2f#circ", locMinAngle, locMaxAngle));
		h1Inc->GetXaxis()->SetRangeUser(0.45, 0.65);
		h1Inc->GetXaxis()->SetTitleSize(0.05);
		h1Inc->GetXaxis()->SetTitleOffset(1.0);
		h1Inc->GetXaxis()->CenterTitle(true);
		
		h1Data->Draw("PE");
		h1Inc->Draw("PE same");
		h1Coh->Draw("hist same");
		h1Bkg->Draw("PE same");
		
		locPad->Update();
		*/
		
		h1Coh->SetTitle(Form("%.2f#circ < #theta < %.2f#circ", locMinAngle, locMaxAngle));
		h1Coh->GetXaxis()->SetRangeUser(0.45, 0.65);
		h1Coh->GetXaxis()->SetTitleSize(0.05);
		h1Coh->GetXaxis()->SetTitleOffset(1.0);
		h1Coh->GetXaxis()->CenterTitle(true);
		h1Coh->SetMinimum(0.);
		
		h1Coh->Draw("hist");
		h1CohCorr->Draw("same hist");
		
		locPad->Update();
		
		TLine *l1 = new TLine(0.547, gPad->GetUymin(), 0.547, gPad->GetUymax());
		l1->SetLineStyle(4);
		l1->SetLineWidth(2);
		l1->SetLineColor(2);
		l1->Draw("same");
		
		if(i==3) {
			c1->cd();
			h1Data->Draw("PE");
			h1Inc->Draw("PE same");
			h1Coh->Draw("hist same");
			h1Bkg->Draw("PE same");
			
			TLegend *leg = new TLegend(0.08, 0.65, 0.45, 0.87);
			leg->AddEntry(h1Data, "Data", "PE");
			leg->AddEntry(h1Coh, "#gamma+^{4}He#rightarrow#eta+^{4}He (MC)", "l");
			leg->AddEntry(h1Inc, "#gamma+^{4}He#rightarrow#eta+p(n)+X (MC)", "PE");
			leg->AddEntry(h1Bkg, "#gamma+^{4}He#rightarrow#eta+#pi+p(n)+X (MC)", "PE");
			leg->Draw();
		}
	}
	
	return;
}
