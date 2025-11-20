void compare_acc(int vetoOption=6)
{
	
	TFile *fQF  = new TFile(Form("acceptance_qf_veto%d.root", vetoOption), "READ");
	
	TH1F *hQF = (TH1F*)fQF->Get("h_Acceptance")->Clone("hQF");
	hQF->SetDirectory(0);
	fQF->Close();
	
	TFile *fCoh = new TFile(Form("acceptance_coh_veto%d.root", vetoOption), "READ");
	
	TH1F *hCoh = (TH1F*)fCoh->Get("h_Acceptance")->Clone("hCoh");
	hCoh->SetDirectory(0);
	fCoh->Close();
	
	hCoh->SetLineColor(kBlue);
	hCoh->SetMarkerColor(kBlue);
	
	TCanvas *cAcc = new TCanvas("cAcc", "cAcc", 1000, 600);
	cAcc->SetTickx(); cAcc->SetTicky();
	cAcc->SetLeftMargin(0.13); cAcc->SetRightMargin(0.07);
	cAcc->SetBottomMargin(0.13); cAcc->SetTopMargin(0.07);
	cAcc->cd();
	hCoh->Draw("hist");
	hQF->Draw("same hist");
	
	
	TH1F *hRatio = (TH1F*)hQF->Clone("hRatio");
	hRatio->Divide(hCoh);
	hRatio->SetLineColor(kBlack);
	hRatio->SetMarkerColor(kBlack);
	hRatio->GetYaxis()->SetRangeUser(0.85,1.15);
	
	TCanvas *cRatio = new TCanvas("cRatio", "cRatio", 1000, 600);
	cRatio->SetTickx(); cRatio->SetTicky();
	cRatio->SetLeftMargin(0.13); cRatio->SetRightMargin(0.07);
	cRatio->SetBottomMargin(0.13); cRatio->SetTopMargin(0.07);
	cRatio->cd();
	hRatio->Draw("hist");
	
	return;
}
