void StyleHist(TH1F*,int);
void StyleCanvas(TCanvas*);

void compareTOF_counts()
{
	gStyle->SetOptStat(0);
	
	TFile *f_noTOF = new TFile("yield_phase3_VetoOption6_noTOF.root", "READ");
	TH1F *h_noTOF = (TH1F*)f_noTOF->Get("Counts")->Clone("cs_noTOF");
	TH1F *h_noTOF_e = (TH1F*)f_noTOF->Get("EmptyCounts")->Clone("cs_noTOF_e");
	//h_noTOF->Add(h_noTOF_e,-1.0);
	h_noTOF->SetDirectory(0);
	f_noTOF->Close();
	
	TFile *f_TOF = new TFile("yield_phase3_VetoOption6_TOF.root", "READ");
	TH1F *h_TOF = (TH1F*)f_TOF->Get("Counts")->Clone("cs_TOF");
	TH1F *h_TOF_e = (TH1F*)f_TOF->Get("EmptyCounts")->Clone("cs_TOF_e");
	//h_TOF->Add(h_TOF_e,-1.0);
	h_TOF->SetDirectory(0);
	f_TOF->Close();
	
	TFile *f_singleTOF = new TFile("yield_phase3_VetoOption6_singleTOF.root", "READ");
	TH1F *h_singleTOF = (TH1F*)f_singleTOF->Get("Counts")->Clone("cs_singleTOF");
	TH1F *h_singleTOF_e = (TH1F*)f_singleTOF->Get("EmptyCounts")->Clone("cs_singleTOF_e");
	//h_noTOF->Add(h_noTOF_e,-1.0);
	h_singleTOF->SetDirectory(0);
	f_singleTOF->Close();
	
	StyleHist(h_noTOF, kBlack);
	StyleHist(h_TOF, kBlue);
	StyleHist(h_singleTOF, kRed);
	
	TCanvas *c = new TCanvas("c","c",950,700);
	StyleCanvas(c);
	
	h_noTOF->Add(h_singleTOF,-1.0);
	
	h_noTOF->Draw("PE1X0");
	//h_TOF->Draw("PE1X0 same");
	//h_singleTOF->Draw("PE1X0 same");
	
	TLegend *leg1 = new TLegend(0.15, 0.7, 0.45, 0.9);
	leg1->AddEntry(h_noTOF,         "No TOF Veto", "PE1X0");
	leg1->AddEntry(h_TOF,        "Loose TOF Veto", "PE1X0");
	leg1->AddEntry(h_singleTOF, "Strict TOF Veto", "PE1X0");
	leg1->Draw();
	
	TH1F *hr_TOF       = (TH1F*)h_TOF->Clone("hr_TOF");
	TH1F *hr_singleTOF = (TH1F*)h_singleTOF->Clone("hr_singleTOF");
	
	hr_TOF->Divide(h_noTOF);
	hr_singleTOF->Divide(h_noTOF);
	hr_TOF->GetYaxis()->SetTitle("cross section ratio");
	
	TCanvas *cr = new TCanvas("cr","cr",950,700);
	StyleCanvas(cr);
	
	hr_TOF->Draw("PE1X0");
	hr_singleTOF->Draw("PE1X0 same");
	
	cr->Update();
	TLine *l1 = new TLine(gPad->GetUxmin(), 1.0, gPad->GetUxmax(), 1.0);
	l1->SetLineStyle(2); 
	l1->SetLineWidth(2);
	l1->SetLineColor(kBlack);
	l1->Draw("same");
	
	TLegend *leg2 = new TLegend(0.59, 0.7, 0.89, 0.9);
	leg2->AddEntry(l1, "No TOF Veto", "l");
	leg2->AddEntry(hr_TOF,        "Loose TOF Veto", "PE1X0");
	leg2->AddEntry(hr_singleTOF, "Strict TOF Veto", "PE1X0");
	leg2->Draw();
	
	return;
}

void StyleHist(TH1F *h1, int lineColor=kBlack)
{
	h1->SetLineColor(lineColor);
	h1->SetMarkerColor(lineColor);
	
	h1->SetLineWidth(2);
	
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetXaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->GetYaxis()->SetMaxDigits(2);
}

void StyleCanvas(TCanvas *c1)
{
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	c1->SetTickx(); c1->SetTicky();
	c1->cd();
}
