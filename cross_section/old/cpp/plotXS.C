void styleHist(TH1F *h1, int color=kBlack, int style=1) {
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitleOffset(1.0);
	h1->GetXaxis()->CenterTitle(true);
	h1->GetYaxis()->SetTitle("d#sigma/d#theta [#mub/rad]");
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	h1->SetTitle("");
	h1->SetLineWidth(2);
	h1->SetLineStyle(style);
	h1->SetLineColor(color);
	h1->GetXaxis()->SetRangeUser(0.0, 6.5);
	h1->GetYaxis()->SetRangeUser(0.0, 2.0);
	h1->SetTitle("");//"#gamma + ^{4}He #rightarrow #eta + X");
}

void plotXS(double egamma1=8.0, double egamma2=11.3)
{
	gStyle->SetOptStat(0);
	
	//TFile *fIn = new TFile("he4-eta-xs-sgevorkyan-v2.root","READ");
	//TFile *fIn = new TFile("/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v1.root","READ");
	
	TFile *fIn = new TFile("/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root", "READ");
	//TFile *fIn = new TFile("/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-sgevorkyan-v2-folded.root", "READ");
	
	
	TFile *fData = new TFile("output/omega0_poly2/yield_phase3_VetoOption6.root", "READ");
	TH1F *hData = (TH1F*)fData->Get("CrossSectionFit");
	
	TH2F *h2_prim_xs = (TH2F*)fIn->Get("xs_prim_vs_egam"); 
	TH2F *h2_coh_xs = (TH2F*)fIn->Get("xs_coh_vs_egam"); 
	//TH2F *h2_inc_xs = (TH2F*)fIn->Get("qf_txs_lab"); 
	TH2F *h2_inc_xs = (TH2F*)fIn->Get("xs_inc_vs_egam"); 
	//TH2F *h2_int_xs = (TH2F*)fIn->Get("xs_int_vs_egam"); 
	
	int ebin1 = h2_prim_xs->GetXaxis()->FindBin(egamma1);
	int ebin2 = h2_prim_xs->GetXaxis()->FindBin(egamma2);
	
	double nbins = (double)(ebin2-ebin1+1);
	
	TH1F *h1_prim = (TH1F*)h2_prim_xs->ProjectionY("prim",ebin1,ebin2);
	TH1F *h1_coh  = (TH1F*)h2_coh_xs->ProjectionY("coh",ebin1,ebin2);
	TH1F *h1_inc  = (TH1F*)h2_inc_xs->ProjectionY("inc",ebin1,ebin2);
	//TH1F *h1_int  = (TH1F*)h2_int_xs->ProjectionY("int",ebin,ebin);
	
	h1_prim->Scale(1.0/nbins);
	h1_coh->Scale(1.0/nbins);
	h1_inc->Scale(1.0/nbins);
	
	/*
	h1_prim->Rebin(4);
	h1_coh->Rebin(4);
	h1_inc->Rebin(4);
	
	h1_prim->Scale(1.0/4.0);
	h1_coh->Scale(1.0/4.0);
	h1_inc->Scale(1.0/4.0);
	*/
	
	double phi_int;
	
	/*
	h1_prim->Scale(0.422818/0.515);
	h1_coh->Scale(0.648731);
	h1_inc->Scale(0.287156);
	phi_int = 28.35;
	*/
	//h1_prim->Scale(0.4/0.515);
	h1_coh->Scale(0.64);
	h1_inc->Scale(0.49);
	phi_int = 55.5;
	
	TH1F *h1_int = (TH1F*)h1_prim->Clone("h1_int");
	for(int ibin=1; ibin<=h1_prim->GetXaxis()->GetNbins(); ibin++) {
		double locPrim = h1_prim->GetBinContent(ibin);
		double locCoh  = h1_coh->GetBinContent(ibin);
		double locInt  = 2.0*sqrt(locPrim*locCoh)*cos(phi_int*TMath::DegToRad());
		h1_int->SetBinContent(ibin, locInt);
	}
	
	h1_prim->SetLineColor(kRed);
	h1_coh->SetLineColor(kBlue);
	h1_int->SetLineColor(kMagenta);
	h1_inc->SetLineColor(kGreen);
	
	//h1_int->Scale(1.6);
	
	TH1F *h1_sum = (TH1F*)h1_prim->Clone("sum");
	h1_sum->Add(h1_coh);
	h1_sum->Add(h1_int);
	h1_sum->Add(h1_inc);
	h1_sum->SetLineColor(kBlack);
	
	TCanvas *c1 = new TCanvas("c1","c1",950,700);
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	c1->cd();
	
	styleHist(h1_sum,  kBlack, 1);
	styleHist(h1_prim, kRed,   9);
	styleHist(h1_coh,  kBlue,  7);
	styleHist(h1_int,  kMagenta, 3);
	styleHist(h1_inc,  kGreen, 5);
	
	h1_inc->Draw("hist");
	h1_sum->Draw("same hist");
	h1_int->Draw("same hist");
	h1_coh->Draw("same hist");
	h1_prim->Draw("same hist");
	
	hData->Draw("same");
	
	TLatex lat;
	lat.SetTextFont(42);
	lat.DrawLatexNDC(0.155,0.873,Form("%.1f GeV < E_{#gamma} = %.1f GeV", egamma1, egamma2));
	//lat.DrawLatexNDC(0.669367,0.59,Form("#Gamma_{#eta#rightarrow#gamma#gamma} = 510 eV"));
	//lat.DrawLatexNDC(0.731603,0.50,Form("#phi = 57.1#circ"));
	
	//TLegend *leg = new TLegend(0.132, 0.701, 0.468, 0.928);
	TLegend *leg = new TLegend(0.629, 0.676, 0.995, 0.932);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(h1_sum,  "Total", "l");
	leg->AddEntry(h1_prim, "Primakoff", "l");
	leg->AddEntry(h1_coh,  "Nuclear Coherent", "l");
	leg->AddEntry(h1_int,  "Interference", "l");
	leg->AddEntry(h1_inc,  "Nuclear Incoherent", "l");
	leg->Draw();
	
	return;
}
