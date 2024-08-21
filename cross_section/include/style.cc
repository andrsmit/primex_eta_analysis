#ifndef _ETAGG_STYLE_
#define _ETAGG_STYLE_

#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"

void styleMggHistogram(TH1F *h1, int line_color=kBlack, int marker_style=25);
void styleCanvas(TCanvas *c1);
void styleCanvas(TPad *p1);

void styleMggHistogram(TH1F *h1, int line_color, int marker_style) {
	
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("M_{#gamma#gamma}^{Constr.} [GeV/c^{2}]");
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitleOffset(1.0);
	h1->GetXaxis()->CenterTitle(true);
	h1->GetYaxis()->SetTitle(Form("counts / %d MeV/c^{2}", (int)(m_mgg_bin_size*1.e3)));
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->GetXaxis()->SetRangeUser(m_min_bkgd_fit, m_max_bkgd_fit);
	
	h1->SetMarkerStyle(marker_style);
	h1->SetMarkerColor(line_color);
	h1->SetLineColor(line_color);
	h1->SetFillColor(line_color-10);
	h1->SetLineWidth(2);
	
	h1->SetMinimum(0.);
	
	return;
}

void styleCanvas(TCanvas *c1) {
	
	c1->SetTickx();
	c1->SetTicky();
	c1->SetLeftMargin(0.13); c1->SetRightMargin(0.07);
	c1->SetBottomMargin(0.13); c1->SetTopMargin(0.07);
	
	return;
}
void styleCanvas(TPad *p1) {
	
	p1->SetTickx();
	p1->SetTicky();
	p1->SetLeftMargin(0.13); p1->SetRightMargin(0.07);
	p1->SetBottomMargin(0.13); p1->SetTopMargin(0.07);
	
	return;
}

#endif
