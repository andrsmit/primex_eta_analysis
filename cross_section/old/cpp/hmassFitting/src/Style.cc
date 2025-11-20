#include "CrossSection.h"

void styleMggHistogram(TH1F *h1, int lineColor, int markerStyle) {
	
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("M_{#gamma#gamma}^{Constr.} [GeV/c^{2}]");
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitleOffset(1.0);
	h1->GetXaxis()->CenterTitle(true);
	
	double locBinSize = h1->GetXaxis()->GetBinWidth(1);
	h1->GetYaxis()->SetTitle(Form("counts / %d MeV/c^{2}", (int)(locBinSize*1.e3)));
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleOffset(1.0);
	h1->GetYaxis()->CenterTitle(true);
	
	h1->SetMarkerStyle(markerStyle);
	h1->SetMarkerColor(lineColor);
	h1->SetLineColor(lineColor);
	h1->SetFillColor(lineColor-10);
	h1->SetLineWidth(1);
	
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
