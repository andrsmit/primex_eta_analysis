void makeTable()
{
	/* table of angular yield vs. bin for:
	   1. Full target counts
	   2. Empty target counts
	   3. Full - Empty counts
	   4. Single-eta yield
	   5. Hadronic bkgd yield
	   6. smooth background yield
	   7. omega yield
	*/
	
	TFile *fIn = new TFile("yield_phase3_VetoOption6.root", "READ");
	TH1F *hFullCounts  = (TH1F*)fIn->Get("Counts");
	TH1F *hEmptyCounts = (TH1F*)fIn->Get("EmptyCounts");
	TH1F *hSubCounts = (TH1F*)hFullCounts->Clone("SubCounts");
	hSubCounts->Add(hEmptyCounts,-1.0);
	
	TH1F *hYield = (TH1F*)fIn->Get("AngularYieldFit");
	TH1F *hHadBkgd = (TH1F*)fIn->Get("HadronicBkgdYield");
	TH1F *hEtaPi = (TH1F*)fIn->Get("EtaPionYield");
	hHadBkgd->Add(hEtaPi);
	TH1F *hOmega = (TH1F*)fIn->Get("OmegaYield");
	TH1F *hBkgd  = (TH1F*)fIn->Get("BkgdYield");
	
	ofstream latexTableStr("latex_table.txt");
	latexTableStr << "		$\\theta_{rec}$ [$^{\\circ}$] & \\specialcell{Total Counts\\\\(Full Target)} & \\specialcell{Total Counts\\\\(Empty Target)} & \\specialcell{Total Counts\\\\(Full-Empty)} & \\specialcell{Hadr.\\\\Bkgd} & \\specialcell{$\\omega$\\\\Bkgd} & \\specialcell{Additional\\\\Bkgd} & \\specialcell{Final\\\\Yield} \\\\\n";
	latexTableStr << "		\\hline\n";
	for(int ibin=1; ibin<=hYield->GetXaxis()->GetNbins(); ibin++) {
		char buf[256];
		double locAngle = hYield->GetXaxis()->GetBinCenter(ibin);
		double locAngleErr = hYield->GetXaxis()->GetBinWidth(ibin);
		double minAngle = locAngle - 0.5*locAngleErr;
		double maxAngle = locAngle + 0.5*locAngleErr;
		
		sprintf(buf, "		%.2f-%.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f  &  %.2f \\\\\n", 
			minAngle, maxAngle, 
			hFullCounts->GetBinContent(ibin),
			hEmptyCounts->GetBinContent(ibin),
			hSubCounts->GetBinContent(ibin),
			hHadBkgd->GetBinContent(ibin),
			hOmega->GetBinContent(ibin),
			hBkgd->GetBinContent(ibin),
			hYield->GetBinContent(ibin)
		);
		latexTableStr << buf;
	}
	latexTableStr << "		\\hline\n";
	latexTableStr.close();
	
	
	return;
}
