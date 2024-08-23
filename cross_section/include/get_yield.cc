#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/style.cc"

void fitAngularYield()
{
	if(DRAW_MGG_FITS) {
		c_fits = new TCanvas("c_fits", "Fit Results", 1200, 800);
		c_fits->Divide(4,3);
		for(int i=0; i<12; i++) {
			p_fits[i] = (TPad*)c_fits->cd(i+1);
			styleCanvas(p_fits[i]);
		}
	}
	
	TF1 *loc_f_bkgd;
	int n_pars = initFitFunction(&loc_f_bkgd, "bkgd");
	
	int counter = 0;
	for(int itbin=0; itbin<m_n_bins_fit; itbin++) {
		
		int loc_bin_lo = m_rebins_theta*(itbin);
		int loc_bin_hi = m_rebins_theta*(itbin+1);
		double loc_min_angle = 0.01*(double)(loc_bin_lo);
		double loc_max_angle = 0.01*(double)(loc_bin_hi);
		double loc_angle     = 0.5*(loc_min_angle+loc_max_angle);
		
		TH1F *loc_h1_full = (TH1F*)h_mgg_vs_theta_full->ProjectionY(
			Form("h1_full_%03d",itbin), loc_bin_lo+1, loc_bin_hi, "e");
		loc_h1_full->Rebin(m_rebins_mgg);
		styleMggHistogram(loc_h1_full);
		
		fit_mgg(loc_h1_full, f_mgg_fit_functions[itbin].second, f_mgg_fit_functions[itbin].first);
		
		if(DRAW_MGG_FITS) {
			p_fits[counter%12]->cd();
			loc_h1_full->Draw("PE1");
			f_mgg_fit_functions[itbin].second->Draw("same");
			if((counter+1)%12==0) c_fits->Update();
		}
		counter++;
		
		// get yield from integrated counts, minus the background:
		
		loc_f_bkgd->SetParameters(f_mgg_fit_functions[itbin].second->GetParameters());
		switch(m_signal_fit_option) {
			case 1:
				loc_f_bkgd->SetParameter("N_{#eta}", 0.0);
				break;
			case 2:
				loc_f_bkgd->SetParameter("N_{#eta,1}", 0.0);
				loc_f_bkgd->SetParameter("N_{#eta,2}", 0.0);
				break;
			case 3:
				loc_f_bkgd->SetParameter("N_{#eta}", 0.0);
				break;
			case 4:
				loc_f_bkgd->SetParameter("N_{#eta,1}", 0.0);
				loc_f_bkgd->SetParameter("N_{#eta,2}", 0.0);
				break;
		}
		
		int min_mgg_bin = loc_h1_full->FindBin(m_min_mgg_cut);
		int max_mgg_bin = loc_h1_full->FindBin(m_max_mgg_cut);
		
		double min_mgg_cut = m_min_mgg_cut - 0.5*(loc_h1_full->GetBinCenter(2)-loc_h1_full->GetBinCenter(1));
		double max_mgg_cut = m_max_mgg_cut + 0.5*(loc_h1_full->GetBinCenter(2)-loc_h1_full->GetBinCenter(1));
		
		double n_signal = loc_h1_full->Integral(min_mgg_bin, max_mgg_bin);
		double n_bkgd   = loc_f_bkgd->Integral(min_mgg_cut, max_mgg_cut) / m_mgg_bin_size;
		
		m_angular_yield[itbin].first  = n_signal - n_bkgd;
		m_angular_yield[itbin].second = sqrt(n_signal + n_bkgd);
	}
	
	return;
}

void plotAngularYield()
{
	double *angle     = new double[m_n_bins_fit];
	double *angle_err = new double[m_n_bins_fit];
	double *yield     = new double[m_n_bins_fit];
	double *yield_err = new double[m_n_bins_fit];
	
	double maximum_yield = 0.;
	
	for(int ibin=0; ibin<m_n_bins_fit; ibin++) {
		angle[ibin]     = m_angular_bin[ibin].first;
		angle_err[ibin] = m_angular_bin[ibin].second;
		yield[ibin]     = m_angular_yield[ibin].first;
		yield_err[ibin] = m_angular_yield[ibin].second;
		if(m_angular_yield[ibin].first > maximum_yield) {
			maximum_yield = m_angular_yield[ibin].first;
		}
	}
	
	gYield = new TGraphErrors(m_n_bins_fit, angle, yield, angle_err, yield_err);
	gYield->SetName("eta_gg_yield");
	gYield->GetXaxis()->SetTitle("Polar Angle, #theta_{#gamma#gamma} [#circ]");
	gYield->GetXaxis()->SetTitleSize(0.05);
	gYield->GetXaxis()->SetTitleOffset(1.0);
	gYield->GetXaxis()->CenterTitle("");
	gYield->GetYaxis()->SetTitle(Form("N#left(#eta#rightarrow#gamma#gamma#right) [counts / %.02f#circ]", m_theta_bin_size));
	gYield->GetYaxis()->SetTitleSize(0.05);
	gYield->GetYaxis()->SetTitleOffset(1.0);
	gYield->GetYaxis()->CenterTitle("");
	gYield->SetTitle("");
	gYield->SetMarkerStyle(4);
	gYield->SetMarkerSize(0.7);
	gYield->SetMarkerColor(kBlue+2);
	gYield->SetLineColor(kBlue+2);
	gYield->SetLineWidth(1);
	
	gYield->GetYaxis()->SetRangeUser(0.0, 1.1*maximum_yield);
	
	cYield = new TCanvas("cYield", "Angular Yield", 700, 500);
	styleCanvas(cYield);
	
	cYield->cd();
	gYield->Draw("APE");
	
	return;
}

int saveAngularYield(TString output_fname)
{
	if(!gSystem->AccessPathName(output_fname.Data())) return 1;
	
	TFile *fOutput = new TFile(output_fname.Data(), "RECREATE");
	fOutput->cd();
	gYield->Write();
	fOutput->Write();
	fOutput->Close();
	
	return 0;
}
