#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/eta_inputs.cc"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/style.cc"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/fit_mgg.cc"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/get_yield.cc"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/plot_fit_results.cc"

char loc_pathName[256] = "/work/halld/home/andrsmit/primex_eta_analysis";

void CrossSection(TString pluginname="primex_eta_analysis/eta_gg", TString hname="mgg_const")
{
	//------------------------------------------------------//
	// Style:
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetStatFont(12);
	
	gStyle->SetStatX(0.89);
	gStyle->SetStatY(0.90);
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.16);
	
	//------------------------------------------------------//
	// Which set of cuts/histogram to use:
	
	TString save_name;
	if(pluginname=="primex_eta_analysis/eta_gg") {
		save_name = hname;
	} else if(pluginname=="primex_eta_analysis_TOF") {
		save_name = "TOF_"+hname;
	} else if(pluginname=="primex_eta_analysis_BCAL") {
		save_name = "BCAL_"+hname;
	} else if(pluginname=="primex_eta_analysis_BEAM") {
		save_name = "BEAM_"+hname;
	} else if(pluginname=="primex_eta_analysis_FCAL") {
		save_name = "FCAL_"+hname;
	} else {
		cout << "invalid plugin name specified." << endl;
		return;
	}
	
	m_hist_name = pluginname+"/"+hname;
	
	//------------------------------------------------------//
	// Which data set to use:
	
	m_phase                   =       3;
	m_luminosity              = 5.88005;
	m_empty_target_flux_ratio = 0.0;
	
	/*
	m_phase                   =       1;
	m_luminosity              = 5.88005;
	m_empty_target_flux_ratio = 6.46002;
	
	m_phase                   =       2;
	m_luminosity              = 1.60250;
	m_empty_target_flux_ratio = 1.94506;
	
	m_phase                   =       3;
	m_luminosity              =    18.0;
	m_empty_target_flux_ratio = 3.08500;
	*/
	
	/*-----------------------------------------------------
	Luminosity Corrections:
		Correct for the cold gas subtracted when empty target is subtracted:
			- Assuming target pressure of 33.8 psia and temperature of 50.4 K
				(measured during empty target runs from phase 3), we get a density of the cold gas of:
				0.00223 g/cm3, which is 1.83% that of the filled target.
				--> Corrected target density should then become: 
					(1.0 - 0.0183) * 0.1217 g/cm3 = 0.1195 g/cm3
				Alternatively, correct the luminosity by multiplying by a factor of 0.9817:
	-------------------------------------------------------*/
	
	m_luminosity *= 0.9817;
	
	//-----------------------------------------------------//
	// Configure invariant mass fitting options:
	
	m_signal_fit_option = 2;
	/*
	Different options for fitting the eta->2gamma signal:
		1. single Gaussian function to describe eta signal
		2. double Gaussian function
		3. Crystal Ball function
		4. Crystal Ball function + Gaussian
	*/
	
	m_background_fit_option = 2;
	/*
	Different options for fitting the eta->2gamma background:
		1. 3rd-order polynomial
		2. exponential function
	*/
	
	m_omega_fit_option = 1;
	/*
	Different options for fitting the omega->pi0+gamma background:
		1. Crystal Ball function
		2. Crystal Ball + Gaussian function
	*/
	
	m_min_bkgd_fit = 0.44;
	m_max_bkgd_fit = 1.05;
	
	initParameters();
	
	//------------------------------------------------------------//
	
	m_root_fname              = Form("%s/data/rootFiles/phase%d/full.root",  loc_pathName, m_phase);
	m_empty_target_root_fname = Form("%s/data/rootFiles/phase%d/empty.root", loc_pathName, m_phase);
	
	if(getDataHistograms()) {
		cout << "Problem accessing histograms from the specified input file. Check path/histogram name." << endl;
		exit(1);
	}
	
	//------------------------------------------------------------//
	
	// do the fits in each angular bin:
	
	DRAW_MGG_FITS = false;
	
	fitAngularYield();
	plotAngularYield();
	
	plotFitResults();
	
	if(saveAngularYield("output2.root")) {
		cout << "Problem saving yield to output ROOT file (filename might already exist)." << endl;
	}
	
	/*	
	TF1 *f_fit;
	int n_pars = initFitFunction(&f_fit, "test_fit");
	
	int min_theta_bin = h_mgg_vs_theta_full->GetXaxis()->FindBin(4.05);
	int max_theta_bin = h_mgg_vs_theta_full->GetXaxis()->FindBin(4.10);
	
	TH1F *h1_full  = (TH1F*)h_mgg_vs_theta_full->ProjectionY( "h1_full",  min_theta_bin, max_theta_bin);
	TH1F *h1_empty = (TH1F*)h_mgg_vs_theta_empty->ProjectionY("h1_empty", min_theta_bin, max_theta_bin);
	
	h1_full->Rebin(m_rebins_mgg);
	h1_empty->Rebin(m_rebins_mgg);
	
	TCanvas *loc_canvas = new TCanvas("loc_canvas", "test", 700, 500);
	
	TH1F *h1_sub = (TH1F*)h1_full->Clone("h1_sub");
	h1_sub->Add(h1_empty, -1.0);
	
	styleMggHistogram(h1_full);
	styleMggHistogram(h1_empty, kBlue);
	styleMggHistogram(h1_sub, kRed, 3);
	styleCanvas(loc_canvas);
	
	loc_canvas->cd();
	h1_full->Draw("PE1");
	h1_empty->Draw("hist same");
	h1_full->Draw("PE1 same");
	h1_sub->Draw("PE1");
	
	fit_mgg(h1_sub, f_fit, n_pars);
	
	f_fit->SetLineColor(kGreen+2);
	f_fit->Draw("same");
	
	TF1 *f_bkgd;
	n_pars = initFitFunction(&f_bkgd, "bkgd");
	f_bkgd->SetParameters(f_fit->GetParameters());
	f_bkgd->SetParameter("N_{#eta}", 0.0);
	f_bkgd->SetParameter("N_{#omega}", 0.0);
	f_bkgd->SetParameter("N_{#eta'}", 0.0);
	f_bkgd->SetLineColor(kBlue);
	f_bkgd->SetLineStyle(2);
	f_bkgd->Draw("same");
	*/
	
	return;
}
