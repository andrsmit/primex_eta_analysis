#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/mgg_fit_function.cc"

int getDataHistograms();
void initParameters();
int initFitFunction(TF1 **f1, TString func_name="f_mgg_fit");

int getDataHistograms() 
{
	if(gSystem->AccessPathName(m_root_fname) || gSystem->AccessPathName(m_empty_target_root_fname)) {
		return 1;
	}
	
	TFile *loc_fIn;
	
	// Get full-target distribution:
	
	loc_fIn = new TFile(m_root_fname.Data(), "READ");
	h_mgg_vs_theta_full = (TH2F*)loc_fIn->Get(m_hist_name.Data())->Clone("mgg_vs_theta_full");
	h_mgg_vs_theta_full->SetDirectory(0);
	loc_fIn->Close();
	
	// Get empty-target distribution:
	
	loc_fIn = new TFile(m_empty_target_root_fname.Data(), "READ");
	h_mgg_vs_theta_empty = (TH2F*)loc_fIn->Get(m_hist_name.Data())->Clone("mgg_vs_theta_empty");
	h_mgg_vs_theta_empty->Scale(m_empty_target_flux_ratio);
	h_mgg_vs_theta_empty->SetDirectory(0);
	loc_fIn->Close();
	
	return 0;
}

void initParameters() {
	
	m_n_bins_fit = (int)(650/m_rebins_theta);
	
	for(int ibin=0; ibin < m_n_bins_fit; ibin++) {
		
		int loc_bin_lo = m_rebins_theta*(ibin);
		int loc_bin_hi = m_rebins_theta*(ibin+1);
		double loc_min_angle = 0.01*(double)(loc_bin_lo);
		double loc_max_angle = 0.01*(double)(loc_bin_hi);
		double loc_angle     = 0.5*(loc_min_angle+loc_max_angle);
		double loc_angle_err = (loc_max_angle-loc_min_angle)/2.;
		
		m_angular_bin.push_back({loc_angle,loc_angle_err});
		m_angular_yield.push_back({0.0, 0.0});
		m_angular_yield_empty.push_back({0.0, 0.0});
		
		TF1 *loc_f1;
		int n_pars = initFitFunction(&loc_f1, Form("f_mgg_fit_%03d",ibin));
		
		// store fit function in a vector:
		
		f_mgg_fit_functions.push_back({n_pars,loc_f1});
	}
	
	return;
}

int initFitFunction(TF1 **f1, TString func_name) {
	
	// get the number of fit parameters based on the specified configuration:
	
	vector<TString> par_names;
	par_names.clear();
	
	// eta fit parameters:
	int n_eta_fit_parameters = 0;
	switch(m_signal_fit_option) {
		case 1:
			n_eta_fit_parameters = 3;
			par_names.push_back("N_{#eta}");
			par_names.push_back("#mu_{#eta}");
			par_names.push_back("#sigma_{#eta}");
			break;
		case 2:
			n_eta_fit_parameters = 6;
			par_names.push_back("N_{#eta,1}");
			par_names.push_back("N_{#eta,2}");
			par_names.push_back("#mu_{#eta,1}");
			par_names.push_back("#mu_{#eta,2}-#mu_{#eta,1}");
			par_names.push_back("#sigma_{#eta,1}");
			par_names.push_back("#sigma_{#eta,2}");
			break;
		case 3:
			n_eta_fit_parameters = 5;
			par_names.push_back("N_{#eta}");
			par_names.push_back("#mu_{#eta}");
			par_names.push_back("#sigma_{#eta}");
			par_names.push_back("#alpha_{#eta}");
			par_names.push_back("n_{#eta}");
			break;
		case 4:
			n_eta_fit_parameters = 8;
			par_names.push_back("N_{#eta,1}");
			par_names.push_back("#mu_{#eta,1}");
			par_names.push_back("#sigma_{#eta,1}");
			par_names.push_back("#alpha_{#eta,1}");
			par_names.push_back("n_{#eta,1}");
			par_names.push_back("N_{#eta,2}");
			par_names.push_back("#mu_{#eta,2}");
			par_names.push_back("#sigma_{#eta,2}");
			break;
	}
	
	// omega->pi0+gamma fit parameters:
	int n_omega_fit_parameters = 5;
	par_names.push_back("N_{#omega}");
	par_names.push_back("#mu_{#omega}");
	par_names.push_back("#sigma_{#omega}");
	par_names.push_back("#alpha_{#omega}");
	par_names.push_back("n_{#omega}");
	
	// background fit parameters:
	int n_background_fit_parameters = 0;
	switch(m_background_fit_option) {
		case 1:
			n_background_fit_parameters = 4;
			break;
		case 2:
			n_background_fit_parameters = 5;
			break;
	}
	for(int ipar=0; ipar<n_background_fit_parameters; ipar++) par_names.push_back(Form("p%d",ipar));
	
	// eta-primex fit parameters:
	int n_etap_fit_parameters = 3;
	par_names.push_back("N_{#eta'}");
	par_names.push_back("#mu_{#eta'}");
	par_names.push_back("#sigma_{#eta'}");
	
	int n_fit_parameters = n_eta_fit_parameters + n_omega_fit_parameters 
		+ n_background_fit_parameters + n_etap_fit_parameters;
	
	// initialize fit function for each angular bin:
	
	*f1 = new TF1(func_name.Data(), mgg_fit, m_min_bkgd_fit, m_max_bkgd_fit, n_fit_parameters);
	
	// set names for each parameter:
	
	for(int ipar = 0; ipar < par_names.size(); ipar++) (*f1)->SetParName(ipar, par_names[ipar]);
	
	return n_fit_parameters;
}
