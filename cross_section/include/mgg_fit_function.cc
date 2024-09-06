#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"

Double_t mgg_fit(Double_t *x, Double_t *par) {
	
	double loc_mgg = x[0];
	
	// for excluding sub-regions of the fit:
	
	for(int iexc = 0; iexc < m_exclude_regions.size(); iexc++) {
		if(m_exclude_regions[iexc].first < loc_mgg && loc_mgg < m_exclude_regions[iexc].second) {
			TF1::RejectPoint();
			return 0;
		}
	}
	
	//----------------------------------------------------------------------------------//
	// Signal:
	
	Double_t f_eta = 0.;
	
	int n_signal_parameters = 0;
	switch(m_signal_fit_option) {
		case 1:
		{
			// single Gaussian:
			
			n_signal_parameters = 3;
			double     N_eta = par[0];
			double    mu_eta = par[1];
			double sigma_eta = par[2];
			double     A_eta = N_eta * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma_eta;
			
			double loc_x_eta = (loc_mgg - mu_eta)/sigma_eta;
			f_eta = A_eta * exp(-0.5*pow(loc_x_eta,2.0));
			
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			n_signal_parameters = 6;
			double     N1_eta = par[0];
			double     N2_eta = par[1];
			double    mu1_eta = par[2];
			double    mu2_eta = par[3] + mu1_eta;
			double sigma1_eta = par[4];
			double sigma2_eta = par[5];
			double     A1_eta = N1_eta * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma1_eta;
			double     A2_eta = N2_eta * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma2_eta;
			
			double loc_x1_eta = (loc_mgg - mu1_eta)/sigma1_eta;
			double loc_x2_eta = (loc_mgg - mu2_eta)/sigma2_eta;
			f_eta = A1_eta*exp(-0.5*pow(loc_x1_eta,2.0)) + A2_eta*exp(-0.5*pow(loc_x2_eta,2.0));
			
			break;
		}
		case 3:
		{
			// Crystal Ball function:
			
			n_signal_parameters = 5;
			
			double     N_eta = par[0];
			double    mu_eta = par[1];
			double sigma_eta = par[2];
			double     a_eta = par[3];
			double     n_eta = par[4];
			double     A_eta = N_eta * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma_eta;
			//cout << "sigma = " << sigma_eta << endl;
			//cout << "A_eta = " << A_eta     << endl;
			
			double   Acb_eta = pow(n_eta/fabs(a_eta), n_eta) * exp(-0.5*pow(fabs(a_eta),2.0));
			double   Bcb_eta = (n_eta/fabs(a_eta)) - fabs(a_eta);
			double loc_x_eta = (mu_eta - loc_mgg)/sigma_eta;
			if(loc_x_eta > -a_eta) {
				f_eta = A_eta * exp(-0.5*pow(loc_x_eta,2.0));
			} else {
				f_eta = A_eta * Acb_eta * pow(Bcb_eta - loc_x_eta, -n_eta);
			}
			
			break;
		}
		case 4:
		{
			// Crystal Ball + Gaussian function:
			
			n_signal_parameters = 8;
			
			break;
		}
	}
	
	//----------------------------------------------------------------------------------//
	// Background:
	
	//-------------------------//
	// omega->pi0+gamma:
	
	Double_t f_omega = 0.;
	
	double     N_omega = par[0+n_signal_parameters];
	double    mu_omega = par[1+n_signal_parameters];
	double sigma_omega = par[2+n_signal_parameters];
	double     a_omega = par[3+n_signal_parameters];
	double     n_omega = par[4+n_signal_parameters];
	
	double Acb_omega = pow(n_omega/fabs(a_omega), n_omega) * exp(-0.5*pow(fabs(a_omega),2.0));
	double Bcb_omega = (n_omega/fabs(a_omega)) - fabs(a_omega);
	
	double loc_x_omega = (loc_mgg - mu_omega)/sigma_omega;
	
	if(loc_x_omega > -a_omega) {
		f_omega = N_omega * exp(-0.5*pow(loc_x_omega,2.0));
	} else {
		f_omega = N_omega * Acb_omega * pow(Bcb_omega - loc_x_omega, -n_omega);
	}
	
	//-------------------------//
	// electromagnetic:
	
	Double_t f_background = 0.;
	
	int n_background_parameters = 0;
	switch(m_background_fit_option) {
		case 1:
		{
			// polynomial background:
			n_background_parameters = 4;
			
			double p0 = par[5+n_signal_parameters];
			double p1 = par[6+n_signal_parameters];
			double p2 = par[7+n_signal_parameters];
			double p3 = par[8+n_signal_parameters];
			
			f_background = p0 + p1*loc_mgg + p2*pow(loc_mgg,2.0) + p3*pow(loc_mgg,3.0);
			break;
		}
		case 2:
		{
			// exponential background:
			n_background_parameters = 5;
			
			double p0 = par[5+n_signal_parameters];
			double p1 = par[6+n_signal_parameters];
			double p2 = par[7+n_signal_parameters];
			double p3 = par[8+n_signal_parameters];
			double p4 = par[9+n_signal_parameters];
			
			f_background = p0 * exp(p1*(loc_mgg - p2) + p3*pow(loc_mgg - p2,2.0)) + p4;
			break;
		}
	}
	
	//-------------------------//
	// enhancement around 0.44 GeV:
	
	double     N_fdc = par[n_signal_parameters+n_background_parameters+5+0];
	double    mu_fdc = par[n_signal_parameters+n_background_parameters+5+1];
	double sigma_fdc = par[n_signal_parameters+n_background_parameters+5+2];
	double     A_fdc = N_fdc * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma_fdc;
	
	Double_t f_fdc = A_fdc * exp(-0.5*pow((loc_mgg-mu_fdc)/sigma_fdc, 2.0));
	
	n_background_parameters += 3;
	
	//-------------------------//
	// eta-prime:
	
	Double_t f_eta_prime;
	
	double     N_etap = par[n_signal_parameters+5+n_background_parameters+0];
	double    mu_etap = par[n_signal_parameters+5+n_background_parameters+1];
	double sigma_etap = par[n_signal_parameters+5+n_background_parameters+2];
	double     A_etap = N_etap * m_mgg_bin_size / sqrt(2.0*TMath::Pi()) / sigma_etap;
	
	f_eta_prime = A_etap*exp(-0.5*pow((loc_mgg-mu_etap)/sigma_etap, 2.0));
	
	//----------------------------------------------------------------------------------//
	
	Double_t f_mgg = f_eta + f_omega + f_background + f_eta_prime + f_fdc;
	return f_mgg;
}
