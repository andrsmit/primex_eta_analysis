#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"

void initializeFitParameters(TF1 *f1);

void   guessOmegaParameters(TF1 *f1, TH1F *h1);
void     fixOmegaParameters(TF1 *f1);
void releaseOmegaParameters(TF1 *f1);

void   guessBackgroundParameters(TF1 *f1, TH1F *h1);
void     fixBackgroundParameters(TF1 *f1);
void releaseBackgroundParameters(TF1 *f1);

void   guessFDCParameters(TF1 *f1, TH1F *h1);
void     fixFDCParameters(TF1 *f1);
void releaseFDCParameters(TF1 *f1);

void   guessEtapParameters(TF1 *f1, TH1F *h1);
void     fixEtapParameters(TF1 *f1);
void releaseEtapParameters(TF1 *f1);

void   guessEtaParameters(TF1 *f1, TH1F *h1);
void     fixEtaParameters(TF1 *f1);
void releaseEtaParameters(TF1 *f1);

void fit_mgg(TH1F *h1, TF1 *f_fit, int n_parameters) {
	
	// initially fix all parameters to zero:
	
	initializeFitParameters(f_fit);
	
	//-----------------------------------------------------------------//
	// release omega parameters and fit the region around the peak:
	
	double min_omega_fit = 0.65;
	double max_omega_fit = 0.85;
	
	// get initial guesses for omega fit parameters:
	
	guessOmegaParameters(f_fit, h1);
	
	// restrict range and fit omega parameters:
	
	f_fit->SetRange(min_omega_fit, max_omega_fit);
	h1->Fit(f_fit, "R0QL");
	
	// fix omega fit parameters and widen fit range:
	
	//fixOmegaParameters(f_fit);
	f_fit->SetRange(m_min_bkgd_fit, 0.90);
	
	// get initial guesses for background fit parameters:
	
	guessBackgroundParameters(f_fit, h1);
	if(!m_SUBTRACT_EMPTY && m_FIT_FDC_ENHANCEMENT) {
		guessFDCParameters(f_fit, h1);
	}
	
	// exlude eta mass region from fit:
	
	m_exclude_regions.clear();
	m_exclude_regions.push_back({0.525, 0.665});
	
	h1->Fit(f_fit, "R0QL");
	
	// fit eta':
	
	f_fit->SetRange(m_min_bkgd_fit, m_max_bkgd_fit);
	if(m_max_bkgd_fit > m_etap) {
		guessEtapParameters(f_fit, h1);
		h1->Fit(f_fit, "R0QL");
	}
	
	//-----------------------------------------------------------------//
	// fix the above parameters and fit the eta peak:
	
	fixOmegaParameters(f_fit);
	fixBackgroundParameters(f_fit);
	if(!m_SUBTRACT_EMPTY && m_FIT_FDC_ENHANCEMENT) fixFDCParameters(f_fit);
	fixEtapParameters(f_fit);
	
	m_exclude_regions.clear();
	guessEtaParameters(f_fit, h1);
	
	f_fit->SetRange(0.5, 0.625);
	h1->Fit(f_fit, "R0QL");
	
	releaseOmegaParameters(f_fit);
	releaseBackgroundParameters(f_fit);
	if(!m_SUBTRACT_EMPTY && m_FIT_FDC_ENHANCEMENT) releaseFDCParameters(f_fit);
	
	f_fit->SetRange(m_min_bkgd_fit, m_max_bkgd_fit);
	h1->Fit(f_fit, "R0QL");
	
	return;
}

void initializeFitParameters(TF1 *f1) {
	
	// eta fit parameters:
	int loc_n_fit_parameters = 0;
	switch(m_signal_fit_option) {
		case 1:
			loc_n_fit_parameters = 3;
			f1->FixParameter(0, 0.0);
			f1->FixParameter(1, m_eta);
			f1->FixParameter(2, 0.02);
			break;
		case 2:
			loc_n_fit_parameters = 6;
			f1->FixParameter(0, 0.0);
			f1->FixParameter(1, 0.0);
			f1->FixParameter(2, m_eta);
			f1->FixParameter(3, 0.0);
			f1->FixParameter(4, 0.02);
			f1->FixParameter(5, 0.02);
			break;
		case 3:
			loc_n_fit_parameters = 5;
			f1->FixParameter(0, 0.0);
			f1->FixParameter(1, m_eta);
			f1->FixParameter(2, 0.02);
			f1->FixParameter(3, 1.0);
			f1->FixParameter(4, 1.0);
			break;
		case 4:
			loc_n_fit_parameters = 8;
			f1->FixParameter(0, 0.0);
			f1->FixParameter(1, m_eta);
			f1->FixParameter(2, 0.02);
			f1->FixParameter(3, 1.0);
			f1->FixParameter(4, 1.0);
			f1->FixParameter(5, 0.0);
			f1->FixParameter(6, m_eta);
			f1->FixParameter(7, 0.02);
			break;
	}
	
	f1->FixParameter(loc_n_fit_parameters+0, 0.0);
	f1->FixParameter(loc_n_fit_parameters+1, m_omega);
	f1->FixParameter(loc_n_fit_parameters+2, 0.02);
	f1->FixParameter(loc_n_fit_parameters+3, 1.0);
	f1->FixParameter(loc_n_fit_parameters+4, 1.0);
	loc_n_fit_parameters += 5;
	
	// background fit parameters:
	switch(m_background_fit_option) {
		case 1:
			for(int ipar=0; ipar<4; ipar++) f1->FixParameter(loc_n_fit_parameters+ipar, 0.0);
			loc_n_fit_parameters += 4;
			break;
		case 2:
			for(int ipar=0; ipar<5; ipar++) f1->FixParameter(loc_n_fit_parameters+ipar, 0.0);
			loc_n_fit_parameters += 5;
			break;
	}
	
	// fdc-enhancement:
	f1->FixParameter(loc_n_fit_parameters+0, 0.0);
	f1->FixParameter(loc_n_fit_parameters+1, 0.45);
	f1->FixParameter(loc_n_fit_parameters+2, 0.015);
	loc_n_fit_parameters += 3;
	
	// eta-prime fit parameters:
	f1->FixParameter(loc_n_fit_parameters+0, 0.0);
	f1->FixParameter(loc_n_fit_parameters+1, m_etap);
	f1->FixParameter(loc_n_fit_parameters+2, 0.02);
	
	return;
}

//--------------------------------------------------------------//
// Omega->pi0+gamma peak:

void guessOmegaParameters(TF1 *f1, TH1F *h1) {
	
	double min_omega_fit = 0.75;
	double max_omega_fit = 0.85;
	
	double     N_omega_guess = h1->Integral(h1->FindBin(min_omega_fit), h1->FindBin(max_omega_fit));
	double    mu_omega_guess = m_omega;
	double sigma_omega_guess = 0.025;
	double alpha_omega_guess = 1.0;
	double     n_omega_guess = 2.0;
	
	double loc_omega_max = 0.0;
	for(int ibin=h1->FindBin(min_omega_fit); ibin<=h1->FindBin(max_omega_fit); ibin++) {
		if(h1->GetBinContent(ibin) > loc_omega_max) {
			loc_omega_max = h1->GetBinContent(ibin);
			mu_omega_guess = h1->GetBinCenter(ibin);
		}
	}
	
	int     N_omega_par = f1->GetParNumber("N_{#omega}");
	int    mu_omega_par = f1->GetParNumber("#mu_{#omega}");
	int sigma_omega_par = f1->GetParNumber("#sigma_{#omega}");
	int alpha_omega_par = f1->GetParNumber("#alpha_{#omega}");
	int     n_omega_par = f1->GetParNumber("n_{#omega}");
	
	f1->SetParameter(    N_omega_par,     N_omega_guess);
	f1->SetParameter(   mu_omega_par,    mu_omega_guess);
	f1->SetParameter(sigma_omega_par, sigma_omega_guess);
	f1->SetParameter(alpha_omega_par, alpha_omega_guess);
	f1->SetParameter(    n_omega_par,     n_omega_guess);
	
	f1->SetParLimits(    N_omega_par, 0.,	1.e6);
	f1->SetParLimits(   mu_omega_par, 0.750, 0.800);
	f1->SetParLimits(sigma_omega_par, 0.015, 0.050);
	f1->SetParLimits(alpha_omega_par, 0.500, 9.999);
	f1->SetParLimits(    n_omega_par, 0.500, 9.999);
	
	return;
}

void fixOmegaParameters(TF1 *f1) {
	
	int     N_omega_par = f1->GetParNumber("N_{#omega}");
	int    mu_omega_par = f1->GetParNumber("#mu_{#omega}");
	int sigma_omega_par = f1->GetParNumber("#sigma_{#omega}");
	int alpha_omega_par = f1->GetParNumber("#alpha_{#omega}");
	int     n_omega_par = f1->GetParNumber("n_{#omega}");
	
	f1->FixParameter(    N_omega_par, f1->GetParameter(    N_omega_par));
	f1->FixParameter(   mu_omega_par, f1->GetParameter(   mu_omega_par));
	f1->FixParameter(sigma_omega_par, f1->GetParameter(sigma_omega_par));
	f1->FixParameter(alpha_omega_par, f1->GetParameter(alpha_omega_par));
	f1->FixParameter(    n_omega_par, f1->GetParameter(    n_omega_par));
	
	return;
}

void releaseOmegaParameters(TF1 *f1) {
	
	int     N_omega_par = f1->GetParNumber("N_{#omega}");
	int    mu_omega_par = f1->GetParNumber("#mu_{#omega}");
	int sigma_omega_par = f1->GetParNumber("#sigma_{#omega}");
	int alpha_omega_par = f1->GetParNumber("#alpha_{#omega}");
	int     n_omega_par = f1->GetParNumber("n_{#omega}");
	
	f1->ReleaseParameter(N_omega_par);
	f1->SetParLimits(N_omega_par, 0., 1.e6);
	
	f1->ReleaseParameter(mu_omega_par);
	f1->SetParLimits(mu_omega_par, 0.750, 0.800);
	
	f1->ReleaseParameter(sigma_omega_par);
	f1->SetParLimits(sigma_omega_par, 0.015, 0.050);
	
	f1->ReleaseParameter(alpha_omega_par);
	f1->SetParLimits(alpha_omega_par, 0.500, 9.999);
	
	f1->ReleaseParameter(n_omega_par);
	f1->SetParLimits(n_omega_par, 0.500, 9.999);
	
	return;
}

//--------------------------------------------------------------//
// Background:

void guessBackgroundParameters(TF1 *f1, TH1F *h1) {
	
	int p0_par = f1->GetParNumber("p0");
	int p1_par = f1->GetParNumber("p1");
	int p2_par = f1->GetParNumber("p2");
	int p3_par = f1->GetParNumber("p3");
	
	double p0_guess = h1->GetBinContent(h1->FindBin(m_min_bkgd_fit)) - f1->Eval(m_min_bkgd_fit);
	double p1_guess = 0., p2_guess = 0., p3_guess = 0., p4_guess = 0.;
	
	if(m_background_fit_option==2) {
		p1_guess = -1.0;
		p2_guess =  m_min_bkgd_fit;
	}
	
	switch(m_background_fit_option) {
		case 1:
		{
			// 3rd-order polynomial:
			
			f1->SetParameter(p0_par, p0_guess);
			f1->SetParameter(p1_par, p1_guess);
			f1->SetParameter(p2_par, p2_guess);
			f1->SetParameter(p3_par, p3_guess);
			
			f1->SetParLimits(p0_par, -1.e4, 1.e4);
			f1->SetParLimits(p1_par, -1.e3, 1.e3);
			f1->SetParLimits(p2_par, -1.e3, 1.e3);
			f1->SetParLimits(p3_par, -1.e3, 1.e3);
			break;
		}
		case 2:
		{
			// Exponential:
			int p4_par = f1->GetParNumber("p4");
			
			f1->SetParameter(p0_par, p0_guess);
			f1->SetParameter(p1_par, p1_guess);
			f1->SetParameter(p2_par, p2_guess);
			//f1->SetParameter(p3_par, p3_guess);
			f1->SetParameter(p4_par, p4_guess);
			
			f1->SetParLimits(p0_par, 0.0, 1.e4);
			f1->SetParLimits(p1_par, -1.e3, 0.0);
			f1->SetParLimits(p2_par, 0.0, 1.0);
			//f1->SetParLimits(p3_par, -1.e3, 1.e3);
			f1->SetParLimits(p4_par, 0.0, 1.e3);
			break;
		}
	}
	
	return;
}

void fixBackgroundParameters(TF1 *f1) {
	
	int p0_par = f1->GetParNumber("p0");
	int p1_par = f1->GetParNumber("p1");
	int p2_par = f1->GetParNumber("p2");
	int p3_par = f1->GetParNumber("p3");
	
	f1->FixParameter(p0_par, f1->GetParameter(p0_par));
	f1->FixParameter(p1_par, f1->GetParameter(p1_par));
	f1->FixParameter(p2_par, f1->GetParameter(p2_par));
	f1->FixParameter(p3_par, f1->GetParameter(p3_par));
	if(m_background_fit_option==2) {
		int p4_par = f1->GetParNumber("p4");
		f1->FixParameter(p4_par, f1->GetParameter(p4_par));
	}
	
	return;
}

void releaseBackgroundParameters(TF1 *f1) {
	
	int p0_par = f1->GetParNumber("p0");
	int p1_par = f1->GetParNumber("p1");
	int p2_par = f1->GetParNumber("p2");
	int p3_par = f1->GetParNumber("p3");
	
	switch(m_background_fit_option) {
		case 1:
		{
			// 3rd-order polynomial:
			
			f1->ReleaseParameter(p0_par);
			f1->ReleaseParameter(p1_par);
			f1->ReleaseParameter(p2_par);
			f1->ReleaseParameter(p3_par);
			
			f1->SetParLimits(p0_par, -1.e4, 1.e4);
			f1->SetParLimits(p1_par, -1.e3, 1.e3);
			f1->SetParLimits(p2_par, -1.e3, 1.e3);
			f1->SetParLimits(p3_par, -1.e3, 1.e3);
			break;
		}
		case 2:
		{
			// Exponential:
			int p4_par = f1->GetParNumber("p4");
			
			f1->ReleaseParameter(p0_par);
			f1->ReleaseParameter(p1_par);
			f1->ReleaseParameter(p2_par);
			//f1->ReleaseParameter(p3_par);
			f1->ReleaseParameter(p4_par);
			
			f1->SetParLimits(p0_par, 0.0, 1.e4);
			f1->SetParLimits(p1_par, -1.e3, 0.0);
			f1->SetParLimits(p2_par, 0.0, 1.0);
			//f1->SetParLimits(p3_par, -1.e3, 1.e3);
			f1->SetParLimits(p4_par, 0.0, 1.e3);
			break;
		}
	}
	
	return;
}

//--------------------------------------------------------------//
// Enhancement from FDC packages around 0.45 GeV:

void guessFDCParameters(TF1 *f1, TH1F *h1) {
	
	int     N_fdc_par = f1->GetParNumber("N_{FDC}");
	int    mu_fdc_par = f1->GetParNumber("#mu_{FDC}");
	int sigma_fdc_par = f1->GetParNumber("#sigma_{FDC}");
	
	// just set N to 0:
	
	double     N_fdc_guess = 0.0;
	double    mu_fdc_guess = 0.45;
	double sigma_fdc_guess = 0.015;
	
	f1->SetParameter(    N_fdc_par,     N_fdc_guess);
	f1->SetParameter(   mu_fdc_par,    mu_fdc_guess);
	f1->SetParameter(sigma_fdc_par, sigma_fdc_guess);
	
	f1->SetParLimits(    N_fdc_par, 0.000, 1.0e4);
	f1->SetParLimits(   mu_fdc_par, 0.420, 0.475);
	f1->SetParLimits(sigma_fdc_par, 0.015, 0.050);
	
	return;
}

void fixFDCParameters(TF1 *f1) {
	
	int     N_fdc_par = f1->GetParNumber("N_{FDC}");
	int    mu_fdc_par = f1->GetParNumber("#mu_{FDC}");
	int sigma_fdc_par = f1->GetParNumber("#sigma_{FDC}");
	
	f1->FixParameter(    N_fdc_par, f1->GetParameter(N_fdc_par));
	f1->FixParameter(   mu_fdc_par, f1->GetParameter(mu_fdc_par));
	f1->FixParameter(sigma_fdc_par, f1->GetParameter(sigma_fdc_par));
	
	return;
}

void releaseFDCParameters(TF1 *f1) {
	
	int     N_fdc_par = f1->GetParNumber("N_{FDC}");
	int    mu_fdc_par = f1->GetParNumber("#mu_{FDC}");
	int sigma_fdc_par = f1->GetParNumber("#sigma_{FDC}");
	
	f1->ReleaseParameter(    N_fdc_par);
	f1->ReleaseParameter(   mu_fdc_par);
	f1->ReleaseParameter(sigma_fdc_par);
	
	f1->SetParLimits(    N_fdc_par, 0.000, 1.0e4);
	f1->SetParLimits(   mu_fdc_par, 0.420, 0.475);
	f1->SetParLimits(sigma_fdc_par, 0.015, 0.050);
	
	return;
}

//--------------------------------------------------------------//
// Eta prime:

void guessEtapParameters(TF1 *f1, TH1F *h1) {
	
	int     N_etap_par = f1->GetParNumber("N_{#eta'}");
	int    mu_etap_par = f1->GetParNumber("#mu_{#eta'}");
	int sigma_etap_par = f1->GetParNumber("#sigma_{#eta'}");
	
	// guess number of eta' by integrating histogram and subtracting background:
	
	double min_etap_fit = 0.91;
	double max_etap_fit = 1.02;
	
	double     N_etap_guess = h1->Integral(h1->FindBin(min_etap_fit), h1->FindBin(max_etap_fit));
	double    mu_etap_guess = m_etap;
	double sigma_etap_guess = 0.025;
	
	f1->SetParameter(    N_etap_par,     N_etap_guess);
	f1->SetParameter(   mu_etap_par,    mu_etap_guess);
	f1->SetParameter(sigma_etap_par, sigma_etap_guess);
	
	f1->SetParLimits(    N_etap_par, 0.,	1.e4);
	f1->SetParLimits(   mu_etap_par, 0.920, 0.990);
	f1->SetParLimits(sigma_etap_par, 0.015, 0.050);
	
	return;
}

void fixEtapParameters(TF1 *f1) {
	
	int     N_etap_par = f1->GetParNumber("N_{#eta'}");
	int    mu_etap_par = f1->GetParNumber("#mu_{#eta'}");
	int sigma_etap_par = f1->GetParNumber("#sigma_{#eta'}");
	
	f1->FixParameter(    N_etap_par, f1->GetParameter(N_etap_par));
	f1->FixParameter(   mu_etap_par, f1->GetParameter(mu_etap_par));
	f1->FixParameter(sigma_etap_par, f1->GetParameter(sigma_etap_par));
	
	return;
}

//--------------------------------------------------------------//
// Eta:

void guessEtaParameters(TF1 *f1, TH1F *h1) {
	
	// guess number of eta' by integrating histogram and subtracting background:
	
	double min_eta_fit = 0.50;
	double max_eta_fit = 0.60;
	
	double     N_eta_guess = h1->Integral(h1->FindBin(min_eta_fit), h1->FindBin(max_eta_fit));
	double    mu_eta_guess = m_eta;
	double sigma_eta_guess = 0.015;
	
	double loc_eta_max = 0.0;
	for(int ibin=h1->FindBin(min_eta_fit); ibin<=h1->FindBin(max_eta_fit); ibin++) {
		if(h1->GetBinContent(ibin) > loc_eta_max) {
			loc_eta_max = h1->GetBinContent(ibin);
			mu_eta_guess = h1->GetBinCenter(ibin);
		}
	}
	
	switch(m_signal_fit_option) {
		case 1:
		{
			// single Gaussian:
			
			int     N_eta_par = f1->GetParNumber("N_{#eta}");
			int    mu_eta_par = f1->GetParNumber("#mu_{#eta}");
			int sigma_eta_par = f1->GetParNumber("#sigma_{#eta}");
			
			f1->SetParameter(    N_eta_par,     N_eta_guess);
			f1->SetParameter(   mu_eta_par,    mu_eta_guess);
			f1->SetParameter(sigma_eta_par, sigma_eta_guess);
			
			f1->SetParLimits(    N_eta_par, 0.,   1.e5);
			f1->SetParLimits(   mu_eta_par, 0.54, 0.62);
			f1->SetParLimits(sigma_eta_par, 0.01, 0.03);
			
			break;
		}
		case 2:
		{
			// double Gaussian:
			
			int     N1_eta_par = f1->GetParNumber("N_{#eta,1}");
			int     N2_eta_par = f1->GetParNumber("N_{#eta,2}");
			int    mu1_eta_par = f1->GetParNumber("#mu_{#eta,1}");
			int    dmu_eta_par = f1->GetParNumber("#mu_{#eta,2}-#mu_{#eta,1}");
			int sigma1_eta_par = f1->GetParNumber("#sigma_{#eta,1}");
			int sigma2_eta_par = f1->GetParNumber("#sigma_{#eta,2}");
			
			f1->SetParameter(    N1_eta_par, 0.9*N_eta_guess);
			f1->SetParameter(    N2_eta_par, 0.1*N_eta_guess);
			f1->SetParameter(   mu1_eta_par, mu_eta_guess);
			f1->SetParameter(   dmu_eta_par, 0.0);
			f1->SetParameter(sigma1_eta_par, sigma_eta_guess);
			f1->SetParameter(sigma2_eta_par, 2.0*sigma_eta_guess);
			
			f1->SetParLimits(    N1_eta_par,  0.00, 1.e5);
			f1->SetParLimits(    N2_eta_par,  0.00, 0.5*N_eta_guess);
			f1->SetParLimits(   mu1_eta_par,  0.54, 0.62);
			f1->SetParLimits(   dmu_eta_par, -0.05, 0.05);
			f1->SetParLimits(sigma1_eta_par,  0.01, 0.03);
			f1->SetParLimits(sigma2_eta_par,  0.01, 0.05);
			
			break;
		}
		case 3:
		{
			// Crsytal ball:
			
			int     N_eta_par = f1->GetParNumber("N_{#eta}");
			int    mu_eta_par = f1->GetParNumber("#mu_{#eta}");
			int sigma_eta_par = f1->GetParNumber("#sigma_{#eta}");
			int alpha_eta_par = f1->GetParNumber("#alpha_{#eta}");
			int     n_eta_par = f1->GetParNumber("n_{#eta}");
			
			f1->SetParameter(    N_eta_par,     N_eta_guess);
			f1->SetParameter(   mu_eta_par,    mu_eta_guess);
			f1->SetParameter(sigma_eta_par, sigma_eta_guess);
			f1->SetParameter(alpha_eta_par,             1.0);
			f1->SetParameter(    n_eta_par,             2.0);
			
			f1->SetParLimits(    N_eta_par, 0.00,  1.e5);
			f1->SetParLimits(   mu_eta_par, 0.54,  0.62);
			f1->SetParLimits(sigma_eta_par, 0.01,  0.03);
			f1->SetParLimits(alpha_eta_par, 0.01,  9.99);
			f1->SetParLimits(    n_eta_par, 0.01, 99.99);
			
			break;
		}
		case 4:
		{
			
			break;
		}
	}
	
	return;
}

void fixEtaParameters(TF1 *f1) {
	
	int     N_etap_par = f1->GetParNumber("N_{#eta'}");
	int    mu_etap_par = f1->GetParNumber("#mu_{#eta'}");
	int sigma_etap_par = f1->GetParNumber("#sigma_{#eta'}");
	
	f1->FixParameter(    N_etap_par, f1->GetParameter(N_etap_par));
	f1->FixParameter(   mu_etap_par, f1->GetParameter(mu_etap_par));
	f1->FixParameter(sigma_etap_par, f1->GetParameter(sigma_etap_par));
	
	return;
}
