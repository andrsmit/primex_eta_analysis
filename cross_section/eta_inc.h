#ifndef _ETAINC_
#define _ETAINC_

//---------------------------------------------------//
// Angular Bin Size [default is 0.08deg]:

int m_rebins_theta = 8;
double m_theta_bin_size = 0.01 * (double)m_rebins_theta;

//---------------------------------------------------//
// Mass Bin Size [default is 6MeV]:

int m_rebins_mgg = 5;
double m_mgg_bin_size = (1.2/600.) * (double)m_rebins_mgg;

//---------------------------------------------------//
// Constants:

int    m_phase, m_pass;
double m_luminosity;
double m_empty_target_flux_ratio;

double m_branching_ratio = 0.3936;

double m_eta   = 0.547862;
double m_omega = 0.78266;
double m_etap  = 0.95778;

/*************************************************************************/
// Switches:

bool USE_SOLID_ANGLE   = false;
bool DRAW_FITS         = false;
bool USE_DOUBLE_GAUS   = false;
bool DRAW_ACCEPTANCE   = true;
bool DRAW_ETA_MC_PARS  = true;
bool SIMPLE_BACKGROUND = false;

/*************************************************************************/
// Data:

TString m_root_fname, m_empty_target_root_fname;
TString m_hist_name;

vector<pair<double,double>> m_angular_bin;
vector<pair<double,double>> m_angular_yield;
vector<pair<double,double>> m_angular_yield_empty;

TH2F *h_mgg_vs_theta_full;
TH2F *h_mgg_vs_theta_empty;

double m_min_mgg_cut = 0.50;
double m_max_mgg_cut = 0.60;

/*************************************************************************/
// fitting:

bool m_SUBTRACT_EMPTY      = true;  // subtract empty target background prior to fit
bool m_FIT_FDC_ENHANCEMENT = false; // use a Gaussian to fit the enhancement seen around 0.45 GeV 
                                    // (only checked if m_SUBTRACT_EMPTY==false)

int m_signal_fit_option = 1, m_omega_fit_option = 1, m_background_fit_option = 1;

vector<pair<double,double>> m_exclude_regions;

vector<pair<int,TF1*>> f_mgg_fit_functions;

vector<pair<TString,double>> m_fit_parameters;
vector<pair<TString,double>> m_fit_errors;

TH1F *h_mgg_mc_ls;

bool DRAW_MGG_FITS = true;

TCanvas *c_fits;
TPad *p_fits[12];

TCanvas *cChi2;
TCanvas *cYield;

TGraphErrors *gYield, *gYieldEmpty;

int    m_n_bins_fit;
double m_min_bkgd_fit, m_max_bkgd_fit;

#endif
