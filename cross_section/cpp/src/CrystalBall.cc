#include "CrossSection.h"

double CrystalBall(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     N_omega1 = par[0];
	double    mu_omega1 = par[1];
	double sigma_omega1 = par[2];
	double     a_omega1 = par[3];
	double     n_omega1 = par[4];
	
	double Acb_omega1 = pow(n_omega1/fabs(a_omega1), n_omega1) * exp(-0.5*pow(fabs(a_omega1),2.0));
	double Bcb_omega1 = (n_omega1/fabs(a_omega1)) - fabs(a_omega1);
	
	double loc_x_omega1 = (loc_mgg - mu_omega1)/sigma_omega1;
	
	double f_omega1;
	if(loc_x_omega1 > -a_omega1) {
		f_omega1 = N_omega1 * exp(-0.5*pow(loc_x_omega1,2.0));
	} else {
		f_omega1 = N_omega1 * Acb_omega1 * pow(Bcb_omega1 - loc_x_omega1, -n_omega1);
	}
	
	return f_omega1;
}

double CrystalBall2(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     N_omega1 = par[0];
	double    mu_omega1 = par[1];
	double sigma_omega1 = par[2];
	double     a_omega1 = par[3];
	double     n_omega1 = par[4];
	
	double Acb_omega1 = pow(n_omega1/fabs(a_omega1), n_omega1) * exp(-0.5*pow(fabs(a_omega1),2.0));
	double Bcb_omega1 = (n_omega1/fabs(a_omega1)) - fabs(a_omega1);
	
	double loc_x_omega1 = (loc_mgg - mu_omega1)/sigma_omega1;
	
	double f_omega1;
	if(loc_x_omega1 > -a_omega1) {
		f_omega1 = N_omega1 * exp(-0.5*pow(loc_x_omega1,2.0));
	} else {
		f_omega1 = N_omega1 * Acb_omega1 * pow(Bcb_omega1 - loc_x_omega1, -n_omega1);
	}
	
	double     N_omega2 = par[5];
	double    mu_omega2 = par[6];
	double sigma_omega2 = par[7];
	double     a_omega2 = par[8];
	double     n_omega2 = par[9];
	
	double Acb_omega2 = pow(n_omega2/fabs(a_omega2), n_omega2) * exp(-0.5*pow(fabs(a_omega2),2.0));
	double Bcb_omega2 = (n_omega2/fabs(a_omega2)) - fabs(a_omega2);
	
	double loc_x_omega2 = (loc_mgg - mu_omega2)/sigma_omega2;
	
	double f_omega2;
	if(loc_x_omega2 > -a_omega2) {
		f_omega2 = N_omega2 * exp(-0.5*pow(loc_x_omega2,2.0));
	} else {
		f_omega2 = N_omega2 * Acb_omega2 * pow(Bcb_omega2 - loc_x_omega2, -n_omega2);
	}
	
	return f_omega1+f_omega2;
}
