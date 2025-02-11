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


double CrystalBall_flip(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     N_omega1 = par[0];
	double    mu_omega1 = par[1];
	double sigma_omega1 = par[2];
	double     a_omega1 = par[3];
	double     n_omega1 = par[4];
	
	double Acb_omega1 = pow(n_omega1/fabs(a_omega1), n_omega1) * exp(-0.5*pow(fabs(a_omega1),2.0));
	double Bcb_omega1 = (n_omega1/fabs(a_omega1)) - fabs(a_omega1);
	
	double loc_x_omega1 = (mu_omega1 - loc_mgg)/sigma_omega1;
	
	double f_omega1;
	if(loc_x_omega1 > -a_omega1) {
		f_omega1 = N_omega1 * exp(-0.5*pow(loc_x_omega1,2.0));
	} else {
		f_omega1 = N_omega1 * Acb_omega1 * pow(Bcb_omega1 - loc_x_omega1, -n_omega1);
	}
	
	return f_omega1;
}

double CrystalBall2_flip(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     N_omega1 = par[0];
	double    mu_omega1 = par[1];
	double sigma_omega1 = par[2];
	double     a_omega1 = par[3];
	double     n_omega1 = par[4];
	
	double Acb_omega1 = pow(n_omega1/fabs(a_omega1), n_omega1) * exp(-0.5*pow(fabs(a_omega1),2.0));
	double Bcb_omega1 = (n_omega1/fabs(a_omega1)) - fabs(a_omega1);
	
	double loc_x_omega1 = (mu_omega1 - loc_mgg)/sigma_omega1;
	
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
	
	double loc_x_omega2 = (mu_omega2 - loc_mgg)/sigma_omega2;
	
	double f_omega2;
	if(loc_x_omega2 > -a_omega2) {
		f_omega2 = N_omega2 * exp(-0.5*pow(loc_x_omega2,2.0));
	} else {
		f_omega2 = N_omega2 * Acb_omega2 * pow(Bcb_omega2 - loc_x_omega2, -n_omega2);
	}
	
	return f_omega1+f_omega2;
}

double CrystalBall3_flip(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     N_omega1 = par[0];
	double    mu_omega1 = par[1];
	double sigma_omega1 = par[2];
	double     a_omega1 = par[3];
	double     n_omega1 = par[4];
	
	double Acb_omega1 = pow(n_omega1/fabs(a_omega1), n_omega1) * exp(-0.5*pow(fabs(a_omega1),2.0));
	double Bcb_omega1 = (n_omega1/fabs(a_omega1)) - fabs(a_omega1);
	
	double loc_x_omega1 = (mu_omega1 - loc_mgg)/sigma_omega1;
	
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
	
	double loc_x_omega2 = (mu_omega2 - loc_mgg)/sigma_omega2;
	
	double f_omega2;
	if(loc_x_omega2 > -a_omega2) {
		f_omega2 = N_omega2 * exp(-0.5*pow(loc_x_omega2,2.0));
	} else {
		f_omega2 = N_omega2 * Acb_omega2 * pow(Bcb_omega2 - loc_x_omega2, -n_omega2);
	}
	
	double     N_omega3 = par[10];
	double    mu_omega3 = par[11];
	double sigma_omega3 = par[12];
	double     a_omega3 = par[13];
	double     n_omega3 = par[14];
	
	double Acb_omega3 = pow(n_omega3/fabs(a_omega3), n_omega3) * exp(-0.5*pow(fabs(a_omega3),2.0));
	double Bcb_omega3 = (n_omega3/fabs(a_omega3)) - fabs(a_omega3);
	
	double loc_x_omega3 = (mu_omega3 - loc_mgg)/sigma_omega3;
	
	double f_omega3;
	if(loc_x_omega3 > -a_omega3) {
		f_omega3 = N_omega3 * exp(-0.5*pow(loc_x_omega3,2.0));
	} else {
		f_omega3 = N_omega3 * Acb_omega3 * pow(Bcb_omega3 - loc_x_omega3, -n_omega3);
	}
	
	return f_omega1+f_omega2+f_omega3;
}

double DoubleGaus(double *x, double *par) {
	
	double loc_mgg = x[0];
	
	double     A1 = par[0];
	double    mu1 = par[1];
	double sigma1 = par[2];
	double     A2 = par[3];
	double    mu2 = par[4];
	double sigma2 = par[5];
	
	double f = A1*exp(-0.5*pow((loc_mgg-mu1)/sigma1,2.0)) + A2*exp(-0.5*pow((loc_mgg-mu2)/sigma2,2.0));
	return f;
}
