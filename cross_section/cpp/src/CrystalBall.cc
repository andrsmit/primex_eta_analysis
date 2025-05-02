#include "CrossSection.h"

double DoubleGausPDF(double *x, double *par)
{
	double mgg = x[0];
	
	//---------------------------------------------//
	// First Gaussian:
	
	double    mu1 = par[0];
	double sigma1 = par[1];
	
	double f1 = NormGaus(mgg, mu1, sigma1);
	
	//---------------------------------------------//
	// Second Gaussian:
	
	double    mu2 = par[2] + mu1;
	double sigma2 = par[3];
	
	double f2 = NormGaus(mgg, mu2, sigma2);
	
	//---------------------------------------------//
	// Combine:
	
	double fraction = par[4];
	double binWidth = par[5];
	
	double f = ((1.0-fraction)*f1 + fraction*f2) * binWidth;
	return f;
}

double CrystalBallPDF(double *x, double *par)
{
	double mgg = x[0];
	
	double    mu = par[0];
	double sigma = par[1];
	double alpha = par[2];
	double     n = par[3];
	
	double binWidth = par[4];
	
	double f = binWidth * NormCrystalBall(mgg, mu, sigma, alpha, n);
	return f;
}

double DoubleCrystalBallPDF(double *x, double *par)
{
	double mgg = x[0];
	
	//---------------------------------------------//
	// First CrystalBall:
	
	double    mu1 = par[0];
	double sigma1 = par[1];
	double alpha1 = par[2];
	double     n1 = par[3];
	
	double f1 = NormCrystalBall(mgg, mu1, sigma1, alpha1, n1);
	
	//---------------------------------------------//
	// Second CrystalBall:
	
	double    mu2 = par[4] + mu1;
	double sigma2 = par[5];
	double alpha2 = par[6];
	double     n2 = par[7];
	
	double f2 = NormCrystalBall(mgg, mu2, sigma2, alpha2, n2);
	
	//---------------------------------------------//
	// Combine:
	
	double fraction = par[8];
	double binWidth = par[9];
	
	double f = ((1.0-fraction)*f1 + fraction*f2) * binWidth;
	return f;
}

double CrystalBallPDF_flip(double *x, double *par)
{
	double mgg = x[0];
	
	double    mu = par[0];
	double sigma = par[1];
	double alpha = par[2];
	double     n = par[3];
	
	double binWidth = par[4];
	
	double f = binWidth * NormCrystalBall(mgg, mu, sigma, alpha, n, 1);
	return f;
}

double DoubleCrystalBallPDF_flip(double *x, double *par)
{
	double mgg = x[0];
	
	//---------------------------------------------//
	// First CrystalBall:
	
	double    mu1 = par[0];
	double sigma1 = par[1];
	double alpha1 = par[2];
	double     n1 = par[3];
	
	double f1 = NormCrystalBall(mgg, mu1, sigma1, alpha1, n1, 1);
	
	//---------------------------------------------//
	// Second CrystalBall:
	
	double    mu2 = par[4] + mu1;
	double sigma2 = par[5];
	double alpha2 = par[6];
	double     n2 = par[7];
	
	double f2 = NormCrystalBall(mgg, mu2, sigma2, alpha2, n2, 1);
	
	//---------------------------------------------//
	// Combine:
	
	double fraction = par[8];
	double binWidth = par[9];
	
	double f = ((1.0-fraction)*f1 + fraction*f2) * binWidth;
	return f;
}

double DoubleCrystalBallPDF_oneflip(double *x, double *par)
{
	double mgg = x[0];
	
	//---------------------------------------------//
	// First CrystalBall:
	
	double    mu1 = par[0];
	double sigma1 = par[1];
	double alpha1 = par[2];
	double     n1 = par[3];
	
	double f1 = NormCrystalBall(mgg, mu1, sigma1, alpha1, n1);
	
	//---------------------------------------------//
	// Second CrystalBall:
	
	double    mu2 = par[4] + mu1;
	double sigma2 = par[5];
	double alpha2 = par[6];
	double     n2 = par[7];
	
	double f2 = NormCrystalBall(mgg, mu2, sigma2, alpha2, n2, 1);
	
	//---------------------------------------------//
	// Combine:
	
	double fraction = par[8];
	double binWidth = par[9];
	
	double f = ((1.0-fraction)*f1 + fraction*f2) * binWidth;
	return f;
}

double DoubleCrystalBallPlusGausPDF(double *x, double *par)
{
	double mgg = x[0];
	
	//---------------------------------------------//
	// First CrystalBall:
	
	double    mu1 = par[0];
	double sigma1 = par[1];
	double alpha1 = par[2];
	double     n1 = par[3];
	
	double f1 = NormCrystalBall(mgg, mu1, sigma1, alpha1, n1);
	
	//---------------------------------------------//
	// Second CrystalBall:
	
	double    mu2 = par[4] + mu1;
	double sigma2 = par[5];
	double alpha2 = par[6];
	double     n2 = par[7];
	
	double f2 = NormCrystalBall(mgg, mu2, sigma2, alpha2, n2, 1);
	
	//---------------------------------------------//
	// Gaussian:
	
	double    mu3 = par[8];
	double sigma3 = par[9];
	
	double f3 = NormGaus(mgg, mu3, sigma3);
	
	//---------------------------------------------//
	// Combine:
	
	double fraction1 = par[10];
	double fraction2 = par[11];
	double binWidth  = par[12];
	
	if((fraction1+fraction2)>1.0) return -1.e9;
	
	double f = ((1.0-fraction1-fraction2)*f1 + fraction1*f2 + fraction2*f3) * binWidth;
	return f;
}

double NormGaus(double x, double mu, double sigma)
{
	// Returns the evaluation of a normalized Gaussian PDF:
	double A = 1.0 / (sqrt(2.0*TMath::Pi()) * sigma);
	double f = A * exp(-0.5*pow((x-mu)/sigma,2.0));
	return f;
}

double NormCrystalBall(double x, double mu, double sigma, double alpha, double n, int doFlip)
{
	// Returns the evaluation of a normalized Crystal Ball PDF:
	
	double alpha_abs = fabs(alpha);
	
	double A = pow(n/alpha_abs, n) * exp(-0.5*pow(alpha_abs,2.0));
	double B = (n/alpha_abs) - alpha_abs;
	
	double locX = doFlip ? (mu - x)/sigma : (x - mu)/sigma;
	
	double f;
	if(locX > -alpha) {
		f = exp(-0.5*pow(locX,2.0));
	} else {
		f = A * pow(B-locX, -n);
	}
	
	// Normalization:
	
	double C = (n/alpha_abs)*(1.0/(n-1.0)) * exp(-0.5*pow(alpha_abs,2.0));
	double D = sqrt(TMath::Pi()/2.0) * (1.0 + erf(alpha_abs/sqrt(2.0)));
	double N = 1.0 / (sigma * (C+D));
	return (N * f);
}
