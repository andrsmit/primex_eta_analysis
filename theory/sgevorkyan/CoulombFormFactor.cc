#include "PrimExCS.h"
#include "FormFactor.h"

double integrationRadius;
double r_eff, coeff_ch, chargeDensityC1, chargeDensityC2;

TF1 *f_CoulombFFIntegrand_B[2];
double CoulombFFIntegrationRange_b;

TF1 *f_CoulombFFIntegrand_Z[2], *f_CoulombFFIntegrand_Z_simple[2];
double CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2;

double CoulombFFIntegrand_B(double *x, double *par) {
	
	double b = x[0];
	
	double q      = par[0];
	double Delta  = par[1];
	double isReal = par[2];
	
	double kinFactor = 2.0*TMath::Pi()*(pow(q,2.0)+pow(Delta,2.0))/q;
	double zIntegral = CoulombFFIntegrate_Z(isReal, q, Delta, b);
	
	double B = zIntegral * kinFactor;
	return B;
}

double CoulombFFIntegrand_Z(double *x, double *par) {
	
	double z = x[0];
	
	double q      = par[0];
	double Delta  = par[1];
	double b      = par[2];
	double isReal = par[3];
	
	double Pch = GetChargeDensityIntegral(b,z) / pow(pow(b,2.0)+pow(z,2.0),1.5);
	std::complex<double> fAbs = IntegrateAbsorption(b, z);
	
	double Z;
	if(isReal>=0.0) {
		Z = Pch * (cos(Delta*z)*real(fAbs) - sin(Delta*z)*imag(fAbs));
	} else {
		Z = Pch * (sin(Delta*z)*real(fAbs) + cos(Delta*z)*imag(fAbs));
	}
	double besselJ1 = pow(b,2.0) * ROOT::Math::cyl_bessel_j(1.0, q*b);
	return Z*besselJ1;
}

double CoulombFFIntegrand_Z_simple(double *x, double *par) {
	
	double z = x[0];
	
	double q      = par[0];
	double Delta  = par[1];
	double b      = par[2];
	double isReal = par[3];
	
	double trig_func = isReal>=0.0 ? cos(Delta*z) : sin(Delta*z);
	
	double Z = trig_func/pow(pow(b,2.0)+pow(z,2.0),1.5);
	double besselJ1 = pow(b,2.0) * ROOT::Math::cyl_bessel_j(1.0, q*b);
	return Z*besselJ1;
}

long double CalcDiff_StruveBessel(long double t) {
	
	// Code from I. Larin to calculate the difference between
	// the modified Struve function, L_-1, and bessel function, I_1.
	
	int n = 0;
	
	const long double PI = 3.141592653589793238;
	
	long double as0 = 2. / PI;
	long double asn = as0;
	long double ai0 = 0.5*t;
	long double ain = ai0;
	long double sum = as0-ai0;
	long double y   =  0.25*t*t;
	long double l05 = 0.5;
	long double l1  = 1.;
	n = 0;
	while(asn/as0>1.e-20) {
		++n;
		asn *= y / ((long double) (n)-l05) / ((long double) (n)+l05);
		ain *= y / (long double) (n) / (long double) (n+1);
		sum += asn - ain;
	}
	//cout << setprecision(16) << n << " " << sum << endl;
	
	return sum;
}

double CoulombFFIntegrate_Z(double isReal, double q, double Delta, double b) {
	
	// isReal >= 0 -> Real part of integration
	// isReal <  0 -> Imaginary part
	
	int index = isReal >= 0.0 ? 0 : 1;
	
	f_CoulombFFIntegrand_Z[index]->SetParameter(0, q);
	f_CoulombFFIntegrand_Z[index]->SetParameter(1, Delta);
	f_CoulombFFIntegrand_Z[index]->SetParameter(2, b);
	
	for(int i=0; i<2; i++) {
		f_CoulombFFIntegrand_Z_simple[i]->SetParameter(0, q);
		f_CoulombFFIntegrand_Z_simple[i]->SetParameter(1, Delta);
		f_CoulombFFIntegrand_Z_simple[i]->SetParameter(2, b);
	}
	
	// Get the z-independent contribution of the absorption:
	
	std::complex<double> fAbs = IntegrateAbsorption(b, -integrationRadius);
	
	// Split the integration into 3 parts:
	//  1. z between -infinity and -R
	//  2. z between -R and +R
	//  3. z between +R and +infinity
	
	double besselJ1 = pow(b,2.0) * ROOT::Math::cyl_bessel_j(1.0, q*b);
	
	double PchInt = GetChargeDensityIntegral(10.0, 10.0);
	
	double ReZ = (Delta/b) * ROOT::Math::cyl_bessel_k(1.0, Delta*b) * besselJ1 
		- f_CoulombFFIntegrand_Z_simple[0]->Integral(0.0, CoulombFFIntegrationRange_z2);
	
	double longArg = Delta*b;
	long double struveDiff;
	if(longArg > 20.0) struveDiff = CalcDiff_StruveBessel(20.0);
	else struveDiff = CalcDiff_StruveBessel(longArg);
	
	double ImZ = (Delta/b) * besselJ1 * struveDiff
		- f_CoulombFFIntegrand_Z_simple[1]->Integral(0.0, CoulombFFIntegrationRange_z2);
	
	double integral1;
	if(isReal>=0.0) integral1 = PchInt * ((1.0 + real(fAbs))*ReZ + imag(fAbs)*ImZ);
	else            integral1 = PchInt * ((1.0 - real(fAbs))*ImZ + imag(fAbs)*ReZ);
	
	double integral2 = f_CoulombFFIntegrand_Z[index]->Integral(CoulombFFIntegrationRange_z1, 
		CoulombFFIntegrationRange_z2, 1.e-4);
	
	return integral1 + integral2;
}

//-----------------------------------------------------------------------//

double GetChargeDensityIntegral(double b, double z) {
	
	// This function returns the analytical expression resulting from integrating
	// r^2 * charge density from 0 to sqrt(b^2+z^2)
	
	double integral = 0.0;
	switch(modelParameters.densityModel) {
		case 0:
		{
			// Harmonic Oscillator Density Model:
			
			double x = sqrt(pow(b,2.0)+pow(z,2.0))/r_eff;
			double erfx = std::erf(x);
			
			double e_x2p = exp( pow(x,2.0));
			double e_x2n = exp(-pow(x,2.0));
			
			double factor1 = SQRT_PI*erfx - 2.0*e_x2n*x;
			double factor2 = 3.0*SQRT_PI*erfx - 4.0*pow(x,3.0)*e_x2n - 6.0*x*e_x2n;
			
			integral = coeff_ch*(chargeDensityC1*factor1 + chargeDensityC2*factor2);
			break;
		}
		case 1:
		{
			// Sum-of-Gaussians parameterization:
			
			double x = sqrt(pow(b,2.0)+pow(z,2.0));
			
			for(int ipar=0; ipar<modelParameters.nParametersSOG; ipar++) {
				
				double Ri = modelParameters.densityParametersSOG[ipar].first;
				double Ai = modelParameters.densityParametersSOG[ipar].second 
					/ (2.0*pow(SQRT_PI*modelParameters.gammaSOG,3.0)*(1.0 + 2.0*pow(Ri/modelParameters.gammaSOG,2.0)));
				if(Ai == 0.0) continue;
				
				double erfp = std::erf((x+Ri)/modelParameters.gammaSOG);
				double erfm = std::erf((x-Ri)/modelParameters.gammaSOG);
				
				double f1 = SQRT_PI*(1.0+2.0*pow(Ri/modelParameters.gammaSOG,2.0))*(erfm+erfp);
				double f2 = exp(-pow((x+Ri)/modelParameters.gammaSOG,2.0))*(x-Ri)/modelParameters.gammaSOG;
				double f3 = exp(-pow((x-Ri)/modelParameters.gammaSOG,2.0))*(x+Ri)/modelParameters.gammaSOG;
				integral += Ai * (f1 - 2.0*(f2+f3));
			}
			integral *= (pow(modelParameters.gammaSOG,3.0)/4.0);
			break;
		}
	}
	
	return integral;
}

//-----------------------------------------------------------------------//

void CalculateCoulombFF(PrimExCS &csObj) {
	
	r_eff           = sqrt(pow(modelParameters.targetRadius,2.0) + pow(PrimExCS::m_protonRadius,2.0));	
	coeff_ch        = 1.0 / pow(TMath::Pi(), 1.5) / modelParameters.targetA;
	chargeDensityC1 = 1.0 + 0.25*(modelParameters.targetA-4.0)*pow(PrimExCS::m_protonRadius/r_eff,2.0);
	chargeDensityC2 = ((modelParameters.targetA-4.0)/12.0)*pow(modelParameters.targetRadius/r_eff,2.0);
	
	// Set the maximum radius of integration over the nuclear density distribution:
	//    After this point, we make the assumption that the density = 0.
	integrationRadius = 20.0 * modelParameters.targetRadius;
	
	CoulombFFIntegrationRange_b  = 3.0e2;
	CoulombFFIntegrationRange_z1 = -integrationRadius;
	CoulombFFIntegrationRange_z2 =  integrationRadius;
	
	//----------------------------------------//
	// Initialize the function which will be integrated over the variable b:
	
	f_CoulombFFIntegrand_B[0] = new TF1("CoulombFFIntegrandB_real", CoulombFFIntegrand_B, 0.0, CoulombFFIntegrationRange_b, 3);
	f_CoulombFFIntegrand_B[1] = new TF1("CoulombFFIntegrandB_imag", CoulombFFIntegrand_B, 0.0, CoulombFFIntegrationRange_b, 3);
	f_CoulombFFIntegrand_B[0]->SetParameter(2, 1.0);
	f_CoulombFFIntegrand_B[1]->SetParameter(2,-1.0);
	
	//----------------------------------------//
	// Initialize the functions which will be integrated over the variable z:
	
	f_CoulombFFIntegrand_Z[0] = new TF1("CoulombFFIntegrand_Z_real", CoulombFFIntegrand_Z, 
		CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2, 4);
	f_CoulombFFIntegrand_Z[1] = new TF1("CoulombFFIntegrand_Z_imag", CoulombFFIntegrand_Z, 
		CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2, 4);
	f_CoulombFFIntegrand_Z[0]->SetParameter(3, 1.0);
	f_CoulombFFIntegrand_Z[1]->SetParameter(3,-1.0);
	
	f_CoulombFFIntegrand_Z_simple[0] = new TF1("CoulombFFIntegrand_Z_real_simple", 
		CoulombFFIntegrand_Z_simple, CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2, 4);
	f_CoulombFFIntegrand_Z_simple[1] = new TF1("CoulombFFIntegrand_Z_imag_simple", 
		CoulombFFIntegrand_Z_simple, CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2, 4);
	f_CoulombFFIntegrand_Z_simple[0]->SetParameter(3, 1.0);
	f_CoulombFFIntegrand_Z_simple[1]->SetParameter(3,-1.0);
	
	//----------------------------------------//
	
	int nThetaBins = (int)csObj.m_CoulombFF.size();
	
	for(int itbin=0; itbin < nThetaBins; itbin++) {
		
		double locTheta = csObj.m_angularBins[itbin];
		
		// longitudinal component of momentum-transferred:
		double locDelta = 0.5 * pow(csObj.getMesonMass(),2.0) / csObj.getBeamEnergy();
		
		// transverse momentum-transferred:
		double locQ2 = 4.0 * csObj.getBeamEnergy() * csObj.getMesonMomentum(locTheta) * 
			pow(sin(0.5*locTheta*TMath::DegToRad()),2.0);
		double locQ  = sqrt(locQ2);
		
		// express q,Delta in fm^-1:
		locDelta /= PrimExCS::m_GeV2fm;
		locQ     /= PrimExCS::m_GeV2fm;
		
		f_CoulombFFIntegrand_B[0]->SetParameter(0, locQ);
		f_CoulombFFIntegrand_B[0]->SetParameter(1, locDelta);
		
		double ReFF = f_CoulombFFIntegrand_B[0]->Integral(0.0, CoulombFFIntegrationRange_b, 1.e-6);
		
		f_CoulombFFIntegrand_B[1]->SetParameter(0, locQ);
		f_CoulombFFIntegrand_B[1]->SetParameter(1, locDelta);
		
		double ImFF = f_CoulombFFIntegrand_B[1]->Integral(0.0, CoulombFFIntegrationRange_b, 1.e-6);
		
		cout << "theta = " << locTheta << "; Re(FF_em) = " << ReFF << ", Im(FF_em) = " << ImFF << endl;
		
		csObj.m_CoulombFF[itbin] = {ReFF, ImFF};
	}
	
	return;
}
