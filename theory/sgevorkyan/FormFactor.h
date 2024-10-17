#ifndef _PRIMEX_FORMFACTOR_
#define _PRIMEX_FORMFACTOR_

#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/GaussIntegrator.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

//--------------------------------//
// Struct to store model parameters:

struct modelParameters_t {
	double sigmaby2 =  0.0;
	double slopeAs  =  0.24;
	double slopeAp  =  0.24;
	double alphaIm  = -0.20;
	double shadowingParameter = 0.0;
	
	double targetA = 0.0;
	double targetRadius = 0.0;
	
	int densityModel = 0; // 0: Harmonic Oscillator
	                      // 1: Sum-of-Gaussians
	
	// Sum-of-Gaussians denisty parameterization:
	int nParametersSOG;
	double gammaSOG;
	vector<pair<double,double>> densityParametersSOG;
};

extern modelParameters_t modelParameters;

//--------------------------------//

int InitializeSOGDensityParameters(PrimExCS);

// global variables used in both form factor calculations:

static const double SQRT_PI = sqrt(TMath::Pi());
static const std::complex<double> ii = {0.0, 1.0};

//--------------------------------//
// Coulomb FF:

extern double integrationRadius;
extern double r_eff, coeff_ch, chargeDensityC1, chargeDensityC2;

void CalculateCoulombFF(PrimExCS &cs_obj);

extern TF1 *f_CoulombFFIntegrand_B[2];
extern double CoulombFFIntegrationRange_b;

extern TF1 *f_CoulombFFIntegrand_Z[2];
extern double CoulombFFIntegrationRange_z1, CoulombFFIntegrationRange_z2;

double CoulombFFIntegrate_Z(double isReal, double q, double Delta, double b);

long double CalcDiff_StruveBessel(long double t);

// ROOT Functions:

double CoulombFFIntegrand_B(double *x, double *par);
double CoulombFFIntegrand_Z(double *x, double *par);
double CoulombFFIntegrand_Z_simple(double *x, double *par);

//--------------------------------//
// Strong FF:

void CalculateStrongFF(PrimExCS &cs_obj);
std::complex<double> IntegrateInitialState(double b, double z, double q, double Delta, double Delta_rho);
void InitializeShadowingMatrix(double q, double Delta, double Delta_rho);

extern double StrongFFIntegrationRange_b;
extern double StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2;
extern double StrongFFIntegrationRange_s;

extern TF2  *f_StrongFFIntegrand_IS_2;
extern TH2F *h_shadowingMatrix[2];

// ROOT Functions:

double StrongFFIntegrand(double *x, double *par);
double StrongFFIntegrand_IS_1(double *x, double *par);
double StrongFFIntegrand_IS_2(double *x, double *par);

//--------------------------------//
// Used for absorption calculation:

std::complex<double> IntegrateAbsorption(double b, double z, double exponent=1.0);
double Interpolate2D(TH2F *h, double x, double y);

extern TH2F *h_absorptionMatrix;
void InitializeAbsorptionMatrix();

double Integrand_Gbz(double *x, double *par);
double GetGbzIntegral_HO(double b, double z);

double GetNuclearDensity(double r);
double density_HO(double r);
double density_SOG(double r);

double GetChargeDensityIntegral(double b, double z);

//--------------------------------//

#endif
