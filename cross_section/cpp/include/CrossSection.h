#ifndef _ETAGG_CS_
#define _ETAGG_CS_

#include "EtaAnalyzer.h"
#include "MggFitter.h"
#include "YieldFitter.h"
#include "MyReadConfig.h"
#include "TApplication.h"

struct configSettings_t {
	int primexPhase    = -1;
	int analysisOption =  1;
};

int LoadConfigSettings(EtaAnalyzer&, TString, int);
void printUsage(configSettings_t, int);

// For styling:

void styleMggHistogram(TH1F *h1, int lineColor=kBlack, int markerStyle=20);
void styleCanvas(TCanvas *c1);
void styleCanvas(TPad *p1);

// Chebshev polynomials:
double Chebyshev0(double, double);
double Chebyshev1(double, double, double);
double Chebyshev2(double, double, double, double);
double Chebyshev3(double, double, double, double, double);
double Chebyshev4(double, double, double, double, double, double);
double Chebyshev5(double, double, double, double, double, double, double);

// Crystal Ball functions for Omega Lineshape fitting:
double CrystalBall(double *x, double *par);
double CrystalBall2(double *x, double *par);
double CrystalBall_flip(double *x, double *par);
double CrystalBall2_flip(double *x, double *par);
double CrystalBall3_flip(double *x, double *par);

double DoubleGaus(double *x, double *par);

// Function to mgg initialize fitter object:
void InitializeMggFitter(MggFitter &fitter, EtaAnalyzer *anaObj, double angle);

// Function to yield initialize fitter object:
int InitializeYieldFitter(YieldFitter &fitter, EtaAnalyzer anaObj);

// Main function of executable:

int main(int argc, char **argv);

#endif
