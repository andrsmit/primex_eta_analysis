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

// PDF functions for Lineshape fitting:

double DoubleGausPDF(double*,double*);

double CrystalBallPDF(double*,double*);
double CrystalBallPDF_flip(double*,double*);

double DoubleCrystalBallPDF(double*,double*);
double DoubleCrystalBallPDF_flip(double*,double*);
double DoubleCrystalBallPDF_oneflip(double*,double*);
double DoubleCrystalBallPlusGausPDF(double*,double*);

double NormGaus(double x, double mu, double sigma);
double NormCrystalBall(double x, double mu, double sigma, double alpha, double n, int doFlip=0);

// Function to mgg initialize fitter object:
void InitializeMggFitter(MggFitter &fitter, EtaAnalyzer *anaObj, double angle, double angleWidth=0.0);

// Function to yield initialize fitter object:
int InitializeYieldFitter(YieldFitter &fitter, EtaAnalyzer &anaObj);

// Main function of executable:

int main(int argc, char **argv);

#endif
