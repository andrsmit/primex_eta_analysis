#ifndef _ETAGG_CS_
#define _ETAGG_CS_

typedef enum
{
	UNKNOWN              = 0,
	AFIX                 = 1,
	SGEVORKYAN           = 2,
	MIXED_V1             = 3, // S.Gevorkyan calculation of Primakoff and Nuclear Coherent, A.Fix Calculation of Incoherent
	MIXED_V2             = 4, // A.Fix calculation of Primakoff and S.Gevorkyan everything else
	SGEVORKYAN_FERMI     = 5, // S.Gevorkyan calculation with Fermi Motion folded into incoherent part
	SGEVORKYAN_UPD_V0    = 6,
	SGEVORKYAN_UPD_V1    = 7,
	SGEVORKYAN_UPD_V2    = 8,
	SGEVORKYAN_UPD_V3    = 9,
	SGEVORKYAN_UPD_FERMI = 10,
	SGEVORKYAN_SIGMA_VAR = 11,
	SGEVORKYAN_AP_VAR    = 12,
	SGEVORKYAN_STRONG_RADIUS_VAR = 13
} ModelType;

#include "EtaAnalyzer.h"
#include "MggFitter.h"
#include "YieldFitter.h"
#include "MyReadConfig.h"
#include "TApplication.h"
#include "TGraphErrors.h"

struct configSettings_t {
	int primexPhase    = -1;
	int analysisOption =  1;
};

ModelType StringToModelType(std::string str);

TString GetOutputFileName(EtaAnalyzer);

int LoadConfigSettings(EtaAnalyzer&, TString, int);
int LoadModelType(YieldFitter&, TString);

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
