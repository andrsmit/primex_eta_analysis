#ifndef _YIELD_FITTER_
#define _YIELD_FITTER_

using namespace std;

#include <stdio.h>
#include <vector>
#include <iostream>

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"

class YieldFitter {
	public:
		YieldFitter() : 
			h_yield(nullptr), 
			f_yield(nullptr), 
			h_matrix(nullptr), 
			h_matrixFine(nullptr), 
			h_fluxWeights(nullptr), 
			c_fit(nullptr)
		{
			m_luminosity  = 0.0;
		}
		
		// Fit Functions:
		double YieldFitFunction(double *x, double *par);
		double YieldDrawFunction(double *x, double *par);
		double YieldDrawFunctionInterference(double *x, double *par);
		
		void SetBeamEnergyBinning(double binSize, double min, double max) {
			m_beamEnergyBinSize  = binSize;
			m_minBeamEnergy      = min;
			m_maxBeamEnergy      = max;
			return;
		}
		void SetReconAngleBinning(double binSize, double min, double max) {
			m_reconAngleBinSize  = binSize;
			m_minReconAngle      = min;
			m_maxReconAngle      = max;
			return;
		}
		void SetThrownAngleBinning(double binSize, double min, double max) {
			m_thrownAngleBinSize = binSize;
			m_minThrownAngle     = min;
			m_maxThrownAngle     = max;
			return;
		}
		void SetLuminosity(double lumi) { m_luminosity = lumi; return; }
		void SetAngularMatrix(TH3F *h3) { h_matrix = h3; return; }
		void SetAngularMatrixFine(TH3F *h3) { h_matrixFine = h3; return; }
		void SetFluxWeights(TH1F *h1) { h_fluxWeights = h1; return; }
		void SetYield(TH1F *h1) { h_yield = h1; return; }
		
		int LoadTheoryHists();
		void FitAngularYield(int drawOption=1);
		
	private:
		// Angular Yield histogram to be fit:
		TH1F *h_yield;
		
		// Fit Function:
		TF1 *f_yield;
		
		// Angular Matrix:
		TH3F *h_matrix;
		
		// more finely binned matrix for drawing purposes:
		TH3F *h_matrixFine;
		
		// Histogram to store the fraction of photon flux in each energy bin:
		TH1F *h_fluxWeights;
		
		// Canvas for drawing fit results:
		TCanvas *c_fit;
		
		// Theory Histograms:
		vector<TH1F*> h_TheoryPrim;
		vector<TH1F*> h_TheoryCoh;
		vector<TH1F*> h_TheoryQFP;
		vector<TH1F*> h_TheoryQFN;
		
		double m_luminosity;
		
		// Binning:
		
		double m_beamEnergyBinSize  =  0.05;
		double m_minBeamEnergy      =  9.00;
		double m_maxBeamEnergy      = 10.90;
		
		double m_reconAngleBinSize  =  0.06;
		double m_minReconAngle      =  0.00;
		double m_maxReconAngle      =  4.50;
		
		double m_thrownAngleBinSize =  0.01;
		double m_minThrownAngle     =  0.00;
		double m_maxThrownAngle     =  5.00;
		
		double GetCrossSection(int, int, double, double, double, double, double);
		double GetCrossSectionInterference(int, int, double, double, double);
		
		void InitializeFitFunction(TF1 **f1, TString funcName);
		void DrawFitResult();
};

#endif
