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
#include "TPaveStats.h"

typedef enum
{
	AFIX              = 0,
	SGEVORKYAN        = 1,
	MIXED_V1          = 2, // S.Gevorkyan calculation of Primakoff and Nuclear Coherent, A.Fix Calculation of Incoherent
	MIXED_V2          = 3, // A.Fix calculation of Primakoff and S.Gevorkyan everything else
	SGEVORKYAN_FERMI  = 4, // S.Gevorkyan calculation with Fermi Motion folded into incoherent part
	SGEVORKYAN_UPD_V0 = 5,
	SGEVORKYAN_UPD_V1 = 6,
	SGEVORKYAN_UPD_V2 = 7,
	SGEVORKYAN_UPD_V3 = 8,
	SGEVORKYAN_UPD_FERMI = 9
} ModelType;

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
			m_model = SGEVORKYAN;
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
		void FitAngularYield(double minFitRange=0.0, double maxFitRange=4.0);
		
		void SetModel(ModelType model) { m_model = model; }
		TString GetModel() {
			switch(m_model) {
				case AFIX:
					return "A. Fix";
				case SGEVORKYAN:
					return "S. Gevorkyan";
				case MIXED_V1:
					return "S. Gevorkyan (Coherent) + A. Fix (Incoherent)";
				case MIXED_V2:
					return "A. Fix (Primakoff) + S. Gevorkyan (Nuclear Coherent+Incoherent)";
				case SGEVORKYAN_FERMI:
					return "S. Gevorkyan (with Fermi Motion)";
				default:
					return "Unknown";
			}
		}
		
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
		vector<vector<TH1F*>> h_Theory;
		vector<TH1F*> h_PrimReal,   h_PrimImag;
		vector<TH1F*> h_StrongReal, h_StrongImag;
		
		double m_luminosity;
		
		vector<TString> m_components;
		int InitializeModelComponents();
		
		ModelType m_model;
		
		// Binning:
		
		double m_beamEnergyBinSize  =  0.050;
		double m_minBeamEnergy      =  9.000;
		double m_maxBeamEnergy      = 10.900;
		
		double m_reconAngleBinSize  =  0.060;
		double m_minReconAngle      =  0.000;
		double m_maxReconAngle      =  4.500;
		
		double m_thrownAngleBinSize =  0.010;
		double m_minThrownAngle     =  0.000;
		double m_maxThrownAngle     =  5.000;
		
		double GetCrossSection(int, int, double, double, double, double, double, double);
		double GetCrossSectionInterference(int, int, double, double, double, double);
		double GetCrossSectionInterferenceUpd(int, int, double, double, double, double);
		
		void InitializeFitFunction(TF1 **f1, TString funcName, int lineColor=kBlack);
		void InitializeDrawFunction(TF1 **f1, TString funcName, int lineColor=kBlack, int lineStyle=1, int lineWidth=2);
		void InitializeInterFunction(TF1 **f1, TString funcName, int lineColor=kBlack, int lineStyle=1, int lineWidth=2);
		void DrawFitResult(double min=0.0, double max=4.0);
};

#endif
