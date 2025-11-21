#ifndef _YIELD_FITTER_
#define _YIELD_FITTER_

using namespace std;

#include "CrossSection.h"

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
#include "TSpline.h"
#include "TAxis.h"
#include "TString.h"

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/ParameterSettings.h"
#include "Math/Functor.h"
#include "TVirtualFitter.h"

class YieldFitter {
	public:
		YieldFitter() : 
			h_matrixFull(nullptr), 
			h_fluxWeightsFull(nullptr),
			h_PrimReal(nullptr),
			h_PrimImag(nullptr),
			h_StrongReal(nullptr),
			h_StrongImag(nullptr)
		{
			m_luminosity  = 0.0;
			m_model = UNKNOWN;
			
			h_dNdTheta.clear();
			h_dNdOmega.clear();
			f_dNdTheta.clear();
			m_energyBins.clear();
			c_dNdTheta.clear();
			h_matrices.clear();
			h_fluxWeights.clear();
		}
		
		//----------------------------------------//
		// Angular Yield distributions to be fit:
		
		vector<TH1F*> h_dNdTheta, h_dNdOmega;
		
		// Fit Functions:
		
		vector<TF1*> f_dNdTheta;
		
		// Energy binning of yield histograms:
		
		vector<pair<double,double>> m_energyBins;
		
		//----------------------------------------//
		
		double minFitRange = 0.0;
		double maxFitRange = 4.0;
		
		// Setter functions for internal parameters:
		
		void SetBeamEnergyBinning(double binSize, double min, double max) {
			m_beamEnergyBinSize  = binSize;
			m_minBeamEnergy      = min;
			m_maxBeamEnergy      = max;
		}
		void SetReconAngleBinning(double binSize, double min, double max) {
			m_reconAngleBinSize  = binSize;
			m_minReconAngle      = min;
			m_maxReconAngle      = max;
		}
		void SetThrownAngleBinning(double binSize, double min, double max) {
			m_thrownAngleBinSize = binSize;
			m_minThrownAngle     = min;
			m_maxThrownAngle     = max;
		}
		void SetLuminosity(double lumi) { m_luminosity = lumi; }
		void SetAngularMatrixFull(TH3F *h3) { h_matrixFull = h3; }
		void SetThrown(TH2F *h2) { h_thrown = h2; }
		void SetFluxWeights(TH1F *h1) { h_fluxWeightsFull = h1; }
		
		void SetModel(ModelType model) { m_model = model; }
		ModelType GetModel() { return m_model; }
		
		void SetModel_apVer(int ver) { m_apVer = ver; }
		void SetModel_sigmaVer(int ver) { m_sigmaVer = ver; }
		void SetModel_strongRadiusStr(TString str) { m_strongRadiusStr = str; }
		
		TString GetModelString() {
			switch(m_model) {
				case AFIX:
					return "A. Fix";
				case SGEVORKYAN:
					return "S. Gevorkyan (original)";
				case MIXED_V1:
					return "S. Gevorkyan (Coherent) + A. Fix (Incoherent)";
				case MIXED_V2:
					return "A. Fix (Primakoff) + S. Gevorkyan (Nuclear Coherent+Incoherent)";
				case SGEVORKYAN_FERMI:
					return "S. Gevorkyan (with Fermi Motion)";
				case SGEVORKYAN_UPD_V0:
					return "S. Gevorkyan (no FSI, no ISI)";
				case SGEVORKYAN_UPD_V1:
					return "S. Gevorkyan (with FSI, no ISI)";
				case SGEVORKYAN_UPD_V2:
					return "S. Gevorkyan (with FSI, with ISI)";
				case SGEVORKYAN_UPD_V3:
					return "S. Gevorkyan (no FSI, with ISI)";
				case SGEVORKYAN_SIGMA_VAR:
					return Form("S. Gevorkyan (sigma index: %d)", m_sigmaVer);
				case SGEVORKYAN_AP_VAR:
					return Form("S. Gevorkyan (ap index: %d)", m_apVer);
				case SGEVORKYAN_STRONG_RADIUS_VAR:
					return Form("S. Gevorkyan (strong radius index: %s)", m_strongRadiusStr.Data());
				default:
					return "Unknown";
			}
		}
		
		// Two Options for initializing yield histograms:
		
		int SetYield(TH1F *h1)
		{
			h_dNdTheta.clear();
			h_dNdTheta.push_back(h1);
			
			m_energyBins.clear();
			m_energyBins.push_back({m_minBeamEnergy, m_maxBeamEnergy});
			return 0;
		}
		
		int SetYield(vector<TH1F*> yieldHists, vector<pair<double,double>> energyBins)
		{
			unsigned int nHists = yieldHists.size();
			unsigned int nBins  = energyBins.size();
			if(nHists!=nBins) {
				printf("Inconsistent number of yield histograms and energy bins provided to YieldFitter\n");
				return 1;
			}
			
			double locMinEnergy = m_minBeamEnergy;
			double locMaxEnergy = m_maxBeamEnergy;
			
			h_dNdTheta.clear(); m_energyBins.clear();
			for(unsigned int ih=0; ih<nHists; ih++)
			{
				h_dNdTheta.push_back(yieldHists[ih]);
				m_energyBins.push_back(energyBins[ih]);
				
				if(ih==0) {
					locMinEnergy = energyBins[ih].first;
					locMaxEnergy = energyBins[ih].second;
				}
				else {
					if(energyBins[ih].first<locMaxEnergy) {
						printf("Conflicting energy bins provided to YieldFitter.\n");
						exit(1);
					}
					locMaxEnergy = energyBins[ih].second;
				}
			}
			m_minBeamEnergy = locMinEnergy;
			m_maxBeamEnergy = locMaxEnergy;
			return 0;
		}
		
		// Fit/Draw Functions:
		
		double GetExpectedYield(double, double, int, double *par, int isolateInter=0);
		
		double YieldFitFunction(double *x, double *par);
		double YieldFitFunction_Interference(double *x, double *par);
		
		int LoadTheoryHists();
		void FitAngularYield(double minFitRange=0.0, double maxFitRange=4.0, TString outputFileName="");
		
		struct CombinedChi2;
		
	private:
		
		vector<TString> m_parameterList;
		vector<vector<int>> m_parIndices;
		
		// Angular Matrices:
		TH3F *h_matrixFull;
		TH2F *h_thrown;
		vector<TH3F*> h_matrices;
		
		// Histogram to store the fraction of photon flux in each energy bin:
		TH1F *h_fluxWeightsFull;
		vector<TH1F*> h_fluxWeights;
		
		// Canvas for drawing fit results:
		vector<TCanvas*> c_dNdTheta;
		
		// Theory Histograms:
		vector<TH2F*> h_Theory;
		TH2F *h_PrimReal, *h_StrongReal;
		TH2F *h_PrimImag, *h_StrongImag;
		
		TH1F *h_TheoryTulio;
		TSpline3 *spline;
		TF1 *f_TheoryTulio;
		
		double m_luminosity;
		
		vector<TString> m_components;
		int InitializeModelComponents();
		
		void InitializeMatrices();
		void InitializeFluxWeights();
		
		void InitializeParameterArrays();
		void InitializeFitParameters(ROOT::Fit::Fitter&);
		void DumpFitParameters(ROOT::Fit::Fitter &fitter);
		
		void UpdateFitFunctions(ROOT::Fit::FitResult);
		
		void WriteOutputASCII(int, ROOT::Fit::FitResult, TString fileName="");
		void DrawFitResult(double min=0.0, double max=4.0, TString fileName="");
		
		ModelType m_model;
		
		int m_sigmaVer = 3, m_apVer = 3;
		TString m_strongRadiusStr = "0";
		
		// Binning:
		
		double m_beamEnergyBinSize  =  0.050;
		double m_minBeamEnergy      =  8.000;
		double m_maxBeamEnergy      = 11.300;
		
		double m_reconAngleBinSize  =  0.060;
		double m_minReconAngle      =  0.000;
		double m_maxReconAngle      =  4.500;
		
		double m_thrownAngleBinSize =  0.010;
		double m_minThrownAngle     =  0.000;
		double m_maxThrownAngle     =  5.000;
		
		double GetCrossSection(double, int, double, double, double, double, double);
		double GetCrossSectionInterference(double, int, double, double, double);
		
		void InitializeFitFunction(TF1 **f1, TString funcName, int lineColor=kBlack, int lineStyle=1, int lineWidth=2);
		void InitializeFitFunction_Interference(TF1 **f1, TString funcName, int lineColor=kMagenta, int lineStyle=4, int lineWidth=2);
};

#endif
