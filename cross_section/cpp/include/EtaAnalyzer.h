#ifndef _ETAGG_ANALYSIS_
#define _ETAGG_ANALYSIS_

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

#include "MggFitter.h"

class EtaAnalyzer {
	public:
		
		EtaAnalyzer() : 
			h_mggVsThetaFull(nullptr), 
			h_mggVsThetaEmpty(nullptr),  
			h_etaLineshape(nullptr), 
			h_omegaLineshape(nullptr), 
			h_etaPionLineshape(nullptr), 
			h_fdcOmegaLineshape(nullptr), 
			h_fluxWeights(nullptr), 
			h_matrix(nullptr), 
			h_matrixFine(nullptr), 
			cFit(nullptr), 
			cYield(nullptr), 
			cCrossSection(nullptr), 
			cAcceptance(nullptr), 
			l0(nullptr), 
			lp(nullptr), 
			lm(nullptr), 
			lx1(nullptr), 
			lx2(nullptr) 
		{
			// Run conditions:
			
			m_phase                =  1;
			m_analysisOption       =  1;
			m_mggHistName          = "mgg_const";
			m_matrixHistName       = "AngularMatrix";
			m_luminosity           =  0.0;
			m_emptyTargetFluxRatio =  1.0;
			
			// Binning:
			
			m_rebinsMgg          =  4;
			m_mggBinSize         =  0.04;
			
			m_rebinsTheta        =  6;
			m_reconAngleBinSize  =  0.06;
			m_minReconAngle      =  0.00;
			m_maxReconAngle      =  4.50;
			
			m_thrownAngleBinSize =  0.01;
			m_minThrownAngle     =  0.00;
			m_maxThrownAngle     =  5.00;
			
			m_beamEnergyBinSize  =  0.05;
			m_minBeamEnergy      =  9.00;
			m_maxBeamEnergy      = 10.90;
			
			// Fitting Configuration:
			
			m_subtractEmpty        = 1;
			m_fitOption_empty      = 0;
			
			m_fitOption_signal     = 1;
			m_fitOption_bkgd       = 3;
			m_fitOption_poly       = 1;
			m_fitOption_omega      = 2;
			m_fitOption_etap       = 1;
			
			m_emptyFitOption_eta   = 0;
			m_emptyFitOption_omega = 0;
			m_emptyFitOption_fdc   = 0; 
			m_emptyFitOption_bkgd  = 2;
			m_emptyFitOption_poly  = 3;
			
			m_minFitRange          = 0.300;
			m_maxFitRange          = 1.100;
			m_minEmptyFitRange     = 0.405;
			m_maxEmptyFitRange     = 0.950;
		};
		
		~EtaAnalyzer(){};
		
		// Constants:
		
		static constexpr double m_massEta        = 0.54786;    // GeV/c2
		static constexpr double m_massOmega      = 0.78266;    // GeV/c2
		static constexpr double m_massEtap       = 0.95778;    // GeV/c2
		static constexpr double m_branchingRatio = 0.3936;     //
		static constexpr double m_targetDensity  = 0.1217;     // g/cm3
		static constexpr double m_targetLength   = 29.5;       // cm
		static constexpr double m_targetMass     = 4.002602;   // g/mol
		static constexpr double m_avogdroNum     = 6.02214e23; // atoms/mol
		
		// Run-conditions:
		
		void SetPhase(int);
		void SetAnalysisOption(int);
		void SetMggHistName(TString);
		void SetMatrixHistName(TString);
		
		int GetPhase() { return m_phase; }
		
		// Binning:
		
		void SetRebinsMgg(int);
		void SetRebinsTheta(int);
		void SetBeamEnergy(double, double);
		
		double GetMggBinSize() { return m_mggBinSize; }
		void  GetBeamEnergyBinning(double&, double&, double&);
		void  GetReconAngleBinning(double&, double&, double&);
		void GetThrownAngleBinning(double&, double&, double&);
		
		// Fitting options:
		
		void SetSubtractEmptyTarget(int);
		void SetFitEmptyTarget(int);
		
		void SetFitOption_signal(int);
		void SetFitOption_bkgd(int, int);
		void SetFitOption_bkgd(int);
		void SetFitOption_omega(int);
		void SetFitOption_etap(int);
		
		void SetEmptyFitOption_eta(int);
		void SetEmptyFitOption_omega(int);
		void SetEmptyFitOption_fdc(int);
		void SetEmptyFitOption_bkgd(int, int);
		void SetEmptyFitOption_bkgd(int);
		
		void SetFitRange(double, double);
		void SetEmptyFitRange(double, double);
		
		int  GetFitOption(int);
		int  GetEmptyFitOption(int);
		int  GetEmptySubtractOption() { return m_subtractEmpty; }
		void GetFitRange(double&, double&);
		void GetEmptyFitRange(double&, double&);
		
		TString GetBkgdFitName();
		
		// For consistency checking:
		
		TString GetFitOptionStr(int);
		void DumpSettings();
		
		// Load data (functions are defined in Data.cc):
		
		int LoadDataHistograms();
		
		int LoadLineshapes();
		int LoadEtaLineshape();
		int LoadEtaPionLineshape();
		int LoadOmegaLineshape();
		int LoadFDCOmegaLineshape();
		
		// Get photon flux (functinos are defined in Flux.cc):
		
		int    LoadLuminosity();
		int    LoadEmptyTargetFluxRatio();
		void   InitializeFluxHist();
		double IntegrateFluxHist(TH1F*, double, double);
		double GetLuminosity() { return m_luminosity; }
		TH1F* GetFluxWeights() { return h_fluxWeights; }
		
		// Angular matrices and acceptance (functions are defined in Acceptance.cc):
		
		int LoadAngularMatrix();
		int CalcAcceptance();
		TH3F* GetAngularMatrix() { return h_matrix; }
		TH3F* GetAngularMatrixFine() { return h_matrixFine; }
		
		// 
		
		void DrawInvariantMass(double minAngle=0.0, double maxAngle=5.0);
		void  FitInvariantMass(double minAngle=0.0, double maxAngle=5.0, int drawFitResult=0, int drawOmegaFit=0);
		
		void ExtractAngularYield(int drawOption=0);
		void    PlotAngularYield();
		
		void PlotCrossSection();
		
		TH1F* GetAngularYield() { return h_YieldFit; }
	
	private:
		
		// Data:
		TH2F *h_mggVsThetaFull, *h_mggVsThetaEmpty;
		
		// MC Lineshapes:
		TH2F *h_etaLineshape, *h_omegaLineshape, *h_etaPionLineshape, *h_fdcOmegaLineshape;
		
		// Photon flux weights:
		TH1F *h_fluxWeights;
		
		// Angular Matrix:
		TH3F *h_matrix, *h_matrixFine;
		
		//-----------------------------------------------------------//
		// Run-specific numbers:
		
		int m_phase, m_analysisOption;
		double m_luminosity, m_emptyTargetFluxRatio;
		TString m_mggHistName, m_matrixHistName;
		
		//-----------------------------------------------------------//
		// Variables defining bin sizes (defaults set in constructor):
		
		int    m_rebinsMgg;
		double m_mggBinSize;
		
		int    m_rebinsTheta;
		double m_reconAngleBinSize;
		double m_minReconAngle,  m_maxReconAngle;
		
		double m_thrownAngleBinSize;
		double m_minThrownAngle, m_maxThrownAngle;
		
		double m_beamEnergyBinSize;
		double m_minBeamEnergy,  m_maxBeamEnergy;
		
		//-----------------------------------------------------------//
		
		bool m_binningSet = false;
		
		// Fitting options:
		
		int m_subtractEmpty, m_fitOption_empty;
		
		int m_fitOption_signal;
		int m_fitOption_bkgd, m_fitOption_poly;
		int m_fitOption_omega;
		int m_fitOption_etap;
		
		int m_emptyFitOption_eta;
		int m_emptyFitOption_omega;
		int m_emptyFitOption_fdc; 
		int m_emptyFitOption_bkgd, m_emptyFitOption_poly;
		
		double      m_minFitRange,      m_maxFitRange;
		double m_minEmptyFitRange, m_maxEmptyFitRange;
		
		// Vectors to store results:
		
		vector<pair<double,double>> m_angularBin;
		vector<pair<double,double>> m_angularYield, m_angularYieldFit, m_angularYieldEmpty;
		vector<pair<double,double>> m_angularSBR;
		
		// Histograms to store results:
		
		TH1F *h_Yield,        *h_YieldFit;
		TH1F *h_CrossSection, *h_CrossSectionFit;
		TH1F *h_Acceptance;
		
		void InitializeBinning();
		
		// Objects needed for drawing:
		
		TCanvas *cFit, *cYield, *cCrossSection, *cAcceptance;
		TPad *pFit, *pRes;
		void InitializeFitCanvas();
		
		TCanvas *cEmpty = NULL;
		void InitializeEmptyCanvas();
		
		TLine *l0, *lp, *lm, *lx1, *lx2;
		void DrawFitResult(TH1F*, TH1F*, TF1*, TF1*, TF1*, TF1*, TF1*, double minAngle, double maxAngle);
};

#endif
