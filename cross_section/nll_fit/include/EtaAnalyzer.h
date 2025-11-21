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
			
			// Data:
			
			h_mggVsThetaFull(nullptr),     h_mggVsThetaEmpty(nullptr), 
			h_mggVsThetaFull_acc(nullptr), h_mggVsThetaEmpty_acc(nullptr), 
			
			// Simulated Lineshapes:
			
			h_incFraction(nullptr),
			h_etaLineshapeCoh(nullptr), 
			h_etaLineshapeQF(nullptr), 
			h_omegaLineshape(nullptr), 
			h_rhoLineshape(nullptr),
			
			h_eta1PionLineshape(nullptr), 
			h_eta2PionLineshape(nullptr), 
			h_eta3PionLineshape(nullptr), 
			h_hadronicBkgdLineshape(nullptr), 
			
			// Fraction of eta yield corresponding to each eta+X background:
			
			h_EtaPionBkgdFraction(nullptr), 
			h_EtaPionBkgdFraction_bggen(nullptr), h_EtaPionBkgdFraction_bggen_cut(nullptr), 
			h_HadronicBkgdFraction(nullptr), 
			h_HadronicBkgdFraction_bggen(nullptr), h_HadronicBkgdFraction_bggen_cut(nullptr), 
			
			// Other inputs needed for cross section
			
			h_fluxWeights(nullptr), 
			h_matrix(nullptr), 
			h_matrixFine(nullptr), 
			h_thrown(nullptr), 
			
			// Histograms to store results:
			
			h_Counts(nullptr), 
			h_EmptyCounts(nullptr), 
			
			h_EmptyEtaRatio(nullptr), 
			h_EmptyYield(nullptr), 
			h_OmegaYield(nullptr), 
			h_RhoYield(nullptr), 
			h_BkgdYield(nullptr), 
			h_EtaPionYield(nullptr), 
			h_HadronicBkgdYield(nullptr), 
			
			h_omega_mu_fit(nullptr),
			h_omega_sigma_fit(nullptr),
			h_omega_alpha_fit(nullptr),
			h_omega_n_fit(nullptr),
			
			// Objects used for drawing:
			
			cFit(nullptr), pFit(nullptr), pRes(nullptr), 
			cEmpty(nullptr), 
			cYield(nullptr), 
			cCrossSection(nullptr), 
			cAcceptance(nullptr), 
			cEmptyRatio(nullptr), 
			cEtaPionFraction(nullptr), 
			cHadronicBkgdFraction(nullptr), 
			cCounts(nullptr), 
			cOmegaFitPars(nullptr), 
			cBackgrounds(nullptr), 
			l0(nullptr), 
			lp(nullptr), 
			lm(nullptr), 
			lx1(nullptr), 
			lx2(nullptr)
		{
			// Run conditions:
			
			m_phase                =  3;
			
			m_analysisOption       =  0;
			m_vetoOption           =  6;
			
			m_mggHistName          = "mgg_const";
			m_matrixHistName       = "AngularMatrix";
			m_luminosity           =  0.0;
			m_emptyTargetFluxRatio =  1.0;
			
			m_IsMatrixLoaded       = false;
			
			// Binning defaults:
			
			m_binningSet         = false;
			
			m_rebinsMgg          =  2;
			m_mggBinSize         =  0.002;
			m_rebinsEmptyMgg     =  5;
			m_emptyMggBinSize    =  0.005;
			
			m_rebinsTheta        =  6;
			m_reconAngleBinSize  =  0.06;
			m_minReconAngle      =  0.00;
			m_maxReconAngle      =  4.50;
			
			m_thrownAngleBinSize =  0.01;
			m_minThrownAngle     =  0.00;
			m_maxThrownAngle     =  5.00;
			
			m_beamEnergyBinSize  =  0.05;
			m_minBeamEnergy      =  8.00;
			m_maxBeamEnergy      = 11.30;
			
			// Fitting Configuration:
			
			m_fitOption_signal     =  2;
			m_fitOption_bkgd       =  2;
			m_fitOption_poly       =  2;
			m_fitOption_omega      =  1;
			m_fitOption_rho        =  0;
			m_fitOption_etap       =  0;
			
			m_lineshapeOption      =  0;
			m_useRawMass           =  0;
			
			m_emptyFitOption_eta   =  2;
			m_emptyFitOption_omega =  1;
			m_emptyFitOption_fdc   =  1; 
			m_emptyFitOption_bkgd  =  3;
			m_emptyFitOption_poly  =  4;
			
			m_minFitRange          = 0.300;
			m_maxFitRange          = 0.950;
			m_minEmptyFitRange     = 0.300;
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
		void SetVetoOption(int);
		void SetMggHistName(TString);
		void SetMatrixHistName(TString);
		
		int GetPhase() { return m_phase; }
		int GetAnalysisOption() { return m_analysisOption; }
		int GetVetoOption() { return m_vetoOption; }
		
		double GetEmptyFluxRatio() { return m_emptyTargetFluxRatio; }
		
		// Binning:
		
		void SetRebinsMgg(int);
		void SetRebinsTheta(int);
		void SetRebinsEmptyMgg(int);
		void SetAngularRange(double min, double max) { m_minReconAngle = min; m_maxReconAngle = max; }
		void SetBeamEnergy(double, double);
		void SetBeamEnergy();
		
		double GetMggBinSize() { return m_mggBinSize; }
		double GetEmptyMggBinSize() { return m_emptyMggBinSize; }
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
		void SetFitOption_rho(int);
		void SetFitOption_etap(int);
		
		void SetEmptyFitOption_eta(int);
		void SetEmptyFitOption_omega(int);
		void SetEmptyFitOption_fdc(int);
		void SetEmptyFitOption_bkgd(int, int);
		void SetEmptyFitOption_bkgd(int);
		
		int  GetRawMassOption() { return m_useRawMass; }
		void SetRawMassOption(int opt) { m_useRawMass = opt; }
		
		int  GetLineshapeOption() { return m_lineshapeOption; }
		void SetLineshapeOption(int opt) { m_lineshapeOption = opt; }
		
		void SetFitRange(double, double);
		void SetEmptyFitRange(double, double);
		
		int  GetFitOption(int);
		int  GetEmptyFitOption(int);
		void GetFitRange(double&, double&);
		void GetEmptyFitRange(double&, double&);
		
		TString GetBkgdFitName();
		
		// For consistency checking:
		
		TString GetFitOptionStr(int);
		TString GetEmptyFitOptionStr();
		void DumpSettings();
		
		bool IsMatrixLoaded() { return m_IsMatrixLoaded; }
		
		// Load data (functions are defined in Data.cc):
		
		int LoadDataHistograms();
		
		int LoadLineshapes();
		int LoadIncoherentFraction();
		int LoadEtaLineshape_Coh();
		int LoadEtaLineshape_QF();
		int LoadEtaXLineshapes();
		int LoadOmegaLineshape();
		int LoadRhoLineshape();
		
		// Get photon flux (functinos are defined in Flux.cc):
		
		int    LoadLuminosity();
		int    LoadEmptyTargetFluxRatio();
		void   InitializeFluxHist();
		void IntegrateFluxHist(TH1F*, double, double, double&, double&);
		double GetLuminosity() { return m_luminosity; }
		TH1F* GetFluxWeights() { return h_fluxWeights; }
		
		// Angular matrices and acceptance (functions are defined in Acceptance.cc):
		
		int LoadAngularMatrix();
		int CalcAcceptance();
		TString GetMatrixFileName();
		TH3F* GetAngularMatrix() { return h_matrix; }
		TH3F* GetAngularMatrixFine() { return h_matrixFine; }
		TH2F* GetThrown() { return h_thrown; }
		
		// 
		
		void DrawInvariantMass(double minAngle=0.0, double maxAngle=5.0);
		void  FitInvariantMass(double minAngle=0.0, double maxAngle=5.0, int drawFitResult=0, int drawOmegaFit=0);
		
		void ExtractAngularYield(MggFitter&, int drawOption=0);
		
		void PlotAngularYield();
		void PlotCrossSection();
		void PlotEmptyEtaRatio();
		void PlotHadronicBkgdFraction();
		void PlotEtaPionFraction();
		void PlotBackgrounds();
		void PlotLineshapeShift();
		void PlotQFFraction();
		void PlotOmegaFitPars();
		void WriteROOTFile(TString fileName="yield.root");
		
		TH1F* GetAngularYield(int opt=0) { 
			if(opt==0)      return h_Yield;
			else if(opt==1) return h_YieldFit;
			else            return h_YieldInclusive;
		}
	
	private:
		
		// Data:
		
		TH2F *h_mggVsThetaFull,     *h_mggVsThetaEmpty;
		TH2F *h_mggVsThetaFull_acc, *h_mggVsThetaEmpty_acc;
		
		// Simulated Lineshapes:
		
		TH1F *h_incFraction;
		TH2F *h_etaLineshapeCoh,   *h_etaLineshapeQF;
		TH2F *h_omegaLineshape,    *h_rhoLineshape;
		TH2F *h_eta1PionLineshape, *h_eta2PionLineshape, *h_eta3PionLineshape;
		TH2F *h_hadronicBkgdLineshape;
		
		// Fraction of eta yield corresponding to each eta+X background:
		
		TH1F *h_EtaPionBkgdFraction;  // from fit result
		TH1F *h_EtaPionBkgdFraction_bggen, *h_EtaPionBkgdFraction_bggen_cut;  // from bggen simulation
		
		TH1F *h_HadronicBkgdFraction; // from fit result
		TH1F *h_HadronicBkgdFraction_bggen, *h_HadronicBkgdFraction_bggen_cut; // from bggen simulation
		
		// Other inputs needed for cross section
		
		TH1F *h_fluxWeights;
		
		TH3F *h_matrix, *h_matrixFine;
		TH2F *h_thrown;
		
		// Histograms to store results:
		
		TH1F *h_Counts, *h_EmptyCounts, *h_EmptyEtaRatio;
		TH1F *h_EmptyYield, *h_OmegaYield, *h_RhoYield, *h_BkgdYield, *h_EtaPionYield, *h_HadronicBkgdYield;
		
		// Objects needed for drawing:
		
		TCanvas *cFit;
		TPad *pFit, *pRes;
		void InitializeFitCanvas();
		
		TCanvas *cEmpty;
		void InitializeEmptyCanvas();
		
		TCanvas *cYield, *cCrossSection, *cAcceptance, *cEmptyRatio, *cEtaPionFraction, *cHadronicBkgdFraction;
		TCanvas *cCounts, *cBackgrounds;
		
		TCanvas *cOmegaFitPars;
		TH1F *h_omega_mu_fit, *h_omega_sigma_fit, *h_omega_alpha_fit, *h_omega_n_fit;
		
		TLine *l0, *lp, *lm, *lx1, *lx2;
		void DrawFitResult(TH1F*, TH1F*, TF1*, TF1*, TF1*, TF1*, TF1*, TF1*, TF1*, TF1*, TF1*, double minAngle, double maxAngle);
		
		//-----------------------------------------------------------//
		// Run-specific numbers:
		
		int m_phase, m_analysisOption, m_vetoOption;
		double m_luminosity, m_emptyTargetFluxRatio;
		TString m_mggHistName, m_matrixHistName;
		
		//-----------------------------------------------------------//
		// Variables defining bin sizes (defaults set in constructor):
		
		int    m_rebinsMgg;
		double m_mggBinSize;
		
		int    m_rebinsTheta;
		double m_reconAngleBinSize;
		double m_minReconAngle,  m_maxReconAngle;
		
		int    m_rebinsEmptyMgg;
		double m_emptyMggBinSize;
		
		double m_thrownAngleBinSize;
		double m_minThrownAngle, m_maxThrownAngle;
		
		double m_beamEnergyBinSize;
		double m_minBeamEnergy,  m_maxBeamEnergy;
		
		//-----------------------------------------------------------//
		
		bool m_binningSet;
		bool m_IsMatrixLoaded;
		
		// Fitting options:
		
		int m_fitOption_signal;
		int m_fitOption_bkgd, m_fitOption_poly;
		int m_fitOption_omega;
		int m_fitOption_rho;
		int m_fitOption_etap;
		
		int m_lineshapeOption;
		int m_useRawMass;
		
		int m_emptyFitOption_eta;
		int m_emptyFitOption_omega;
		int m_emptyFitOption_fdc;
		int m_emptyFitOption_bkgd, m_emptyFitOption_poly;
		
		double      m_minFitRange,      m_maxFitRange;
		double m_minEmptyFitRange, m_maxEmptyFitRange;
		
		// Vectors to store results:
		
		vector<pair<double,double>> m_angularBin;               // central angle and error of each bin
		vector<pair<double,double>> m_angularCounts;            // simple integration of histogram within mgg cut
		vector<pair<double,double>> m_angularCountsEmpty;       // same as above, but for empty target
		
		vector<pair<double,double>> m_angularYield;             // eta yield extracted as counts minus background fit functions
		vector<pair<double,double>> m_angularYieldFit;          // eta yield extracted as integrated signal fit function
		vector<pair<double,double>> m_angularYieldInc;
		
		vector<pair<double,double>> m_angularEtaPionBkgdFraction;
		vector<pair<double,double>> m_angularHadronicBkgdFraction;
		
		vector<pair<double,double>> m_angularYield_Empty;
		vector<pair<double,double>> m_angularYield_EtaPion;
		vector<pair<double,double>> m_angularYield_HadronicBkgd;
		vector<pair<double,double>> m_angularYield_Omega;
		vector<pair<double,double>> m_angularYield_Rho;
		vector<pair<double,double>> m_angularYield_Bkgd;
		
		vector<pair<double,double>> m_omegaFitMu, m_omegaFitSigma, m_omegaFitAlpha, m_omegaFitN;
		
		vector<pair<double,double>> m_fitResult_shift;
		vector<pair<double,double>> m_fitResult_bkgdShift;
		vector<pair<double,double>> m_fitResult_zqf;
		
		// Histograms to store results:
		
		TH1F *h_Yield,        *h_YieldFit,        *h_YieldInclusive;
		TH1F *h_CrossSection, *h_CrossSectionFit, *h_CrossSectionInclusive;
		TH1F *h_Acceptance;
		
		void IntegrateHistogram(TH1F*, double&, double&, double, double);
		void StyleYieldHistogram(TH1F *h1, int markerStyle=4, int markerColor=kBlack);
		void InitializeBinning();
};

#endif
