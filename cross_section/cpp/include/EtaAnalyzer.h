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

#include "MggFitter.h"

class EtaAnalyzer {
	public:
		EtaAnalyzer(int);
		~EtaAnalyzer(){};
		
		// Constants:
		
		static constexpr double m_massEta        = 0.54786;
		static constexpr double m_massOmega      = 0.78266;
		static constexpr double m_massEtap       = 0.95778;
		static constexpr double m_branchingRatio = 0.3936;
		
		// defaults for beam energy range:
		double minBeamEnergy     =  9.00;
		double maxBeamEnergy     = 10.90;
		double beamEnergyBinSize =  0.05;
		
		void SetRebinsMgg(int);
		void SetRebinsTheta(int);
		
		// Fitting options:
		void SetFitOption_signal(int);
		void SetFitOption_bkgd(int, int);
		void SetFitOption_bkgd(int);
		void SetFitOption_omega(int);
		void SetFitOption_etap(int);
		void SetFitRange(double, double);
		
		int  GetFitOption(int);
		void GetFitRange(double&, double&);
		double GetMggBinSize() { return m_mggBinSize; }
		
		// For consistency checking:
		TString GetFitOptionStr(int);
		void DumpSettings();
		
		// Read in Data:
		int GetDataHistograms(int anaOption, TString histName="");
		
		int LoadLineshapes(int anaOption, TString histName="");
		int LoadEtaLineshape(int anaOption, TString histName="");
		int LoadOmegaLineshape(int anaOption, TString histName="");
		
		void DrawInvariantMass(double minAngle=0.0, double maxAngle=5.0);
		void  FitInvariantMass(double minAngle=0.0, double maxAngle=5.0, int drawFitResult=0, int drawOmegaFit=0);
		
		double IntegrateFluxHist(TH1F *hFlux);
		void CalcLuminosity();
		void CalcEmptyTargetFluxRatio();
		
		void ExtractAngularYield(int drawOption=0);
		void    PlotAngularYield();
		
		void PlotCrossSection();
		int GetAcceptance(TString matrixName="AngularMatrix");
	
	private:
		
		// Run-specific numbers:
		
		int m_phase; // supported options: 1-3
		double m_luminosity = 0., m_emptyTargetFluxRatio = 1.;
		
		// Variables defining bin sizes (defaults set in constructor):
		
		int m_rebinsMgg, m_rebinsTheta;
		double m_mggBinSize, m_thetaBinSize;
		
		vector<double> m_fluxWeights;
		
		bool m_binningSet = false;
		
		// Fitting options:
		
		int m_fitOption_signal = 1;
		int m_fitOption_bkgd   = 3;
		int m_fitOption_poly   = 1;
		int m_fitOption_omega  = 2;
		int m_fitOption_etap   = 1;
		
		double m_minFitRange=0.30, m_maxFitRange=1.10;
		
		TH2F *h_etaLineshape;
		TH2F *h_omegaLineshape;
		
		// Data:
		
		TH2F *h_mggVsThetaFull, *h_mggVsThetaEmpty;
		
		vector<pair<double,double>> m_angularBin;
		vector<pair<double,double>> m_angularYield, m_angularYieldFit, m_angularYieldEmpty;
		vector<pair<double,double>> m_angularSBR;
		
		TH1F *h_Yield,        *h_YieldFit;
		TH1F *h_CrossSection, *h_CrossSectionFit;
		TH1F *h_Acceptance;
		
		void InitializeBinning();
		
		TCanvas *cFit, *cYield, *cCrossSection, *cAcceptance;
		TPad *pFit, *pRes;
		void InitializeFitCanvas();
		
		TLine *l0, *lp, *lm, *lx1, *lx2;
		void DrawFitResult(TH1F*, TH1F*, TF1*, TF1*, TF1*);
};

#endif
