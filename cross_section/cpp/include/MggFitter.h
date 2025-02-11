#ifndef _MGG_FITTER_
#define _MGG_FITTER_

using namespace std;

#include <stdio.h>
#include <vector>
#include <iostream>

#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"

#include "EtaAnalyzer.h"

class MggFitter {
	public:
		MggFitter() 
			: h_data(nullptr), h_empty(nullptr), f_fit(nullptr), f_empty(nullptr), 
				h_etaLineshape(nullptr), f_etaLineshape(nullptr),
				h_omegaLineshape(nullptr), f_omegaLineshape(nullptr),
				h_etaPionLineshape(nullptr), f_etaPionLineshape(nullptr), 
				h_fdcOmegaLineshape(nullptr), f_fdcOmegaLineshape(nullptr)
		{
			f_chebyshev = new TF1("chebyshev", "cheb5", 0.0, 1.0);
		}
		
		double angle   = 0.0;
		double binSize = 0.001;
		
		// Fitting options:
		
		int fitOption_signal = 1;
		int fitOption_bkgd   = 3;
		int fitOption_poly   = 1;
		int fitOption_omega  = 2;
		int fitOption_etap   = 1;
		
		// Fitting options for empty target background:
		
		int fitOption_empty = 0;
		
		int emptyFitOption_eta   = 0;
		int emptyFitOption_omega = 0;
		int emptyFitOption_fdc   = 0;
		int emptyFitOption_bkgd  = 2;
		int emptyFitOption_poly  = 3; 
		
		double      minFitRange=0.300,      maxFitRange=1.100;
		double minEmptyFitRange=0.405, maxEmptyFitRange=0.950;
		
		double minMggCut = 0.50;
		double maxMggCut = 0.60;
		
		// Set full and empty target histograms:
		void SetData(TH1F* h1) { h_data = h1; return; }
		void SetEmpty(TH1F *h1, double ratio=1.0, double ratioErr = 0.05) { 
			h_empty = h1;
			m_emptyRatio    = ratio;
			m_emptyRatioErr = ratioErr;
			return;
		}
		
		// Set/Fit lineshape of eta signal:
		void SetEtaLineshape(TH1F *h1, int drawOption=0) { 
			h_etaLineshape = h1;
			h_etaLineshape->Scale(1.0/h_etaLineshape->Integral());
			FitEtaLineshape(drawOption);
			return;
		}
		void FitEtaLineshape(int drawOption=0);
		
		// Set/fit lineshape of eta+pion background:
		void SetEtaPionLineshape(TH1F *h1, int drawOption=0) {
			h_etaPionLineshape = h1;
			h_etaPionLineshape->Scale(1.0/h_etaPionLineshape->Integral());
			FitEtaPionLineshape(drawOption);
			return;
		}
		void FitEtaPionLineshape(int drawOption=0);
		
		// Set/fit lineshape of omega background:
		void SetOmegaLineshape(TH1F *h1, int drawOption=0) { 
			h_omegaLineshape = h1;
			h_omegaLineshape->Scale(1.0/h_omegaLineshape->Integral());
			FitOmegaLineshape(drawOption);
			return;
		}
		void FitOmegaLineshape(int drawOption=0);
		
		// Set/fit lineshape of omegas form 1st FDC package:
		void SetFDCOmegaLineshape(TH1F *h1, int drawOption=0) { 
			h_fdcOmegaLineshape = h1;
			h_fdcOmegaLineshape->Scale(1.0/h_fdcOmegaLineshape->Integral());
			FitFDCOmegaLineshape(drawOption);
			return;
		}
		void FitFDCOmegaLineshape(int drawOption=0);
		
		void FitData();
		void FitEmpty();
		void FitDataWithEmpty();
		void DrawFitResults(TCanvas&);
		
		TF1* GetFitFunction() { 
			//ZeroSignalPars(f_fit);
			f_fit->SetLineColor(kGreen);
			f_fit->SetNpx(1000);
			return f_fit;
		}
		TF1* GetSignalFunction() {
			TF1 *f1;
			InitializeFitFunction(&f1,"signalFit");
			f1->SetParameters(f_fit->GetParameters());
			ZeroBkgdPars(f1);
			f1->SetLineColor(kCyan+2);
			f1->SetNpx(1000);
			return f1;
		}
		TF1* GetBkgdFunction() {
			TF1 *f1;
			InitializeFitFunction(&f1,"bkgdFit");
			f1->SetParameters(f_fit->GetParameters());
			ZeroSignalPars(f1);
			f1->SetParameter("N_{#omega}", 0.0);
			f1->SetParameter("N_{empty}", 0.0);
			f1->SetLineColor(kMagenta);
			f1->SetLineStyle(4);
			f1->SetNpx(1000);
			return f1;
		}
		TF1* GetEtaPionFunction() {
			TF1 *f1;
			InitializeFitFunction(&f1,"etaPionFit");
			f1->SetParameters(f_fit->GetParameters());
			f1->SetParameter("N_{#eta}", 0.0);
			ZeroBkgdPars(f1);
			f1->SetParameter("N_{#eta#pi}", f_fit->GetParameter("N_{#eta#pi}"));
			f1->SetLineColor(kCyan);
			f1->SetNpx(1000);
			return f1;
		}
		TF1* GetEmptyFitFunction() {
			TF1 *f1;
			InitializeEmptyFitFunction(&f1,"emptyFitFunction");
			f1->SetParameters(f_empty->GetParameters());
			f1->SetLineColor(kRed);
			f1->SetLineStyle(4);
			f1->SetLineWidth(2);
			f1->SetNpx(1000);
			/*
			f1->SetParameter("N_{#eta}",   0.0);//m_nEmptyEtaPar);
			f1->SetParameter("N_{#omega}", 0.0);//m_nEmptyOmegaPar);
			*/
			return f1;
		}
		double GetEmptyEtaFitPar() { return m_nEmptyEtaPar; }
		double GetEmptyOmegaFitPar() { return m_nEmptyOmegaPar; }
		
		void FillPull(TH1F*);
		void FillEmptyPull(TH1F*);
		
		void GetYield(double&,double&,int useSignalPars=0);
		
		double IntegrateGaussian(double A, double mu, double sigma, double x1, double x2);
		
		void DumpFitParameters();
		void DumpEmptyFitParameters();
		
		
		// Fit Functions:
		double      MggFitFunction(double *x, double *par);
		double EmptyMggFitFunction(double *x, double *par);
		
	
	private:
		// Data histogram to be fit:
		TH1F *h_data;
		
		// Empty target background:
		TH1F *h_empty;
		double m_emptyRatio = 1.0, m_emptyRatioErr = 0.05;
		
		// Fit Functions:
		TF1 *f_fit, *f_empty;
		
		double m_nEmptyEtaPar = 0.0, m_nEmptyOmegaPar = 0.0;
		
		int m_nParameters = 0, m_nEmptyParameters = 0;
		
		Option_t *fitOption = "R0Q";
		
		vector<pair<double,double>> excludeRegions;
		vector<double> m_muFDC = {0.61, 0.55, 0.45};
		
		// Lineshapes:
		
		TH1F *h_etaLineshape;
		TF1  *f_etaLineshape;
		
		TH1F *h_etaPionLineshape;
		TF1  *f_etaPionLineshape;
		
		TH1F *h_omegaLineshape;
		TF1  *f_omegaLineshape;
		
		TH1F *h_fdcOmegaLineshape;
		TF1  *f_fdcOmegaLineshape;
		
		TF1  *f_chebyshev;
		
		void ZeroSignalPars(TF1*);
		void ZeroBkgdPars(TF1*);
		
		int InitializeFitFunction(TString funcName="");
		int InitializeFitFunction(TF1**, TString funcName="");
		int InitializeFitParameters();
		
		void   GuessOmegaParameters();
		void     FixOmegaParameters();
		void ReleaseOmegaParameters();
		
		void    GuessBkgdParameters();
		void      FixBkgdParameters();
		void  ReleaseBkgdParameters();
		
		void    GuessEtapParameters();
		/*
		void      FixEtapParameters();
		void  ReleaseEtapParameters();
		*/
		
		void     GuessEtaParameters();
		/*
		void       FixEtaParameters();
		void   ReleaseEtaParameters();
		*/
		
		void   GuessEmptyParameters();
		void     FixEmptyParameters();
		
		int InitializeEmptyFitFunction(TF1**, TString funcName="");
		int InitializeEmptyFitParameters();
		
		void   GuessEmptyBkgdParameters();
		void     FixEmptyBkgdParameters();
		void ReleaseEmptyBkgdParameters();
		
		void    GuessEmptyEtaParameters();
		void      FixEmptyEtaParameters();
		
		void  GuessEmptyOmegaParameters();
		void    FixEmptyOmegaParameters();
		
		void    GuessEmptyFDCParameters();
		
};

#endif
