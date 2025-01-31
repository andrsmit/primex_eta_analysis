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
			: h_data(nullptr), f_fit(nullptr), h_etaLineshape(nullptr), h_omegaLineshape(nullptr), 
				f_omegaLineshape(nullptr), f_chebyshev(nullptr) 
		{
			f_chebyshev = new TF1("chebyshev", "cheb5", 0.0, 1.0);
		}
		
		double binSize = 0.001;
		
		// Fitting options:
		
		int fitOption_signal = 1;
		int fitOption_bkgd   = 3;
		int fitOption_poly   = 1;
		int fitOption_omega  = 2;
		int fitOption_etap   = 1;
		
		double minFitRange=0.30, maxFitRange=1.10;
		double minMggCut = 0.50;
		double maxMggCut = 0.60;
		
		vector<pair<double,double>> excludeRegions;
		
		double mggFitFunction(double *x, double *par);
		
		void SetData(TH1F* h1) { h_data = h1; }
		void SetEtaLineshape(TH1F *h1) { 
			h_etaLineshape = h1;
			h_etaLineshape->Scale(1.0/h_etaLineshape->Integral());
		}
		void SetOmegaLineshape(TH1F *h1, int drawOption=0) { 
			h_omegaLineshape = h1;
			h_omegaLineshape->Scale(1.0/h_omegaLineshape->Integral());
			FitOmegaLineshape(drawOption);
		}
		void FitOmegaLineshape(int drawOption=0);
		void FitData();
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
			f1->SetLineColor(kGreen+2);
			f1->SetLineStyle(4);
			f1->SetNpx(1000);
			return f1;
		}
		
		//TH1F* GetPull();
		void FillPull(TH1F*);
		
		void GetYield(double&,double&,int useSignalPars=0);
		
		double IntegrateGaussian(double A, double mu, double sigma, double x1, double x2);
	
	private:
		// Data histogram to be fit:
		TH1F *h_data;
		
		// Function:
		TF1 *f_fit;
		
		Option_t *fitOption = "R0Q";
		
		// Lineshapes:
		
		TH1F *h_etaLineshape;
		TH1F *h_omegaLineshape;
		TF1  *f_omegaLineshape;
		TF1  *f_chebyshev;
		
		void ZeroSignalPars(TF1*);
		void ZeroBkgdPars(TF1*);
		
		int  InitializeFitFunction(TString funcName="");
		void InitializeFitFunction(TF1**, TString funcName="");
		void InitializeFitParameters();
		
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
		
};

#endif
