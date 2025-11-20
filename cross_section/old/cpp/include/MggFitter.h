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
				h_hadronicBkgdLineshape(nullptr), f_hadronicBkgdLineshape(nullptr), 
				h_etaLineshape(nullptr),          f_etaLineshape(nullptr),
				h_omegaLineshape(nullptr),        f_omegaLineshape(nullptr),
				h_rhoLineshape(nullptr), 
				h_etaPionLineshape(nullptr),      f_etaPionLineshape(nullptr), 
				h_fdcOmegaLineshape(nullptr),     f_fdcOmegaLineshape(nullptr)
		{
			f_chebyshev = new TF1("chebyshev", "cheb5", 0.0, 1.0);
		}
		
		double angle   = 0.0,   angleWidth = 0.0;
		double binSize = 0.001, emptyBinSize = 0.001;
		
		// Fitting options:
		
		int fitOption_signal = 1;
		int fitOption_bkgd   = 3;
		int fitOption_poly   = 1;
		int fitOption_omega  = 2;
		int fitOption_etap   = 1;
		
		int useRawMass = 0;
		int vetoOption = 6;
		
		// Fitting options for empty target background:
		
		int fitOption_empty = 0;
		
		int emptyFitOption_eta   = 0;
		int emptyFitOption_omega = 0;
		int emptyFitOption_fdc   = 0;
		int emptyFitOption_bkgd  = 2;
		int emptyFitOption_poly  = 3; 
		
		double      minFitRange=0.300,      maxFitRange=1.100;
		double minEmptyFitRange=0.305, maxEmptyFitRange=0.950;
		
		double minMggCut = 0.5;
		double maxMggCut = 0.6;
		
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
		void SetHadronicBkgdLineshape(TH1F *h1, int drawOption=0) {
			m_hadronicBkgdYieldBGGEN = h1->Integral();
			h_hadronicBkgdLineshape = h1;
			if(m_hadronicBkgdYieldBGGEN>1.0) {
				h_hadronicBkgdLineshape->Scale(1.0/m_hadronicBkgdYieldBGGEN);
			}
			FitHadronicBkgdLineshape(drawOption);
			return;
		}
		void FitHadronicBkgdLineshape(int drawOption=0);
		
		// Set/fit lineshape of eta+pion background:
		void SetEtaPionLineshape(TH1F *h1, int drawOption=0) {
			m_etaPionYieldBGGEN = h1->Integral();
			h_etaPionLineshape = h1;
			if(m_etaPionYieldBGGEN>1.0) {
				h_etaPionLineshape->Scale(1.0/m_etaPionYieldBGGEN);
			}
			FitEtaPionLineshape(drawOption);
			return;
		}
		void FitEtaPionLineshape(int drawOption=0);
		
		// Set/fit lineshape of omega background:
		void SetOmegaLineshape(TH1F *h1, int drawOption=0) { 
			h_omegaLineshape = h1;
			if((fitOption_omega<3) || (fitOption_omega>4)) {
				h_omegaLineshape->Scale(1.0/h_omegaLineshape->Integral());
				FitOmegaLineshape(drawOption);
			}
			return;
		}
		void FitOmegaLineshape(int drawOption=0);
		
		// Set lineshape of rho background:
		void SetRhoLineshape(TH1F *h1, int drawOption=0) { 
			h_rhoLineshape = h1;
			h_rhoLineshape->Scale(1.0/h_rhoLineshape->Integral());
			return;
		}
		
		// Set/fit lineshape of omegas form 1st FDC package:
		void SetFDCOmegaLineshape(TH1F *h1, int drawOption=0) { 
			h_fdcOmegaLineshape = h1;
			h_fdcOmegaLineshape->Scale(1.0/h_fdcOmegaLineshape->Integral());
			FitFDCOmegaLineshape(drawOption);
			return;
		}
		void FitFDCOmegaLineshape(int drawOption=0);
		
		void SetHadronicBkgdFraction(double frac, double fracErr) { 
			m_hadronicBkgdFrac    = frac;
			m_hadronicBkgdFracErr = fracErr;
			return;
		}
		void SetEtaPionBkgdFraction(double frac, double fracErr) { 
			m_etaPionBkgdFrac    = frac;
			m_etaPionBkgdFracErr = fracErr;
			return;
		}
		
		void FitData();
		void FitEmpty();
		void FitDataWithEmpty();
		void DrawFitResults(TCanvas&);
		
		TF1* GetFitFunction() { 
			f_fit->SetLineColor(kGreen);
			f_fit->SetNpx(1000);
			return f_fit;
		}
		void GetSignalFunction(TF1** f1, TString fname="signalFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroBkgdPars(locf1);
			locf1->SetLineColor(kCyan+2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetBkgdFunction(TF1** f1, TString fname="bkgdFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroSignalPars(locf1);
			locf1->SetParameter("N_{#omega}", 0.0);
			locf1->SetParameter("N_{empty}", 0.0);
			locf1->SetLineColor(kMagenta);
			locf1->SetLineStyle(4);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetOmegaFunction(TF1** f1, TString fname="omegaFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroSignalPars(locf1);
			switch(fitOption_bkgd) {
				case 1:
					for(int ipar=0; ipar<=fitOption_poly; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
				case 2:
					for(int ipar=0; ipar<5; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
				case 3:
					for(int ipar=0; ipar<=fitOption_poly; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
			}
			locf1->SetParameter("N_{empty}", 0.0);
			locf1->SetLineColor(kGreen);
			locf1->SetLineStyle(4);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetHadronicBkgdFunction(TF1** f1, TString fname="hadronicBkgdFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroBkgdPars(locf1);
			switch(fitOption_signal) {
				case 6:
					locf1->SetParameter("N_{#eta}",0.0);
					locf1->SetParameter("N_{#eta,bkgd}", f_fit->GetParameter("N_{#eta,bkgd}"));
					break;
				case 7:
					locf1->SetParameter("frac_{#eta}", 0.0);
					locf1->SetParameter("frac_{bkgd}", f_fit->GetParameter("frac_{bkgd}"));
					break;
				case 8:
					locf1->SetParameter("frac_{#eta}", 0.0);
					locf1->SetParameter("frac_{bkgd}", f_fit->GetParameter("frac_{bkgd}"));
					break;
				case 9:
					locf1->SetParameter("frac_{#eta}",    0.0);
					locf1->SetParameter("frac_{#eta#pi}", 0.0);
					locf1->SetParameter("frac_{bkgd}",    f_fit->GetParameter("frac_{bkgd}"));
					break;
				case 10:
					locf1->SetParameter("frac_{#eta}",    0.0);
					locf1->SetParameter("frac_{#eta#pi}", 0.0);
					locf1->SetParameter("frac_{bkgd}",    f_fit->GetParameter("frac_{bkgd}"));
					break;
				case 11:
					locf1->SetParameter("N_{#eta}",    0.0);
					locf1->SetParameter("A_{#eta#pi}", 0.0);
					locf1->SetParameter("A_{#eta#pi#pi}",    f_fit->GetParameter("A_{#eta#pi#pi}"));
					break;
				case 12:
					locf1->SetParameter("N_{#eta}",    0.0);
					locf1->SetParameter("A_{#eta#pi}", 0.0);
					locf1->SetParameter("A_{#eta#pi#pi}",    f_fit->GetParameter("A_{#eta#pi#pi}"));
					break;
			}
			locf1->SetLineColor(kBlue+2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetEtaPionFunction(TF1** f1, TString fname="etaPionFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroBkgdPars(locf1);
			switch(fitOption_signal) {
				case 9:
					locf1->SetParameter("frac_{#eta}",    0.0);
					locf1->SetParameter("frac_{#eta#pi}", f_fit->GetParameter("frac_{#eta#pi}"));
					locf1->SetParameter("frac_{bkgd}",    0.0);
					break;
				case 10:
					locf1->SetParameter("frac_{#eta}",    0.0);
					locf1->SetParameter("frac_{#eta#pi}", f_fit->GetParameter("frac_{#eta#pi}"));
					locf1->SetParameter("frac_{bkgd}",    0.0);
					break;
				case 11:
					locf1->SetParameter("N_{#eta}",    0.0);
					locf1->SetParameter("A_{#eta#pi}", f_fit->GetParameter("A_{#eta#pi}"));
					locf1->SetParameter("A_{#eta#pi#pi}",    0.0);
					break;
				case 12:
					locf1->SetParameter("N_{#eta}",    0.0);
					locf1->SetParameter("A_{#eta#pi}", f_fit->GetParameter("A_{#eta#pi}"));
					locf1->SetParameter("A_{#eta#pi#pi}",    0.0);
					break;
			}
			locf1->SetLineColor(kCyan);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetEmptyFitFunction(TF1** f1, TString fname="emptyFit", double binWidth=0.0) {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_fit->GetParameters());
			ZeroSignalPars(locf1);
			switch(fitOption_bkgd) {
				case 1:
					for(int ipar=0; ipar<=fitOption_poly; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
				case 2:
					for(int ipar=0; ipar<5; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
				case 3:
					for(int ipar=0; ipar<=fitOption_poly; ipar++) {
						locf1->SetParameter(Form("p%d",ipar), 0.0);
					}
					break;
			}
			locf1->SetParameter("N_{#omega}", 0.0);
			locf1->SetLineColor(kRed);
			locf1->SetLineStyle(4);
			locf1->SetLineWidth(2);
			locf1->SetNpx(1000);
			
			if(binWidth==0.0) locf1->SetParameter(locf1->GetNpar()-1, binSize);
			else              locf1->SetParameter(locf1->GetNpar()-1, binWidth);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetEmptyBGFitFunction(TF1** f1, TString fname="emptyFit", double binWidth=0.0) {
			TF1 *locf1;
			InitializeEmptyFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_empty->GetParameters());
			locf1->SetLineColor(kMagenta);
			locf1->SetLineStyle(4);
			locf1->SetLineWidth(2);
			locf1->SetNpx(1000);
			
			if(binWidth==0.0) locf1->SetParameter(locf1->GetNpar()-1, binSize);
			else              locf1->SetParameter(locf1->GetNpar()-1, binWidth);
			
			locf1->SetParameter("N_{#eta}", 0.0);
			locf1->SetParameter("N_{#omega}", 0.0);
			locf1->SetParameter("N_{fdc}", 0.0);
			for(int ifdc=0; ifdc<m_muFDC.size();       ifdc++) locf1->SetParameter(Form("N_{fdc,%d}",       ifdc+1), 0.0);
			for(int ifdc=0; ifdc<m_muFDC_omega.size(); ifdc++) locf1->SetParameter(Form("N_{fdc,#omega,%d}",ifdc+1), 0.0);
			for(int ifdc=0; ifdc<m_muFDC_eta.size();   ifdc++) locf1->SetParameter(Form("N_{fdc,#eta,%d}",  ifdc+1), 0.0);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		
		double GetEmptyEtaFitPar() { return m_nEmptyEtaPar; }
		double GetEmptyOmegaFitPar() { return m_nEmptyOmegaPar; }
		
		void FillPull(TH1F*);
		void FillEmptyPull(TH1F*);
		
		void GetYield(double&, double&, int useSignalPars=0, int subtractEtaPion=0);
		void GetEmptyYield(double&, double&, int excludeNonPeaking=0);
		void GetOmegaYield(double&, double&);
		void GetBkgdYield(double&, double&);
		void GetHadronicBkgdYield(double&, double&);
		void GetEtaPionYield(double&, double&);
		
		void GetOmegaFitPars(double&, double&, double&, double&, double&, double&, double&, double&);
		
		void GetEmptyEtaFraction(double&, double&);
		void GetHadronicBkgdFraction(double&, double&);
		void GetEtaPionBkgdFraction(double&, double&);
		
		void GetLineshapeShift(double&, double&);
		
		double IntegrateGaussian(double mu, double sigma, double x1, double x2);
		
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
		
		double m_hadronicBkgdFrac = 0.0, m_hadronicBkgdFracErr = 0.0;
		double m_etaPionBkgdFrac  = 0.0, m_etaPionBkgdFracErr  = 0.0;
		
		double m_hadronicBkgdYieldBGGEN = 0.0, m_etaPionYieldBGGEN = 0.0;
		
		int m_nParameters = 0, m_nEmptyParameters = 0;
		
		Option_t *fitOption = "R0Q";
		
		vector<pair<double,double>> excludeRegions;
		vector<double> m_muFDC       = {0.619, 0.532, 0.454, 0.399, 0.339};
		vector<double> m_muFDC_omega = {0.613, 0.535, 0.462};
		vector<double> m_muFDC_eta   = {0.445, 0.406, 0.367};
		// Lineshapes:
		
		TH1F *h_etaLineshape;
		TF1  *f_etaLineshape;
		
		TH1F *h_hadronicBkgdLineshape;
		TF1  *f_hadronicBkgdLineshape;
		
		TH1F *h_etaPionLineshape;
		TF1  *f_etaPionLineshape;
		
		TH1F *h_omegaLineshape;
		TF1  *f_omegaLineshape;
		
		TH1F *h_rhoLineshape;
		
		TH1F *h_fdcOmegaLineshape;
		TF1  *f_fdcOmegaLineshape;
		
		TF1  *f_chebyshev;
		
		void CheckBinSize(TH1F *h1, TString histTitle="Hist");
		
		void ZeroSignalPars(TF1*, int subtractEtaPion=0);
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
		void  ReleaseEmptyEtaParameters();
		
		void  GuessEmptyOmegaParameters();
		void    FixEmptyOmegaParameters();
		
		void    GuessEmptyFDCParameters();
		void      FixEmptyFDCParameters();
		void  ReleaseEmptyFDCParameters();
};

#endif
