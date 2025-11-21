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

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "HFitInterface.h"
#include "Fit/Chi2FCN.h"
#include "Fit/ParameterSettings.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TVirtualFitter.h"
#include "TROOT.h"

#include "EtaAnalyzer.h"

class MggFitter {
	public:
		MggFitter();
		
		int phase = 3;
		double angle   = 0.0,   angleWidth = 0.0;
		double binSize = 0.001, emptyBinSize = 0.001;
		
		double empty_flux_ratio = 1.0;
		
		double incFraction_theory = 1.0;
		
		// Fitting options:
		
		int fitOption_signal =  2;
		int fitOption_bkgd   =  2;
		int fitOption_poly   =  2;
		int fitOption_omega  =  1;
		int fitOption_rho    =  0;
		int fitOption_etap   =  0;
		
		int useRawMass = 0;
		int vetoOption = 6;
		
		// Fitting options for empty target background:
		
		int emptyFitOption_eta   = 2;
		int emptyFitOption_omega = 1;
		int emptyFitOption_fdc   = 1;
		int emptyFitOption_bkgd  = 3;
		int emptyFitOption_poly  = 4;
		
		double lineshapeOffset = 0.0;
		
		double      minFitRange=0.300,      maxFitRange=0.950;
		double minEmptyFitRange=0.300, maxEmptyFitRange=0.950;
		
		double minMggCut = 0.5;
		double maxMggCut = 0.6;
		
		int combinedFit = 1;
		
		// histograms to store the fitted values of expected accidental counts:
		TH1F *h_acc_fit_full, *h_acc_fit_empty;
		
		// Set full and empty target histograms:
		
		void SetData(TH1F* h1, TH1F *h2)
		{
			h_full[0]           = h1;
			h_full[1]           = h2;
			h_acc_fit_full      = (TH1F*)h2->Clone("acc_full");
		}
		void SetEmpty(TH1F *h1, TH1F *h2, double ratio, double ratioErr)
		{ 
			h_empty[0]          = h1;
			h_empty[1]          = h2;
			m_emptyFluxRatio    = ratio;
			m_emptyFluxRatioErr = ratioErr;
			h_acc_fit_empty     = (TH1F*)h2->Clone("acc_empty");
		}
		void SetEmptyWide(TH1F *h1, TH1F *h2, double ratio)
		{
			h_emptyWide[0]      = h1;
			h_emptyWide[1]      = h2;
			m_emptyWideRatio    = ratio;
		}
		
		// Set/Fit lineshape of eta signal:
		void SetCohLineshape(TH1F *h1, int drawOption=0) { 
			h_cohLineshape = h1;
			h_cohLineshape->Scale(1.0/h_cohLineshape->Integral());
			FitCohLineshape(drawOption);
			return;
		}
		void FitCohLineshape(int drawOption=0);
		
		void SetQFLineshape(TH1F *h1, int drawOption=0) { 
			h_qfLineshape = h1;
			h_qfLineshape->Scale(1.0/h_qfLineshape->Integral());
			FitQFLineshape(drawOption);
			return;
		}
		void FitQFLineshape(int drawOption=0);
		
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
			h_omegaLineshape->Scale(1.0/h_omegaLineshape->Integral());
			if(fitOption_omega>2) FitOmegaLineshape(drawOption);
			return;
		}
		void FitOmegaLineshape(int drawOption=0);
		
		// Set lineshape of rho background:
		void SetRhoLineshape(TH1F *h1, int drawOption=0) { 
			h_rhoLineshape = h1;
			h_rhoLineshape->Scale(1.0/h_rhoLineshape->Integral());
			return;
		}
		
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
		void FitEmptyWide();
		
		TF1* GetFitFunction() { 
			f_full->SetLineColor(kGreen);
			f_full->SetNpx(1000);
			return f_full;
		}
		void GetFitFunction_Signal(TF1** f1, TString fname="signalFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroHadronicBkgdPars(locf1);
			ZeroOmegaPars(locf1);
			ZeroRhoPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			
			locf1->SetLineColor(kCyan+2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_Empty(TF1** f1, TString fname="signalFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1);
			ZeroOmegaPars(locf1);
			ZeroRhoPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			
			locf1->SetLineColor(kRed);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_Bkgd(TF1** f1, TString fname="bkgdFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1);
			ZeroOmegaPars(locf1);
			ZeroRhoPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			
			locf1->SetLineColor(kMagenta);
			locf1->SetLineStyle(4);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_Omega(TF1** f1, TString fname="omegaFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1);
			ZeroRhoPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			ZeroRhoPars(locf1);
			
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			
			locf1->SetLineColor(kGreen);
			locf1->SetLineStyle(7);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_Rho(TF1** f1, TString fname="rhoFit") {
			if(fitOption_rho==0) return;
			
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1);
			ZeroOmegaPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			
			locf1->SetLineColor(kGreen+2);
			locf1->SetLineStyle(2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_HadronicBkgd(TF1** f1, TString fname="hadronicBkgdFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1, 1);
			ZeroOmegaPars(locf1);
			ZeroRhoPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			locf1->SetParameter("A_{#eta#pi}",0.0);
			
			locf1->SetLineColor(kBlue+2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetFitFunction_EtaPion(TF1** f1, TString fname="etaPionFit") {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1, 1);
			ZeroOmegaPars(locf1);
			ZeroRhoPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			ZeroAccPars(locf1);
			//locf1->SetParameter("#alpha_{acc,switch}",0.0);
			locf1->SetParameter("A_{#eta#pi#pi}",0.0);
			
			locf1->SetLineColor(kCyan);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetEmptyFitFunction(TF1** f1, TString fname="emptyFit", double binWidth=0.0, int removeAcc=1) {
			TF1 *locf1;
			InitializeEmptyFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_empty->GetParameters());
			
			locf1->SetLineColor(kRed);
			locf1->SetLineStyle(2);
			locf1->SetLineWidth(2);
			locf1->SetNpx(1000);
			
			if(binWidth==0.0) locf1->SetParameter(locf1->GetNpar()-1, binSize);
			else              locf1->SetParameter(locf1->GetNpar()-1, binWidth);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetEmptyBGFitFunction(TF1** f1, TString fname="emptyFit", double binWidth=0.0, int removeAcc=1) {
			TF1 *locf1;
			InitializeEmptyFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_empty->GetParameters());
			
			locf1->SetParameter("A_{empty}",0.0);
			ZeroEmptyFDCPars(locf1);
			
			locf1->SetLineColor(kMagenta);
			locf1->SetLineStyle(4);
			locf1->SetLineWidth(2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetAccFunction(TF1** f1, TString fname="accFit", double binWidth=0.0) {
			TF1 *locf1;
			InitializeFitFunction(&locf1,"locf1");
			locf1->SetParameters(f_full->GetParameters());
			
			ZeroSignalPars(locf1);
			ZeroOmegaPars(locf1);
			ZeroBkgdPars(locf1);
			ZeroEtaPrimePars(locf1);
			ZeroEmptyPars(locf1);
			
			locf1->SetLineColor(kBlack);
			locf1->SetLineStyle(4);
			locf1->SetLineWidth(2);
			locf1->SetNpx(1000);
			
			*f1 = (TF1*)locf1->Clone(fname.Data());
			delete locf1;
			return;
		}
		void GetAccFitRatio(TH1F** h1, int isEmpty=0, TString hname="accFitRatio") {
			switch(isEmpty) {
				case 0:
				{
					*h1 = (TH1F*)h_acc_fit_full->Clone(hname.Data());
					for(int ibin=1; ibin<(*h1)->GetXaxis()->GetNbins(); ibin++) {
						(*h1)->SetBinContent(ibin, 1.0);
						if(h_full[1]->GetBinContent(ibin)>0.0) {
							(*h1)->SetBinContent(ibin, h_acc_fit_full->GetBinContent(ibin)/h_full[1]->GetBinContent(ibin));
						}
					}
					(*h1)->GetXaxis()->SetRangeUser(minFitRange, maxFitRange);
					(*h1)->SetLineColor(kBlue);
					break;
				}
				case 1:
				{
					*h1 = (TH1F*)h_acc_fit_empty->Clone(hname.Data());
					for(int ibin=1; ibin<(*h1)->GetXaxis()->GetNbins(); ibin++) {
						(*h1)->SetBinContent(ibin, 1.0);
						if(h_empty[1]->GetBinContent(ibin)>0.0) {
							(*h1)->SetBinContent(ibin, h_acc_fit_empty->GetBinContent(ibin)/h_empty[1]->GetBinContent(ibin));
						}
					}
					(*h1)->GetXaxis()->SetRangeUser(minFitRange, maxFitRange);
					(*h1)->SetLineColor(kRed);
					break;
				}
				default:
				{
					h1 = nullptr;
					break;
				}
			}
			return;
		}
		
		void FillPull(TH1F*);
		void FillEmptyPull(TH1F*);
		
		void GetYield(double&, double&, int useSignalPars=0, int subtractEtaPion=0, int verbose=0);
		void GetEmptyYield(double&, double&);
		void GetEtaPionYield(double&, double&);
		void GetHadronicBkgdYield(double&, double&);
		void GetOmegaYield(double&, double&);
		void GetRhoYield(double&, double&);
		void GetBkgdYield(double&, double&);
		
		void GetOmegaFitPars(double&, double&, double&, double&, double&, double&, double&, double&);
		
		void GetEmptyEtaFraction(double&, double&);
		void GetHadronicBkgdFraction(double&, double&);
		void GetEtaPionBkgdFraction(double&, double&);
		
		void GetLineshapeShift(double&, double&);
		void GetQFFraction(double&, double&);
		
		double IntegrateGaussian(double mu, double sigma, double x1, double x2);
		
		void DumpFitParameters();
		void DumpEmptyFitParameters();
		
		// Fit Functions:
		double      MggFitFunction(double *x, double *par);
		double EmptyMggFitFunction(double *x, double *par);
		double        EtaLineshape(double *x, double *par);
		
		//-----------------------------------------------------------------//
		// Models for each fit component (defined in Models.cc):
		
		double model_eta(double, double*);
		double model_omega(double, double*);
		double model_em(double, double*);
		double model_etap(double, double*);
		double model_beamline(double, double*);
		
		struct CombinedNLL;
	
	private:
		
		// Full and empty target histograms to be fit (0: prompt, 1: out-of-time):
		TH1F *h_full[2];
		TH1F *h_empty[2];
		TH1F *h_emptyWide[2];
		
		// Fit Functions:
		TF1 *f_full, *f_empty, *f_emptyWide;
		
		ROOT::Fit::FitResult result;
		
		// Scaling Factors:
		
		double m_emptyFluxRatio = 1.0, m_emptyFluxRatioErr = 0.05;
		double m_emptyWideRatio = 1.0;
		
		
		// Fit parameters and index arrays:
		vector<TString> m_parametersFull;
		vector<int> m_parIndexFull, m_parIndexEmpty;
		
		// should be called as soon as object is configured (fitting options set):
		void InitializeParameterArrays();
		
		int GetSignalParameters(vector<TString>&);
		int GetOmegaParameters(vector<TString>&);
		int GetEMParameters(vector<TString>&);
		int GetEtaPrimeParameters(vector<TString>&);
		int GetBeamlineParameters(vector<TString>&);
		
		// Predictions from BGGEN:
		double m_hadronicBkgdYieldBGGEN = 0.0, m_etaPionYieldBGGEN   = 0.0;
		double m_hadronicBkgdFrac       = 0.0, m_hadronicBkgdFracErr = 0.0;
		double m_etaPionBkgdFrac        = 0.0, m_etaPionBkgdFracErr  = 0.0;
		
		Option_t *fitOption = "R0Q";
		
		vector<pair<double,double>> excludeRegions;
		vector<double> m_muFDC       = {0.750, 0.619, 0.532, 0.454, 0.399, 0.339};
		vector<double> m_muFDC_omega = {0.613, 0.535, 0.462};
		vector<double> m_muFDC_eta   = {0.445, 0.406, 0.367};
		
		// MC Lineshapes:
		
		TH1F *h_cohLineshape, *h_qfLineshape;
		TF1  *f_cohLineshape, *f_qfLineshape;
		TF1  *f_etaLineshape;
		
		TH1F *h_hadronicBkgdLineshape;
		TF1  *f_hadronicBkgdLineshape;
		
		TH1F *h_etaPionLineshape;
		TF1  *f_etaPionLineshape;
		
		TH1F *h_omegaLineshape;
		TF1  *f_omegaLineshape;
		
		TH1F *h_rhoLineshape;
		TF1  *f_rhoLineshape;
		
		TF1  *f_chebyshev;
		
		void CheckBinSize(TH1F *h1, TString histTitle="Hist");
		
		// Functions to "turn off" each component of the fit function for purposes of
		// drawing, or estimating yields:
		
		void ZeroSignalPars(TF1*, int excludeHadBkgd=0);
		void ZeroHadronicBkgdPars(TF1*);
		void ZeroOmegaPars(TF1*);
		void ZeroRhoPars(TF1*);
		void ZeroBkgdPars(TF1*);
		void ZeroEtaPrimePars(TF1 *f1);
		void ZeroEmptyPars(TF1*);
		void ZeroAccPars(TF1*);
		void ZeroEmptyBkgdPars(TF1 *f1);
		void ZeroEmptyFDCPars(TF1 *f1);
		void ZeroEmptyEtaPars(TF1 *f1);
		void ZeroEmptyOmegaPars(TF1 *f1);
		
		// Routines to setup fit functions and initialize parameter settings:
		
		int InitializeFitFunction(TString funcName="");
		int InitializeFitFunction(TF1**, TString funcName="");
		
		void InitializeFitParameters(ROOT::Fit::Fitter&);
		void SetFitParameters(ROOT::Fit::Fitter&, ROOT::Fit::Fitter&);
		
		void GuessEtaParameters(vector<double>&);
		void GuessOmegaParameters(vector<double>&);
		void GuessEMParameters(vector<double>&);
		void GuessEtaPrimeParameters(vector<double>&);
		void GuessBeamlineParameters(vector<double>&);
		
		//-----------------------------------------------------------------//
		// Functions for modifying which parameters can float in combined fit:
		
		void FixOmegaParameters(ROOT::Fit::Fitter&);
		void ReleaseOmegaParameters(ROOT::Fit::Fitter&);
		
		void FixEMParameters(ROOT::Fit::Fitter&);
		void ReleaseEMParameters(ROOT::Fit::Fitter&);
		
		void FixBeamlineParameters(ROOT::Fit::Fitter&);
		void ReleaseBeamlineParameters(ROOT::Fit::Fitter&);
		void ReleaseFDCParameters(ROOT::Fit::Fitter&);
		
		void FixEtaParameters(ROOT::Fit::Fitter&);
		void ReleaseEtaParameters(ROOT::Fit::Fitter&);
		
		//-----------------------------------------------------------------//
		
		void UpdateFitFunctions(ROOT::Fit::FitResult);
		void UpdateFitFunctions(ROOT::Fit::Fitter&);
		
		int InitializeEmptyFitFunction(TF1**, TString funcName="");
		int InitializeEmptyWideFitParameters();
		
		void   GuessEmptyBkgdParameters();
		void     FixEmptyBkgdParameters();
		void ReleaseEmptyBkgdParameters();
		
		void    GuessEmptyFDCParameters();
		void      FixEmptyFDCParameters();
		void  ReleaseEmptyFDCParameters();
};

#endif
