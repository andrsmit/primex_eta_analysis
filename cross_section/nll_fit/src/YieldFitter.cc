#include "CrossSection.h"
#include "YieldFitter.h"

int InitializeYieldFitter(YieldFitter &fitter, EtaAnalyzer &anaObj)
{
	double locBinSize, locMin, locMax;
	
	anaObj.GetBeamEnergyBinning(locBinSize, locMin, locMax);
	fitter.SetBeamEnergyBinning(locBinSize, locMin, locMax);
	
	anaObj.GetReconAngleBinning(locBinSize, locMin, locMax);
	fitter.SetReconAngleBinning(locBinSize, locMin, locMax);
	
	anaObj.GetThrownAngleBinning(locBinSize, locMin, locMax);
	fitter.SetThrownAngleBinning(locBinSize, locMin, locMax);
	
	fitter.SetLuminosity(anaObj.GetLuminosity());
	
	// check that the angular matrices have been loaded. If not, load them:
	if(!anaObj.IsMatrixLoaded()) {
		if(anaObj.LoadAngularMatrix()) {
			printf("\n\nProblem loading angular matrices.\n\n");
			return 1;
		}
	}
	
	fitter.SetAngularMatrixFull((TH3F*)anaObj.GetAngularMatrixFine());
	fitter.SetThrown((TH2F*)anaObj.GetThrown());
	fitter.SetFluxWeights((TH1F*)anaObj.GetFluxWeights());
	
	fitter.SetYield((TH1F*)anaObj.GetAngularYield(1));
	
	return 0;
}

struct YieldFitter::CombinedChi2 {
	
	YieldFitter *fitter;
	
	double x1, x2;
	int nHists;
	vector<vector<int>> parIndices;
	
	CombinedChi2(YieldFitter *f, const vector<vector<int>> &indices, double a, double b) 
		: fitter(f), parIndices(indices), x1(a), x2(b) 
	{
		nHists = (int)indices.size();
		//printf("%d HISTS INCLUDED IN FIT\n", nHists);
	}
	CombinedChi2() {};
	
	double operator()(const double* p) const
	{
		vector<vector<double>> pars(nHists);
		for(int ih=0; ih<nHists; ih++) {
			pars[ih].clear();
			for(int ipar=0; ipar<parIndices[ih].size(); ipar++) {
				pars[ih].push_back(p[parIndices[ih][ipar]]);
			}
		}
		
		double chi2 = 0.0;
		
		// loop over all histograms and compute chi2:
		
		for(int ih=0; ih<nHists; ih++) {
			
			// loop over all bins within fit range:
			
			TAxis *locAxis = fitter->h_dNdTheta[ih]->GetXaxis();
			
			int locNbins = locAxis->GetNbins(); 
			double locBinSize = locAxis->GetBinWidth(1);
			
			double minEnergy = fitter->m_energyBins[ih].first;
			double maxEnergy = fitter->m_energyBins[ih].second;
			
			int minBin = locAxis->FindBin(x1);
			int maxBin = locAxis->FindBin(x2);
			
			for(int ibin=minBin; ibin<=maxBin; ibin++) {
				double minAngle = locAxis->GetBinCenter(ibin) - 0.5*locAxis->GetBinWidth(ibin);
				double maxAngle = locAxis->GetBinCenter(ibin) + 0.5*locAxis->GetBinWidth(ibin);
				
				double measYield    = fitter->h_dNdTheta[ih]->GetBinContent(ibin);
				double measYieldErr = fitter->h_dNdTheta[ih]->GetBinError(ibin);
				if(measYieldErr < sqrt(measYield)) measYieldErr = sqrt(measYield);
				
				double expYield = fitter->GetExpectedYield(minAngle, maxAngle, ih, pars[ih].data());
				/*
				printf("%.2f-%.2f  (%.1f GeV - %.1f GeV):\n", minAngle, maxAngle, minEnergy, maxEnergy);
				printf("  Measured Yield: %f +/- %f\n", measYield, measYieldErr);
				printf("  Expected Yield: %f\n", expYield);
				*/
				chi2 += pow((measYield - expYield)/measYieldErr, 2.0);
			}
		}
		
		return chi2;
	}
};


void YieldFitter::InitializeParameterArrays()
{
	/*------------------------------------------------------------------------------------------------------//
	
	This function initializes the private member vectors:
		- 'm_parameterList': std::vector storing the full list of parameter names
		- 'm_parIndices'   : std::vector<vector<int>> - Each element is a vector asociated with a specific yield distribution.
								That vector then lists the index of the relevant parameters for that particular histogram
								inside of the full list (m_parameterList).
	
	The full list of parameters includes:
		
		- Gamma (the same across all distributions)
		
		- A_coh (unique values for each histogram)
		- Phi (unique values for each histogram)
		- A_inc (unique values for each histogram)
		
		- Depending on the model used, there could be an additional paramter (associated with quasifree-neutron)
			(unique for all histograms)
	
	------------------------------------------------------------------------------------------------------*/
	
	int nParameters = 0;
	
	m_parameterList.push_back("#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	nParameters++;
	
	m_parameterList.push_back("#phi[#circ]");
	nParameters++;
	
	int nHists = (int)h_dNdTheta.size();
	for(int ih=0; ih<nHists; ih++)
	{
		vector<int> locParIndices;
		locParIndices.clear();
		
		locParIndices.push_back(0);
		locParIndices.push_back(1);
		
		locParIndices.push_back(nParameters+0);
		m_parameterList.push_back(Form("A_{Coh,%d}",ih));
		
		if(m_components.size()==4)
		{
			locParIndices.push_back(nParameters+1);
			locParIndices.push_back(nParameters+2);
			
			m_parameterList.push_back(Form("A_{QFP,%d}",ih));
			m_parameterList.push_back(Form("A_{QFN,%d}",ih));
			nParameters += 3;
		}
		else
		{
			locParIndices.push_back(nParameters+1);
			
			m_parameterList.push_back(Form("A_{Inc,%d}",ih));
			nParameters += 2;
		}
		
		m_parIndices.push_back(locParIndices);
	}
}

void YieldFitter::InitializeFitParameters(ROOT::Fit::Fitter &fitter)
{
	int nParameters      = 0;
	int nTotalParameters = (int)m_parameterList.size();
	
	// Gamma:
	fitter.Config().ParSettings(nParameters).Set(m_parameterList[nParameters].Data(), f_dNdTheta[0]->GetParameter(1));
	fitter.Config().ParSettings(nParameters).Release();
	fitter.Config().ParSettings(nParameters).SetLimits(0.1, 1.0);
	nParameters++;
	
	// Phi:
	fitter.Config().ParSettings(nParameters).Set(m_parameterList[nParameters].Data(), f_dNdTheta[0]->GetParameter(2));
	fitter.Config().ParSettings(nParameters).Release();
	fitter.Config().ParSettings(nParameters).SetLimits(-90.0, 360.0);
	nParameters++;
	
	for(int ih=0; ih<(int)f_dNdTheta.size(); ih++) {
		for(int ipar=3; ipar<f_dNdTheta[ih]->GetNpar(); ipar++) {
			fitter.Config().ParSettings(nParameters).Set(m_parameterList[nParameters].Data(), f_dNdTheta[ih]->GetParameter(ipar));
			
			int minLimit = 0.0, maxLimit = 2.5;
			fitter.Config().ParSettings(nParameters).Release();
			fitter.Config().ParSettings(nParameters).SetLimits(minLimit, maxLimit);
			nParameters++;
		}
	}
}

void YieldFitter::DumpFitParameters(ROOT::Fit::Fitter &fitter)
{
	for(int ipar=0; ipar<(int)m_parameterList.size(); ipar++) {
		printf("  p%d: %f (%s)\n", ipar, fitter.Config().ParSettings(ipar).Value(), m_parameterList[ipar].Data());
	}
}

void YieldFitter::UpdateFitFunctions(ROOT::Fit::FitResult result)
{
	//=======================================================================================================//
	// Copy the parameters & associated errors from combined fit over to the 'f_dNdTheta' objects:
	
	int nParameters = 0;
	
	int nHists = (int)h_dNdTheta.size();
	for(int ih=0; ih<nHists; ih++)
	{
		f_dNdTheta[ih]->SetParameter(1, result.Parameter(0));
		f_dNdTheta[ih]->SetParError(1, result.ParError(0));
		
		f_dNdTheta[ih]->SetParameter(2, result.Parameter(1));
		f_dNdTheta[ih]->SetParError(2, result.ParError(1));
		
		if(ih==0) nParameters += 2;
		
		for(int ipar=3; ipar<(int)f_dNdTheta[ih]->GetNpar(); ipar++) {
			//printf("f_dNdTheta[%d]->SetParameter(%d, result.GetParameter(%d))  (%s)\n", ih, ipar, nParameters, m_parameterList[nParameters].Data());
			f_dNdTheta[ih]->SetParameter(ipar, result.Parameter(nParameters));
			f_dNdTheta[ih]->SetParError(ipar, result.ParError(nParameters));
			nParameters++;
		}
	}
}

void YieldFitter::InitializeMatrices()
{
	h_matrices.clear();
	
	int nHists = (int)h_dNdTheta.size();
	
	for(int ih=0; ih<nHists; ih++) {
		
		TH3F *loch3 = (TH3F*)h_matrixFull->Clone(Form("h_matrix_%d\n",ih));
		
		int minEBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[ih].first  + (1.e-6));
		int maxEBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[ih].second - (1.e-6));
		
		for(int iEnergyBin=1; iEnergyBin<=h_matrixFull->GetZaxis()->GetNbins(); iEnergyBin++) {
			
			if((iEnergyBin<minEBin) || (iEnergyBin>maxEBin)) {
				for(int iThetaBin=1; iThetaBin<=h_matrixFull->GetXaxis()->GetNbins(); iThetaBin++) {
					for(int iRecBin=1; iRecBin<=h_matrixFull->GetYaxis()->GetNbins(); iRecBin++) {
						loch3->SetBinContent(iThetaBin, iRecBin, iEnergyBin, 0.0);
						loch3->SetBinError(iThetaBin, iRecBin, iEnergyBin, 0.0);
					}
				}
				continue;
			}
			
			double locEnergy = h_matrixFull->GetZaxis()->GetBinCenter(iEnergyBin);
			
			for(int iThetaBin=1; iThetaBin<=h_matrixFull->GetXaxis()->GetNbins(); iThetaBin++) {
				double locTheta = h_matrixFull->GetXaxis()->GetBinCenter(iThetaBin);
				
				double locThrown = h_thrown->GetBinContent(h_thrown->GetXaxis()->FindBin(locEnergy), 
					h_thrown->GetYaxis()->FindBin(locTheta));
				
				for(int iRecBin=1; iRecBin<=h_matrixFull->GetYaxis()->GetNbins(); iRecBin++) {
					double locAcc = 0.0, locAccErr = 0.0;
					if(locThrown > 0.0) {
						locAcc    = h_matrixFull->GetBinContent(iThetaBin, iRecBin, iEnergyBin) / locThrown;
						locAccErr = sqrt(locThrown * locAcc * (1.0-locAcc)) / locThrown;
					}
					loch3->SetBinContent(iThetaBin, iRecBin, iEnergyBin, locAcc);
					loch3->SetBinError(iThetaBin, iRecBin, iEnergyBin, locAccErr);
				}
			}
		}
		
		h_matrices.push_back(loch3);
	}
	
	int minEBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[0].first  + (1.e-6));
	int maxEBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[nHists-1].second - (1.e-6));
	
	for(int iEnergyBin=1; iEnergyBin<=h_matrixFull->GetZaxis()->GetNbins(); iEnergyBin++) {
		
		if((iEnergyBin<minEBin) || (iEnergyBin>maxEBin)) {
			for(int iThetaBin=1; iThetaBin<=h_matrixFull->GetXaxis()->GetNbins(); iThetaBin++) {
				for(int iRecBin=1; iRecBin<=h_matrixFull->GetYaxis()->GetNbins(); iRecBin++) {
					h_matrixFull->SetBinContent(iThetaBin, iRecBin, iEnergyBin, 0.0);
					h_matrixFull->SetBinError(iThetaBin, iRecBin, iEnergyBin, 0.0);
				}
			}
			continue;
		}
		
		double locEnergy = h_matrixFull->GetZaxis()->GetBinCenter(iEnergyBin);
		
		for(int iThetaBin=1; iThetaBin<=h_matrixFull->GetXaxis()->GetNbins(); iThetaBin++) {
			double locTheta = h_matrixFull->GetXaxis()->GetBinCenter(iThetaBin);
			
			double locThrown = h_thrown->GetBinContent(h_thrown->GetXaxis()->FindBin(locEnergy), 
				h_thrown->GetYaxis()->FindBin(locTheta));
			
			for(int iRecBin=1; iRecBin<=h_matrixFull->GetYaxis()->GetNbins(); iRecBin++) {
				double locAcc = 0.0, locAccErr = 0.0;
				if(locThrown > 0.0) {
					locAcc    = h_matrixFull->GetBinContent(iThetaBin, iRecBin, iEnergyBin) / locThrown;
					locAccErr = sqrt(locThrown * locAcc * (1.0-locAcc)) / locThrown;
				}
				h_matrixFull->SetBinContent(iThetaBin, iRecBin, iEnergyBin, locAcc);
				h_matrixFull->SetBinError(iThetaBin, iRecBin, iEnergyBin, locAccErr);
			}
		}
	}
}

void YieldFitter::InitializeFluxWeights()
{
	h_fluxWeights.clear();
	
	int nHists = (int)h_dNdTheta.size();
	
	for(int ih=0; ih<nHists; ih++) {
		
		TH1F *loch1 = (TH1F*)h_fluxWeightsFull->Clone(Form("h_fluxWeights_%d\n",ih));
		
		int minEBin = h_fluxWeightsFull->GetXaxis()->FindBin(m_energyBins[ih].first  + (1.e-6));
		int maxEBin = h_fluxWeightsFull->GetXaxis()->FindBin(m_energyBins[ih].second - (1.e-6));
		
		for(int iEnergyBin=1; iEnergyBin<=h_fluxWeightsFull->GetXaxis()->GetNbins(); iEnergyBin++) {
			
			if((iEnergyBin<minEBin) || (iEnergyBin>maxEBin)) {
				loch1->SetBinContent(iEnergyBin, 0.0);
				continue;
			}
		}
		
		double integral = loch1->Integral();
		loch1->Scale(1.0/integral);
		
		h_fluxWeights.push_back(loch1);
	}
}

void YieldFitter::FitAngularYield(double minFitRange, double maxFitRange, TString outputFileName)
{
	if(LoadTheoryHists()) return;
	
	// Split h_matrixFull into appropriate energy bins:
	InitializeMatrices();
	
	// Split h_fluxWeightsFull into appropriate energy bins:
	InitializeFluxWeights();
	
	if(TVirtualFitter::GetFitter()) {
		delete TVirtualFitter::GetFitter();
		TVirtualFitter::SetFitter(nullptr);
	}
	
	/*-----------------------------*/
	// First thing to do is set up vector's storing parameters for each yield distribution:
	
	InitializeParameterArrays();
	
	int locNPars = (int)m_parameterList.size();
	vector<double> dummyPars(locNPars, 0.0);
	
	/*-----------------------------*/
	// Next, initialize the fit functions themselves:
	
	int nHists = (int)h_dNdTheta.size();
	
	for(int ih=0; ih<nHists; ih++) {
		TF1 *locFunc;
		InitializeFitFunction(&locFunc, Form("f_dNdTheta_%d",ih));
		
		// Initial Guesses for parameter values:
		locFunc->SetParameter(0, ih);
		
		locFunc->SetParameter(1, 0.515);
		locFunc->SetParameter(2, 0.5*TMath::Pi()*TMath::RadToDeg());
		locFunc->SetParameter(3, 1.0);
		locFunc->SetParameter(4, 1.0);
		if(m_components.size()==4) locFunc->SetParameter(5, 1.0);
		
		f_dNdTheta.push_back(locFunc);
	}
	
	// Some de-bugging printouts:
	printf("\n");
	printf("Fitting %d Yield Distributions:\n", nHists);
	for(int ih=0; ih<nHists; ih++) {
		printf("  energy bin %d: %.1f GeV - %.1f GeV\n", ih+1, m_energyBins[ih].first, m_energyBins[ih].second);
	}
	printf("\n");
	
	//=======================================================================================================//
	// Do the fit:
	
	CombinedChi2 locChiSq(this, m_parIndices, minFitRange, maxFitRange);
	
	ROOT::Math::Functor locFCN(locChiSq, locNPars);
	
	ROOT::Fit::Fitter locFitter;
	locFitter.Config().SetMinimizer("Minuit2", "Migrad");
	locFitter.Config().MinimizerOptions().SetPrintLevel(0);
	locFitter.Config().MinimizerOptions().SetDefaultErrorDef(0.5);
	locFitter.SetFCN(locFCN);
	locFitter.Config().SetParamsSettings(locNPars, dummyPars.data());
	
	InitializeFitParameters(locFitter);
	
	printf("Initial parameter settings:\n");
	DumpFitParameters(locFitter);
	
	//---------------------------------------------//
	// fix interference angle:
	/*
	for(int ih=0; ih<nHists; ih++) {
		int phiPar = (int)(find(m_parameterList.begin(), m_parameterList.end(), Form("#phi_{%d}[#circ]",ih)) - m_parameterList.begin());
		if(phiPar >= m_parameterList.size()) continue;
		
		locFitter.Config().ParSettings(phiPar).Set(locFitter.Config().ParSettings(phiPar).Name(), 45.5);
		locFitter.Config().ParSettings(phiPar).Fix();
	}
	*/
	//---------------------------------------------//
	
	bool ok = locFitter.FitFCN();
	
	auto result = locFitter.Result();
	result.Print(std::cout);
	
	if(m_model==SGEVORKYAN_SIGMA_VAR) {
		outputFileName = Form("fit_results_sigma_%d",m_sigmaVer);
	}
	else if(m_model==SGEVORKYAN_AP_VAR) {
		outputFileName = Form("fit_results_ap_%d",m_apVer);
	}
	else if(m_model==SGEVORKYAN_STRONG_RADIUS_VAR) {
		outputFileName = Form("fit_results_strongRadius_%s",m_strongRadiusStr.Data());
	}
	
	int nBinsTotal = 0;
	for(int ih=0; ih<nHists; ih++) {
		int locNbins = h_dNdTheta[ih]->GetXaxis()->GetNbins();
		for(int ibin=1; ibin<locNbins; ibin++) {
			double binCenter = h_dNdTheta[ih]->GetBinCenter(ibin);
			if((minFitRange<=binCenter) && (binCenter<=maxFitRange)) nBinsTotal++;
		}
	}
	int ndf = nBinsTotal - result.NPar();
	WriteOutputASCII(ndf, result, outputFileName);
	
	UpdateFitFunctions(result);
	/*
	f_dNdTheta[0]->SetParameter(1, 0.408761);
	f_dNdTheta[0]->SetParameter(2, 0.884369);
	f_dNdTheta[0]->SetParameter(3, 55.1587);
	f_dNdTheta[0]->SetParameter(4, 0.544789);
	f_dNdTheta[1]->SetParameter(1, 0.408761);
	f_dNdTheta[1]->SetParameter(2, 0.750513);
	f_dNdTheta[1]->SetParameter(3, 0.100165);
	f_dNdTheta[1]->SetParameter(4, 0.601827);
	f_dNdTheta[2]->SetParameter(1, 0.408761);
	f_dNdTheta[2]->SetParameter(2, 0.733869);
	f_dNdTheta[2]->SetParameter(3, 16.6884);
	f_dNdTheta[2]->SetParameter(4, 0.617338);
	*/
	DrawFitResult(0.0, 4.0, outputFileName);
}

void YieldFitter::DrawFitResult(double min, double max, TString fileName)
{
	int nHists = (int)h_dNdTheta.size();
	
	for(int ih=0; ih<nHists; ih++)
	{
		//-----------------------------------------//
		// Initialize canvas:
		
		TCanvas *locCanvas  = new TCanvas(Form("c_dNdTheta_%d",ih), 
			Form("%.1f GeV - %.1f GeV", m_energyBins[ih].first, m_energyBins[ih].second), 950, 700);
		
		TPad* pFit = new TPad(Form("padFit_%d",ih), Form("Fit %d",ih), 0.005, 0.3025, 0.995, 0.9950);
		pFit->SetLeftMargin(0.10);
		pFit->SetRightMargin(0.02);
		pFit->SetTopMargin(0.075);
		pFit->SetBottomMargin(0.015);
		pFit->SetTickx(); pFit->SetTicky();
		pFit->SetFrameLineWidth(2);
		pFit->SetGrid();
		
		TPad* pRes = new TPad(Form("padPull_%d",ih), Form("Pull %d",ih), 0.005, 0.0050, 0.995, 0.2975);
		pRes->SetLeftMargin(0.10);
		pRes->SetRightMargin(0.02);
		pRes->SetTopMargin(0.005);
		pRes->SetBottomMargin(0.325);
		pRes->SetTickx(); pRes->SetTicky();
		pRes->SetFrameLineWidth(2);
		pRes->SetGrid();
		
		pFit->Draw();
		pRes->Draw();
		
		//-----------------------------------------//
		// Style histogram:
		
		h_dNdTheta[ih]->GetXaxis()->SetRangeUser(min, max);
		h_dNdTheta[ih]->SetMinimum(0.0);
		
		h_dNdTheta[ih]->GetXaxis()->SetTitleSize(0.06);
		h_dNdTheta[ih]->GetXaxis()->SetLabelSize(0.05);
		h_dNdTheta[ih]->GetXaxis()->SetTitleOffset(0.9);
		h_dNdTheta[ih]->GetXaxis()->CenterTitle(true);
		
		h_dNdTheta[ih]->GetYaxis()->SetTitleSize(0.06);
		h_dNdTheta[ih]->GetYaxis()->SetLabelSize(0.05);
		h_dNdTheta[ih]->GetYaxis()->SetTitleOffset(0.8);
		h_dNdTheta[ih]->GetYaxis()->CenterTitle(true);
		
		h_dNdTheta[ih]->GetYaxis()->SetMaxDigits(2);
		
		double locMax = 0.0;
		for(int ibin=2; ibin<h_dNdTheta[ih]->GetXaxis()->GetNbins(); ibin++) {
			if(h_dNdTheta[ih]->GetXaxis()->GetBinCenter(ibin)>4.0) continue;
			double locC  = h_dNdTheta[ih]->GetBinContent(ibin);
			double locC1 = h_dNdTheta[ih]->GetBinContent(ibin-1);
			double locC2 = h_dNdTheta[ih]->GetBinContent(ibin+1);
			if(locC > (locC1+locC2)) continue;
			if(locC > locMax) locMax = locC;
		}
		h_dNdTheta[ih]->GetYaxis()->SetRangeUser(0.0, 1.2*locMax);
		
		//-----------------------------------------//
		// Create pull distribution:
		
		TH1F *h_pull = (TH1F*)h_dNdTheta[ih]->Clone(Form("yield_fit_pull_%d",ih));
		for(int ibin=1; ibin<=h_pull->GetXaxis()->GetNbins(); ibin++) {
			double y = h_dNdTheta[ih]->GetBinContent(ibin);
			double e = h_dNdTheta[ih]->GetBinError(ibin);
			double f = f_dNdTheta[ih]->Eval(h_dNdTheta[ih]->GetBinCenter(ibin));
			
			if(e!=e || e<1.0) e = sqrt(y);
			h_pull->SetBinContent(ibin, (y-f)/e);
			h_pull->SetBinError(ibin, 1.0);
		}
		h_pull->GetYaxis()->SetRangeUser(-7.0,7.0);
		
		h_pull->GetXaxis()->SetTitleSize(0.15);
		h_pull->GetXaxis()->SetLabelSize(0.10);
		h_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma}");
		h_pull->GetYaxis()->SetTitleSize(0.125);
		h_pull->GetYaxis()->SetTitleOffset(0.3);
		h_pull->GetYaxis()->SetLabelSize(0.10);
		
		//-----------------------------------------//
		// Create 'draw' functions:
		
		TF1 *fDraw, *fPrim, *fCoh, *fInc, *fInt, *fQFN;
		
		InitializeFitFunction(&fDraw, Form("f_Draw_%d",ih),         kBlack,   2, 3);
		InitializeFitFunction(&fPrim, Form("f_Primakoff_%d",ih),    kRed,     2, 2);
		InitializeFitFunction(&fCoh,  Form("f_Coherent_%d",ih),     kBlue,    2, 2);
		InitializeFitFunction_Interference(&fInt,  Form("f_Interference_%d",ih), kMagenta, 2, 2);
		if(m_components.size()==4) {
			InitializeFitFunction(&fInc, Form("f_QuasifreeProton_%d",ih),  kGreen,   2, 2);
			InitializeFitFunction(&fQFN, Form("f_QuasifreeNeutron_%d",ih), kGreen-7, 2, 2);
		}
		else {
			InitializeFitFunction(&fInc, Form("f_Incoherent_%d",ih), kGreen, 2, 2);
		}
		
		fDraw->SetParameters(f_dNdTheta[ih]->GetParameters());
		fPrim->SetParameters(ih, f_dNdTheta[ih]->GetParameter(1), f_dNdTheta[ih]->GetParameter(2), 0.0, 0.0, 0.0);
		fCoh->SetParameters(ih, 0.0, f_dNdTheta[ih]->GetParameter(2), f_dNdTheta[ih]->GetParameter(3), 0.0, 0.0);
		fInt->SetParameters(f_dNdTheta[ih]->GetParameters());
		
		fInc->SetParameters(ih, 0.0, 0.0, 0.0, f_dNdTheta[ih]->GetParameter(4), 0.0);
		if(m_components.size()>3) fQFN->SetParameters(ih, 0.0, 0.0, 0.0, 0.0, f_dNdTheta[ih]->GetParameter(5));
		
		//-----------------------------------------//
		// Draw everything:
		
		pFit->cd();
		h_dNdTheta[ih]->Draw("PE1X0");
		if(m_components.size()>3) fQFN->Draw("same");
		fInc->Draw("same");
		fInt->Draw("same");
		fCoh->Draw("same");
		fPrim->Draw("same");
		fDraw->Draw("same");
		
		//TLegend *locLeg = new TLegend(0.154, 0.657, 0.454, 0.907); // version used when just drawing yield on a single canvas (no pull)
		TLegend *locLeg = new TLegend(0.126, 0.631, 0.426, 0.881);
		locLeg->AddEntry(fDraw, "Total Fit",          "l");
		locLeg->AddEntry(fPrim, "Primakoff",          "l");
		locLeg->AddEntry(fCoh,  "Nuclear Coherent",   "l");
		locLeg->AddEntry(fInt,  "Interference",       "l");
		locLeg->AddEntry(fInc,  "Nuclear Incoherent", "l");
		locLeg->Draw();
		
		pRes->cd();
		
		h_pull->Draw("PE1X0");
		pRes->Update();
		
		TLine *l0 = new TLine(gPad->GetUxmin(),  0.0, gPad->GetUxmax(),  0.0);
		l0->SetLineColor(kBlack);
		l0->Draw("same");
		
		TLine *lp = new TLine(gPad->GetUxmin(), +2.0, gPad->GetUxmax(), +2.0);
		lp->SetLineColor(kBlack);
		lp->SetLineWidth(2);
		lp->Draw("same");
		
		TLine *lm = new TLine(gPad->GetUxmin(), -2.0, gPad->GetUxmax(), -2.0);
		lm->SetLineColor(kBlack);
		lm->SetLineWidth(2);
		lm->Draw("same");
		
		if(fileName!="") {
			TString outputFileName = Form("%s_%d.pdf", fileName.Data(), ih);
			locCanvas->SaveAs(Form("%s", outputFileName.Data()));
		}
	}
}

void YieldFitter::WriteOutputASCII(int ndf, ROOT::Fit::FitResult result, TString fileName)
{
	double chi2 = result.MinFcnValue();
	//int ndf = result.Ndf();
	
	printf("Fit Result:\n");
	printf("  Chi-squared: %f\n", chi2);
	printf("  NDF:         %d\n", ndf);
	printf("  Reduced:     %f\n", chi2/((double)ndf));
	
	if(fileName=="") return;
	
	ofstream outf(Form("%s.txt",fileName.Data()));
	char buf[256];
	
	for(int ipar=0; ipar<result.NPar(); ipar++) {
		
		int precisionVal = 5;
		if(result.Parameter(ipar)>100.0) {
			precisionVal = 3;
		}
		else if(result.Parameter(ipar)>10.0) {
			precisionVal = 4;
		}
		
		int precisionErr = 5;
		if(result.ParError(ipar)>100.0) {
			precisionErr = 3;
		}
		else if(result.ParError(ipar)>10.0) {
			precisionErr = 4;
		}
		
		outf << std::left  << std::setw(45) << result.ParName(ipar) << " " 
			 << std::right << std::setw(7)  << std::fixed << std::setprecision(precisionVal) << result.Parameter(ipar) << "   "
			 << std::right << std::setw(7)  << std::fixed << std::setprecision(precisionErr) << result.ParError(ipar) << endl;
	}
	
	int precisionVal = 5;
	if(chi2>1000.0) {
		precisionVal = 2;
	}
	else if(chi2>100.0) {
		precisionVal = 3;
	}
	else if(chi2>10.0) {
		precisionVal = 4;
	}
	
	outf << std::left  << std::setw(45) << "Chi2 / Ndf" << " " 
		 << std::right << std::setw(7)  << std::fixed << std::setprecision(precisionVal) << chi2 << "   "
		 << std::right << std::setw(7)  << ndf << endl;
	outf.close();
}

void YieldFitter::InitializeFitFunction(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	int nParameters = 1 + m_components.size() + 1;
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldFitFunction, m_minReconAngle, m_maxReconAngle, nParameters);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "HistIndex");
	
	(*f1)->SetParName(1, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(2, "#phi[#circ]");
	(*f1)->SetParName(3, "A_{Coh}");
	if(m_components.size()==4) {
		(*f1)->SetParName(4, "A_{QFP}");
		(*f1)->SetParName(5, "A_{QFN}");
	}
	else {
		(*f1)->SetParName(4, "A_{Inc}");
	}
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

void YieldFitter::InitializeFitFunction_Interference(TF1 **f1, TString funcName, int lineColor, int lineStyle, int lineWidth)
{
	// initialize fit function for each angular bin:
	
	int nParameters = 1 + m_components.size() + 1;
	
	*f1 = new TF1(funcName.Data(), this, &YieldFitter::YieldFitFunction_Interference, m_minReconAngle, m_maxReconAngle, nParameters);
	
	// set names for each parameter:
	
	(*f1)->SetParName(0, "HistIndex");
	
	(*f1)->SetParName(1, "#Gamma(#eta#rightarrow#gamma#gamma)[keV]");
	(*f1)->SetParName(2, "#phi[#circ]");
	(*f1)->SetParName(3, "A_{Coh}");
	if(m_components.size()==4) {
		(*f1)->SetParName(4, "A_{QFP}");
		(*f1)->SetParName(5, "A_{QFN}");
	}
	else {
		(*f1)->SetParName(4, "A_{Inc}");
	}
	
	(*f1)->SetLineColor(lineColor);
	(*f1)->SetLineStyle(lineStyle);
	(*f1)->SetLineWidth(lineWidth);
	
	return;
}

int YieldFitter::LoadTheoryHists()
{
	printf("\nREADING THEORY CALCULATIONS...\n");
	printf("  Model: %s\n", GetModelString().Data());
	
	int nComponents = InitializeModelComponents();
	
	vector<TString> theoryFileNames;
	for(int i=0; i<nComponents; i++) {
		theoryFileNames.push_back("");
	}
	
	switch(m_model) {
		case UNKNOWN:
		{
			printf("UNKNOWN ModelType in YieldFitter Object.\n");
			return 1;
		}
		case AFIX:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			break;
		}
		case SGEVORKYAN:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			break;
		}
		case MIXED_V1:
		{
			for(int i=1; i<2; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			for(int i=2; i<4; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			break;
		}
		case MIXED_V2:
		{
			theoryFileNames[0] = "/work/halld/home/ijaegle/afix_calculation/root-files/he4-eta-xs-theory-AFix-v5.root";
			for(int i=1; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			break;
		}
		case SGEVORKYAN_FERMI:
		{
			for(int i=0; i<(nComponents-1); i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-theory-SGevorkyan-v4.root";
			theoryFileNames[nComponents-1] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/test-sergey-gevorkyanb.root";
			break;
		}
		case SGEVORKYAN_UPD_V0:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v0.root";
			break;
		}
		case SGEVORKYAN_UPD_V1:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v1.root";
			break;
		}
		case SGEVORKYAN_UPD_V2:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root";
			break;
		}
		case SGEVORKYAN_UPD_V3:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v3.root";
			break;
		}
		case SGEVORKYAN_UPD_FERMI:
		{
			for(int i=0; i<nComponents; i++) 
				theoryFileNames[i] = "/work/halld/home/ijaegle/sgevorkyan_calculation/root-files/he4-eta-xs-sgevorkyan-v2-folded.root";
			break;
		}
		case SGEVORKYAN_SIGMA_VAR:
		{
			for(int i=0; i<nComponents; i++) {
				theoryFileNames[i] = Form(
					"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/sigma_variations/he4-eta-xs-sgevorkyan-v%d.root", 
					m_sigmaVer);
			}
			break;
		}
		case SGEVORKYAN_AP_VAR:
		{
			for(int i=0; i<(nComponents-1); i++) {
				theoryFileNames[i] = Form(
					"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/ap_variations/he4-eta-xs-sgevorkyan-v%d.root", 
					m_apVer);
			}
			theoryFileNames[nComponents-1] = 
				"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root";
			break;
		}
		case SGEVORKYAN_STRONG_RADIUS_VAR:
		{
			theoryFileNames[0] = 
				"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root";
			
			theoryFileNames[1] = Form(
				"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/radius_variations/he4-eta-xs-sgevorkyan-%s.root", 
				m_strongRadiusStr.Data());
			
			theoryFileNames[2] = 
				"/work/halld/home/andrsmit/primex_eta_analysis/theory/sgevorkyan/farm/rootFiles/he4-eta-xs-sgevorkyan-v2.root";
			break;
		}
	}
	
	for(int i=0; i<nComponents; i++) {
		
		printf("  component %d: %s\n", i, theoryFileNames[i].Data());
		
		TFile *fTheory  = new TFile(theoryFileNames[i].Data(), "READ");
		TH2F *hTheory2D = (TH2F*)fTheory->Get(Form("%s",m_components[i].Data()));
		
		/*
		NOTE:
		I'm hard-coding the 2-d theory histogram (theta vs. energy) from E=7 to E=12 GeV
		Need to place multiple checks in the code to ensure that I don't actually set the cross section for something to 0
		when it shouldn't be. For example, if later on, I input a histogram for the yield at 7 GeV. 
		*/
		
		if(m_beamEnergyBinSize != 0.05) {
			printf("BEAM ENERGY BINNING IS WRONG. QUITING.\n");
			exit(1);
		}
		
		double locMinEnergy =  7.0;
		double locMaxEnergy = 12.0;
		int locNEnergyBins = (locMaxEnergy-locMinEnergy)/m_beamEnergyBinSize;
		
		if(m_thrownAngleBinSize != 0.01) {
			printf("THROWN ANGULAR BINNING IS WRONG. QUITING.\n");
			exit(1);
		}
		
		double locMinTheta = 0.0;
		double locMaxTheta = 5.0;
		int locNThetaBins = (locMaxTheta-locMinTheta)/m_thrownAngleBinSize;
		
		TH2F *hTheory_rebinned = new TH2F(Form("h_theory_%s",m_components[i].Data()), "",
			locNEnergyBins, locMinEnergy, locMaxEnergy,
			locNThetaBins,  locMinTheta,  locMaxTheta);
		hTheory_rebinned->SetDirectory(0);
		
		double locBeamEnergy = locMinEnergy + 0.5*m_beamEnergyBinSize;
		while(locBeamEnergy < locMaxEnergy)
		{
			int locMinEnergyBin = hTheory2D->GetXaxis()->FindBin(locBeamEnergy - 0.5*m_beamEnergyBinSize + (1.e-6));
			int locMaxEnergyBin = hTheory2D->GetXaxis()->FindBin(locBeamEnergy + 0.5*m_beamEnergyBinSize - (1.e-6));
			
			double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin+1.0);
			
			TH1F *loch1 = (TH1F*)hTheory2D->ProjectionY(Form("%s_%f_%f", m_components[i].Data(), 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin);
			loch1->Scale(1.0/nEnergyBins);
			
			// rebin histograms to match acceptance matrices:
			
			double theoryAngleBinSize = loch1->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				loch1->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				loch1->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			
			// check that the binning of loch1 overlaps with hTheory_rebinned:
			
			if(fabs(loch1->GetBinWidth(1) - hTheory_rebinned->GetYaxis()->GetBinWidth(1)) > 1.e-6) {
				printf("ISSUE WITH THROWN ANGLE BINNING IN THEORY HISTS (flag 1). QUITING.\n");
				printf("  loch1->GetBinWidth(1) = %f; h_Theory_rebinned->GetYaxis()->GetBinWidth(1) = %f\n", 
					loch1->GetBinWidth(1), hTheory_rebinned->GetYaxis()->GetBinWidth(1));
				exit(1);
			}
			if(fabs(loch1->GetBinCenter(1) - hTheory_rebinned->GetYaxis()->GetBinCenter(1)) > 1.e-6) {
				printf("ISSUE WITH THROWN ANGLE BINNING IN THEORY HISTS (flag 2). QUITING.\n");
				exit(1);
			}
			
			for(int thBin=1; thBin<hTheory_rebinned->GetYaxis()->GetNbins(); thBin++) {
				hTheory_rebinned->SetBinContent(hTheory_rebinned->GetXaxis()->FindBin(locBeamEnergy), thBin,
					loch1->GetBinContent(thBin));
			}
			
			delete loch1;
			locBeamEnergy += m_beamEnergyBinSize;
		}
		
		h_Theory.push_back(hTheory_rebinned);
		
		fTheory->Close();
	}
	
	if(m_model>=6) {
		
		TFile *fTheory = new TFile(theoryFileNames[0].Data(), "READ");
		
		TH2F *hPrimReal2D = (TH2F*)fTheory->Get("amp_prim_real_vs_egam");
		TH2F *hPrimImag2D = (TH2F*)fTheory->Get("amp_prim_imag_vs_egam");
		TH2F *hStrongReal2D = (TH2F*)fTheory->Get("amp_coh_real_vs_egam");
		TH2F *hStrongImag2D = (TH2F*)fTheory->Get("amp_coh_imag_vs_egam");
		
		/*
		NOTE:
		I'm hard-coding the 2-d theory histogram (theta vs. energy) from E=7 to E=12 GeV
		Need to place multiple checks in the code to ensure that I don't actually set the cross section for something to 0
		when it shouldn't be. For example, if later on, I input a histogram for the yield at 7 GeV. 
		*/
		
		double locMinEnergy =  7.0;
		double locMaxEnergy = 12.0;
		int locNEnergyBins = (locMaxEnergy-locMinEnergy)/m_beamEnergyBinSize;
		
		double locMinTheta = 0.0;
		double locMaxTheta = 5.0;
		int locNThetaBins = (locMaxTheta-locMinTheta)/m_thrownAngleBinSize;
		
		h_PrimReal = new TH2F("h_PrimReal", "Real Part of Primakoff Amplitude",
			locNEnergyBins, locMinEnergy, locMaxEnergy,
			locNThetaBins,  locMinTheta,  locMaxTheta);
		h_PrimReal->SetDirectory(0);
		
		h_PrimImag = new TH2F("h_PrimImag", "Imaginary Part of Primakoff Amplitude",
			locNEnergyBins, locMinEnergy, locMaxEnergy,
			locNThetaBins,  locMinTheta,  locMaxTheta);
		h_PrimImag->SetDirectory(0);
		
		h_StrongReal = new TH2F("h_StrongReal", "Real Part of Nuclear Coherent Amplitude",
			locNEnergyBins, locMinEnergy, locMaxEnergy,
			locNThetaBins,  locMinTheta,  locMaxTheta);
		h_StrongReal->SetDirectory(0);
		
		h_StrongImag = new TH2F("h_StrongImag", "Imaginary Part of Nuclear Coherent Amplitude",
			locNEnergyBins, locMinEnergy, locMaxEnergy,
			locNThetaBins,  locMinTheta,  locMaxTheta);
		h_StrongImag->SetDirectory(0);
		
		double locBeamEnergy = locMinEnergy + 0.5*m_beamEnergyBinSize;
		while(locBeamEnergy < locMaxEnergy)
		{
			int locMinEnergyBin = hPrimReal2D->GetXaxis()->FindBin(locBeamEnergy - 0.5*m_beamEnergyBinSize + (1.e-6));
			int locMaxEnergyBin = hPrimReal2D->GetXaxis()->FindBin(locBeamEnergy + 0.5*m_beamEnergyBinSize - (1.e-6));
			
			double nEnergyBins = (double)(locMaxEnergyBin-locMinEnergyBin+1.0);
			
			TH1F *locPrimReal = (TH1F*)hPrimReal2D->ProjectionY(Form("prim_real_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin);
			locPrimReal->Scale(1.0/nEnergyBins);
			
			TH1F *locPrimImag = (TH1F*)hPrimImag2D->ProjectionY(Form("prim_imag_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin);
			locPrimImag->Scale(1.0/nEnergyBins);
			
			TH1F *locStrongReal = (TH1F*)hStrongReal2D->ProjectionY(Form("coh_real_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin);
			locStrongReal->Scale(1.0/nEnergyBins);
			
			TH1F *locStrongImag = (TH1F*)hStrongImag2D->ProjectionY(Form("coh_imag_%f_%f", 
				locBeamEnergy, locBeamEnergy+m_beamEnergyBinSize), locMinEnergyBin, locMaxEnergyBin);
			locStrongImag->Scale(1.0/nEnergyBins);
			
			
			// rebin histograms to match acceptance matrices:
			
			double theoryAngleBinSize;
			
			theoryAngleBinSize = locPrimReal->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locPrimReal->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locPrimReal->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			
			// check that the binning of locPrimReal overlaps with h_PrimReal:
			
			if(fabs(locPrimReal->GetBinWidth(1) - h_PrimReal->GetYaxis()->GetBinWidth(1)) > 1.e-6) {
				printf("ISSUE WITH THROWN ANGLE BINNING IN THEORY HISTS (flag 3). QUITING.\n");
				exit(1);
			}
			if(fabs(locPrimReal->GetBinCenter(1) - h_PrimReal->GetYaxis()->GetBinCenter(1)) > 1.e-6) {
				printf("ISSUE WITH THROWN ANGLE BINNING IN THEORY HISTS (flag 4). QUITING.\n");
				exit(1);
			}
			
			for(int thBin=1; thBin<h_PrimReal->GetYaxis()->GetNbins(); thBin++) {
				h_PrimReal->SetBinContent(h_PrimReal->GetXaxis()->FindBin(locBeamEnergy), thBin,
					locPrimReal->GetBinContent(thBin));
			}
			delete locPrimReal;
			
			
			theoryAngleBinSize = locPrimImag->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locPrimImag->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locPrimImag->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			for(int thBin=1; thBin<h_PrimImag->GetYaxis()->GetNbins(); thBin++) {
				h_PrimImag->SetBinContent(h_PrimImag->GetXaxis()->FindBin(locBeamEnergy), thBin,
					locPrimImag->GetBinContent(thBin));
			}
			delete locPrimImag;
			
			
			theoryAngleBinSize = locStrongReal->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locStrongReal->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locStrongReal->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			for(int thBin=1; thBin<h_StrongReal->GetYaxis()->GetNbins(); thBin++) {
				h_StrongReal->SetBinContent(h_StrongReal->GetXaxis()->FindBin(locBeamEnergy), thBin,
					locStrongReal->GetBinContent(thBin));
			}
			delete locStrongReal;
			
			
			theoryAngleBinSize = locStrongImag->GetXaxis()->GetBinWidth(1);
			if(fabs(theoryAngleBinSize-m_thrownAngleBinSize) > 1.e-6) {
				locStrongImag->Rebin((int)(m_thrownAngleBinSize/theoryAngleBinSize));
				locStrongImag->Scale(theoryAngleBinSize/m_thrownAngleBinSize);
			}
			for(int thBin=1; thBin<h_StrongImag->GetYaxis()->GetNbins(); thBin++) {
				h_StrongImag->SetBinContent(h_StrongImag->GetXaxis()->FindBin(locBeamEnergy), thBin,
					locStrongImag->GetBinContent(thBin));
			}
			delete locStrongImag;
			
			locBeamEnergy += m_beamEnergyBinSize;
		}
		
		delete hPrimReal2D;
		delete hPrimImag2D;
		delete hStrongReal2D;
		delete hStrongImag2D;
		fTheory->Close();
	}
	
	/*
	h_TheoryTulio = new TH1F("TulioInc", "; #theta_{#eta} [#circ]; d#sigma/d#theta_{#eta} [#mub/rad]", 50, 0.0, 10.0);
	
	ifstream inf("/work/halld/home/andrsmit/primex_eta_analysis/theory/tulio/test_180_0.90.txt");
	double aa, bb, cc;
	for(int iline=0; iline<50; iline++) {
		inf >> aa >> bb >> cc;
		int locBin = h_TheoryTulio->FindBin(aa);
		h_TheoryTulio->SetBinContent(locBin, cc);
	}
	inf.close();
	
	// Get Spline fit to Tulio's distribution:
	spline = new TSpline3(h_TheoryTulio);
	f_TheoryTulio = new TF1("f_TheoryTulio", [=](double * x, double * p) { 
		return p[0] * spline->Eval(x[0]);
	}, 0, 1, 1);
	
	f_TheoryTulio->SetParameter(0, 1);
	h_TheoryTulio->Fit("f_TheoryTulio");
	*/
	return 0;
}

int YieldFitter::InitializeModelComponents() 
{
	int nComponents = 0;
	switch(m_model) {
		case AFIX:
			nComponents = 4;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_qfp_vs_egam");
			m_components.push_back("xs_qfn_vs_egam");
			break;
		case SGEVORKYAN:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case MIXED_V1:
			nComponents = 4;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_qfp_vs_egam");
			m_components.push_back("xs_qfn_vs_egam");
			break;
		case MIXED_V2:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_FERMI:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("qf_txs_lab");
			break;
		case SGEVORKYAN_UPD_V0:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V1:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V2:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_V3:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_UPD_FERMI:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("qf_txs_lab");
			break;
		case SGEVORKYAN_SIGMA_VAR:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_AP_VAR:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
		case SGEVORKYAN_STRONG_RADIUS_VAR:
			nComponents = 3;
			m_components.push_back("xs_prim_vs_egam");
			m_components.push_back("xs_coh_vs_egam");
			m_components.push_back("xs_inc_vs_egam");
			break;
	}
	
	return nComponents;
}


double YieldFitter::GetExpectedYield(double minAngle, double maxAngle, int histIndex, double *par, int isolateInter)
{
	int locMinReconAngleBin = h_matrixFull->GetYaxis()->FindBin(minAngle+(1.e-6));
	int locMaxReconAngleBin = h_matrixFull->GetYaxis()->FindBin(maxAngle-(1.e-6));
	
	double angleBinSize = m_thrownAngleBinSize * TMath::DegToRad();
	
	double Gamma = par[0];
	double Phi   = par[1];
	double Acoh  = par[2];
	double AincP = par[3];
	double AincN = (m_components.size()>3) ? par[4] : 0.0;
	/*
	printf("    Gamma: %f\n", Gamma);
	printf("    Acoh:  %f\n", Acoh);
	printf("    Phi:   %f\n", Phi);
	printf("    Ainc:  %f\n", AincP);
	*/
	int locMinEnergyBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[histIndex].first  + 0.5*m_beamEnergyBinSize);
	int locMaxEnergyBin = h_matrixFull->GetZaxis()->FindBin(m_energyBins[histIndex].second - 0.5*m_beamEnergyBinSize);
	
	double dNdTheta = 0.0;
	for(int iReconBin=locMinReconAngleBin; iReconBin<=locMaxReconAngleBin; iReconBin++) {
		
		for(int iEnergyBin=locMinEnergyBin; iEnergyBin<=locMaxEnergyBin; iEnergyBin++) {
			double locEnergy = h_matrixFull->GetZaxis()->GetBinCenter(iEnergyBin);
			
			for(int iThetaBin=1; iThetaBin<=h_matrixFull->GetXaxis()->GetNbins(); iThetaBin++) {
				double locMatrix = h_matrices[histIndex]->GetBinContent(iThetaBin, iReconBin, iEnergyBin);
				double locCS;
				if(isolateInter) locCS = GetCrossSectionInterference(locEnergy, iThetaBin, Gamma, Acoh, Phi);
				else             locCS = GetCrossSection(locEnergy, iThetaBin, Gamma, Acoh, AincP, AincN, Phi);
				
				dNdTheta += (locMatrix * locCS * h_fluxWeights[histIndex]->GetBinContent(h_fluxWeights[histIndex]->FindBin(locEnergy)));
			}
		}
	}
	
	int locMinFluxBin = h_fluxWeightsFull->GetXaxis()->FindBin(m_energyBins[histIndex].first  + 0.5*m_beamEnergyBinSize);
	int locMaxFluxBin = h_fluxWeightsFull->GetXaxis()->FindBin(m_energyBins[histIndex].second - 0.5*m_beamEnergyBinSize);
	
	double fractionalLumi = m_luminosity * h_fluxWeightsFull->Integral(locMinFluxBin, locMaxFluxBin);
	//printf("  fracitonalLumi: %f\n", fractionalLumi/m_luminosity);
	
	double yield = dNdTheta * fractionalLumi * EtaAnalyzer::m_branchingRatio * (m_thrownAngleBinSize*TMath::DegToRad());
	return yield;
}

double YieldFitter::YieldFitFunction(double *x, double *par)
{
	int HistIndex = (int)(par[0]+1.e-6);
	
	int reconAngleBin = h_matrixFull->GetYaxis()->FindBin(x[0]);
	double reconAngleBinWidth = h_matrixFull->GetYaxis()->GetBinWidth(reconAngleBin);
	
	double minAngle = h_matrixFull->GetYaxis()->GetBinCenter(reconAngleBin) - 0.5*reconAngleBinWidth;
	double maxAngle = h_matrixFull->GetYaxis()->GetBinCenter(reconAngleBin) + 0.5*reconAngleBinWidth;
	
	vector<double> locPars;
	locPars.push_back(par[1]); // Gamma
	locPars.push_back(par[2]); // phi
	locPars.push_back(par[3]); // A_coh
	locPars.push_back(par[4]); // A_inc
	locPars.push_back(par[5]);
	
	double yield = GetExpectedYield(minAngle, maxAngle, HistIndex, locPars.data());
	
	//printf("expected yield between %f-%f deg: %f\n", minAngle, maxAngle, yield);
	//printf(" correction needed: %f\n", h_dNdTheta[HistIndex]->GetBinWidth(1) / reconAngleBinWidth);
	
	// we need to apply a correction based on differences in bin widths:
	
	yield *= (h_dNdTheta[HistIndex]->GetBinWidth(1) / reconAngleBinWidth);
	
	return yield;
}

double YieldFitter::YieldFitFunction_Interference(double *x, double *par)
{
	int HistIndex = (int)(par[0]+1.e-6);
	
	int reconAngleBin = h_matrixFull->GetYaxis()->FindBin(x[0]);
	double reconAngleBinWidth = h_matrixFull->GetYaxis()->GetBinWidth(reconAngleBin);
	
	double minAngle = h_matrixFull->GetYaxis()->GetBinCenter(reconAngleBin) - 0.5*reconAngleBinWidth;
	double maxAngle = h_matrixFull->GetYaxis()->GetBinCenter(reconAngleBin) + 0.5*reconAngleBinWidth;
	
	vector<double> locPars;
	locPars.push_back(par[1]); // Gamma
	locPars.push_back(par[2]); // phi
	locPars.push_back(par[3]); // A_coh
	locPars.push_back(par[4]); // A_inc
	locPars.push_back(par[5]);
	
	double yield = GetExpectedYield(minAngle, maxAngle, HistIndex, locPars.data(), 1);
	
	// we need to apply a correction based on differences in bin widths:
	
	yield *= (h_dNdTheta[HistIndex]->GetBinWidth(1) / reconAngleBinWidth);
	
	return yield;
}

double YieldFitter::GetCrossSection(double beamEnergy, int thetaBin, 
	double Gamma, double Acoh, double AincP, double AincN, double Phi) 
{
	double gammaGen = (m_model>=6) ? 0.515 : 0.510;
	
	int locEnergyBin = h_Theory[0]->GetXaxis()->FindBin(beamEnergy);
	if((locEnergyBin < 1) || (locEnergyBin > h_Theory[0]->GetXaxis()->GetNbins())) {
		printf("BEAM ENERGY BIN OUTSIDE OF HISTOGRAM LIMITS (flag 1). QUITING.\n");
		exit(1);
	}
	
	double locPrim = (Gamma/gammaGen) * h_Theory[0]->GetBinContent(locEnergyBin, thetaBin);
	double locCoh  = Acoh * h_Theory[1]->GetBinContent(locEnergyBin, thetaBin);
	double locInc  = AincP * h_Theory[2]->GetBinContent(locEnergyBin, thetaBin);
	if(m_components.size()>3) {
		locInc += (AincN * h_Theory[3]->GetBinContent(locEnergyBin, thetaBin));
	}
	
	// use tulio's incoherent:
	
	//double locTheta = h_Theory[0][energyBin]->GetXaxis()->GetBinCenter(thetaBin);
	//locInc = AincP * f_TheoryTulio->Eval(locTheta);
	
	// TEST:
	//locCoh /= pow(beamEnergy,0.2);
	
	double locInt  = GetCrossSectionInterference(beamEnergy, thetaBin, Gamma, Acoh, Phi);
	return (locPrim+locCoh+locInc+locInt);
}

double YieldFitter::GetCrossSectionInterference(double beamEnergy, int thetaBin, double Gamma, double Acoh, double Phi) 
{
	double gammaGen = (m_model>=6) ? 0.515 : 0.510;
	
	double locInt = 0.0;
	if(m_model>=6) {
		
		int locEnergyBin = h_PrimReal->GetXaxis()->FindBin(beamEnergy);
		if((locEnergyBin < 1) || (locEnergyBin > h_PrimReal->GetXaxis()->GetNbins())) {
			printf("BEAM ENERGY BIN OUTSIDE OF HISTOGRAM LIMITS (flag 2). QUITING.\n");
			exit(1);
		}
		
		double locPrimReal = h_PrimReal->GetBinContent(locEnergyBin, thetaBin);
		double locPrimImag = h_PrimImag->GetBinContent(locEnergyBin, thetaBin);
		double locStrongReal = h_StrongReal->GetBinContent(locEnergyBin, thetaBin);
		double locStrongImag = h_StrongImag->GetBinContent(locEnergyBin, thetaBin);
		
		// TEST:
		//locStrongReal /= pow(beamEnergy,0.2);
		//locStrongImag /= pow(beamEnergy,0.2);
		
		double int1 = (locPrimReal*locStrongReal + locPrimImag*locStrongImag) * cos(Phi*TMath::DegToRad());
		double int2 = (locPrimImag*locStrongReal - locPrimReal*locStrongImag) * sin(Phi*TMath::DegToRad());
		locInt = 2.0*sqrt(Gamma/gammaGen)*sqrt(Acoh)*(int1+int2);
	} else {
		double locPrim = (Gamma/gammaGen) * h_Theory[0]->GetBinContent(h_Theory[0]->GetXaxis()->FindBin(beamEnergy), thetaBin);
		double locCoh  = Acoh * h_Theory[1]->GetBinContent(h_Theory[1]->GetXaxis()->FindBin(beamEnergy), thetaBin);
		
		// TEST:
		//locCoh /= pow(beamEnergy,0.4);
		
		locInt  = 2.0*sqrt(locPrim*locCoh)*cos(Phi*TMath::DegToRad());
	}
	return locInt;
}
