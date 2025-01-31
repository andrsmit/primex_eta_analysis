#include "EtaAnalyzer.h"

TString fluxDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/rootFiles";

void EtaAnalyzer::CalcLuminosity() {
	
	m_fluxWeights.clear();
	
	double targetDensity = 0.1217;    // g/cm3
	double targetLength  = 29.535;    // cm
	double targetMass    = 4.002602;  // g/mol
	double avogdroNum    = 6.02214e23;
	
	double targetThickness = targetDensity * targetLength * (1.0/targetMass) * avogdroNum * (1.e-30); // atoms/ub
		
	TString fluxFileName = Form("%s/phase%d/full.root", fluxDirectory.Data(), m_phase);
	if(gSystem->AccessPathName(fluxFileName.Data())) {
		cout << "Flux filename doesn't exist." << endl;
		return;
	}
	
	// read in flux histogram from file:
	
	TFile *locFluxFile = new TFile(fluxFileName.Data(), "READ");
	TH1F  *locFluxHist = (TH1F*)locFluxFile->Get("flux_vs_egamma");
	
	double integratedFlux = IntegrateFluxHist(locFluxHist);
	m_luminosity = integratedFlux * targetThickness;
	
	//------------------------------------------------//
	// Store fraction of flux in each energy bin in m_fluxWeights vector:
	
	double fluxBinSize = locFluxHist->GetXaxis()->GetBinWidth(1); // assumes fixed bin size
	double locBeamE = minBeamEnergy + 0.5*beamEnergyBinSize;
	
	while(locBeamE < maxBeamEnergy) {
		
		double locMinBeamE = locBeamE - 0.5*beamEnergyBinSize;
		double locMaxBeamE = locBeamE + 0.5*beamEnergyBinSize;
		double locFlux = locFluxHist->Integral(locFluxHist->GetXaxis()->FindBin(locMinBeamE+0.5*fluxBinSize),
			locFluxHist->GetXaxis()->FindBin(locMaxBeamE-0.5*fluxBinSize));
		m_fluxWeights.push_back(locFlux);
		locBeamE += beamEnergyBinSize;
	}
	for(int ibin=0; ibin<m_fluxWeights.size(); ibin++) {
		m_fluxWeights[ibin] = m_fluxWeights[ibin] / integratedFlux;
		//cout << m_fluxWeights[ibin] << endl;
	}
	//------------------------------------------------//
	
	locFluxFile->Close();
	
	printf("Integrated Luminosity [%.2f GeV - %.2f GeV]: %f pb-1\n", minBeamEnergy, maxBeamEnergy, (1.e-6)*m_luminosity);
	
	return;
}

void EtaAnalyzer::CalcEmptyTargetFluxRatio() {
	
	//
	// If we assume the beamline background we want to subtract originates from downstream the target, then
	// we need to calculate the ratio of photon flux at the downstream exit window of the target cell.
	// For the 3.8% RL Helium target, roughly 3% of photons are absorbed within the target. 
	// The ratio we use to scale our empty target data should be reduced by this factor.
	//
	double scaleFactor = 0.97;
	
	TString fluxFileNameFull  = Form("%s/phase%d/full.root",  fluxDirectory.Data(), m_phase);
	TString fluxFileNameEmpty = Form("%s/phase%d/empty.root", fluxDirectory.Data(), m_phase);
	
	if(gSystem->AccessPathName(fluxFileNameFull.Data()) || gSystem->AccessPathName(fluxFileNameEmpty.Data())) {
		cout << "Flux filenames don't exist." << endl;
		return;
	}
	
	// read in flux histograms from files:
	
	TFile *fFull  = new TFile(fluxFileNameFull.Data(),  "READ");
	TFile *fEmpty = new TFile(fluxFileNameEmpty.Data(), "READ");
	
	TH1F  *hFull  = (TH1F*)fFull->Get("flux_vs_egamma")->Clone("hFull");
	TH1F  *hEmpty = (TH1F*)fEmpty->Get("flux_vs_egamma")->Clone("hEmpty");
	
	double integratedFluxFull  = IntegrateFluxHist(hFull);
	double integratedFluxEmpty = IntegrateFluxHist(hEmpty);
	
	m_emptyTargetFluxRatio = scaleFactor * (integratedFluxFull / integratedFluxEmpty);
	
	printf("Photon Flux Ratio (Full/Empty): %f\n", m_emptyTargetFluxRatio);
	
	return;
}

double EtaAnalyzer::IntegrateFluxHist(TH1F *hFlux) {
	
	// find the bin numbers to integrate flux over:
	
	double locBinSize = hFlux->GetXaxis()->GetBinWidth(1);
	int locMinimumBin = hFlux->GetXaxis()->FindBin(minBeamEnergy + 0.5*locBinSize);
	int locMaximumBin = hFlux->GetXaxis()->FindBin(maxBeamEnergy - 0.5*locBinSize);
	
	// check that the bin edges exactly correspond to specified energy range:
	
	double locMinimumEnergy = hFlux->GetXaxis()->GetBinCenter(locMinimumBin) - 0.5*locBinSize;
	double locMaximumEnergy = hFlux->GetXaxis()->GetBinCenter(locMaximumBin) + 0.5*locBinSize;
	
	if((fabs(locMinimumEnergy-minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-maxBeamEnergy)>1.e6)) {
		printf("\n\nWarning! Bin edges of flux histogram are inconsistent with specified energy range:\n");
		printf("   Specified energy range: %.5f GeV - %.5f GeV\n", minBeamEnergy-locMinimumEnergy, maxBeamEnergy-locMaximumEnergy);
		printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
	}
	
	// integrate photon flux:
	
	double integratedFlux = hFlux->Integral(locMinimumBin, locMaximumBin);
	return integratedFlux;
}
