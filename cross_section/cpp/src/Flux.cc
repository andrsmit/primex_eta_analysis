#include "EtaAnalyzer.h"

TString fluxDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/rootFiles";

int EtaAnalyzer::LoadLuminosity() {
	
	printf("\nREADING PHOTON FLUX...\n");
	
	double targetThickness = m_targetDensity * m_targetLength * (1.0/m_targetMass) * m_avogdroNum * (1.e-30); // atoms/ub
	
	TString fluxFileName = Form("%s/phase%d/full.root", fluxDirectory.Data(), m_phase);
	if(gSystem->AccessPathName(fluxFileName.Data())) {
		cout << "Flux filename doesn't exist." << endl;
		return 1;
	}
	
	// read in flux histogram from file:
	
	TFile *locFluxFile = new TFile(fluxFileName.Data(), "READ");
	TH1F  *locFluxHist = (TH1F*)locFluxFile->Get("flux_vs_egamma");
	
	double integratedFlux = IntegrateFluxHist(locFluxHist, m_minBeamEnergy, m_maxBeamEnergy);
	m_luminosity = integratedFlux * targetThickness;
	
	//
	// To account for the residual-gas contribution when subtracting empty background:
	// Previously I only applied this correction factor if the empty target bkgd was subtracted 
	// prior to fitting the invariant mass spectrum. However, there is no way to distinguish etas
	// that were produced off the gas and etas produced off the target walls.
	// So even when fitting and obtaining the pdf of the empty bkgd, I still subtract the full peak
	// around the eta mass region and apply the correction to the target density accordingly.
	// 
	m_luminosity *= 0.9817;
	
	//------------------------------------------------//
	// Store fraction of flux in each energy bin in h_fluxWeights:
	
	InitializeFluxHist();
	
	for(int ibin=1; ibin<=h_fluxWeights->GetXaxis()->GetNbins(); ibin++) {
		double locMinEnergy = h_fluxWeights->GetBinCenter(ibin) - 0.5*m_beamEnergyBinSize;
		double locMaxEnergy = h_fluxWeights->GetBinCenter(ibin) + 0.5*m_beamEnergyBinSize;
		
		double locFlux = IntegrateFluxHist(locFluxHist, locMinEnergy, locMaxEnergy);
		if(integratedFlux>0.0) {
			h_fluxWeights->SetBinContent(ibin, locFlux/integratedFlux);
			h_fluxWeights->SetBinError(ibin, sqrt(locFlux)/integratedFlux);
		}
	}
	
	//------------------------------------------------//
	
	locFluxFile->Close();
	printf("  Integrated Photon Flux: %f\n", integratedFlux);
	printf("  Integrated Luminosity [%.2f GeV - %.2f GeV]: %f pb-1\n", m_minBeamEnergy, m_maxBeamEnergy, (1.e-6)*m_luminosity);
	
	return 0;
}

int EtaAnalyzer::LoadEmptyTargetFluxRatio() {
	
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
		return 1;
	}
	
	// read in flux histograms from files:
	
	TFile *fFull  = new TFile(fluxFileNameFull.Data(),  "READ");
	TFile *fEmpty = new TFile(fluxFileNameEmpty.Data(), "READ");
	
	TH1F  *hFull  = (TH1F*)fFull->Get("flux_vs_egamma")->Clone("hFull");
	TH1F  *hEmpty = (TH1F*)fEmpty->Get("flux_vs_egamma")->Clone("hEmpty");
	
	double integratedFluxFull  = IntegrateFluxHist(hFull,  m_minBeamEnergy, m_maxBeamEnergy);
	double integratedFluxEmpty = IntegrateFluxHist(hEmpty, m_minBeamEnergy, m_maxBeamEnergy);
	
	m_emptyTargetFluxRatio = scaleFactor * (integratedFluxFull / integratedFluxEmpty);
	
	printf("  Photon Flux Ratio (Full/Empty): %f\n", m_emptyTargetFluxRatio);
	
	return 0;
}

double EtaAnalyzer::IntegrateFluxHist(TH1F *hFlux, double minEnergy, double maxEnergy) {
	
	// find the bin numbers to integrate flux over:
	
	double locBinSize = hFlux->GetXaxis()->GetBinWidth(1);
	int locMinimumBin = hFlux->GetXaxis()->FindBin(minEnergy + 0.5*locBinSize);
	int locMaximumBin = hFlux->GetXaxis()->FindBin(maxEnergy - 0.5*locBinSize);
	
	// check that the bin edges exactly correspond to specified energy range:
	
	double locMinimumEnergy = hFlux->GetXaxis()->GetBinCenter(locMinimumBin) - 0.5*locBinSize;
	double locMaximumEnergy = hFlux->GetXaxis()->GetBinCenter(locMaximumBin) + 0.5*locBinSize;
	
	if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
		printf("\n\nWarning! Bin edges of flux histogram are inconsistent with specified energy range:\n");
		printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
		printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
	}
	
	// integrate photon flux:
	
	double integratedFlux = hFlux->Integral(locMinimumBin, locMaximumBin);
	return integratedFlux;
}

void EtaAnalyzer::InitializeFluxHist() {
	
	double locNBins = (int)((m_maxBeamEnergy-m_minBeamEnergy)/m_beamEnergyBinSize);
	
	// consistency check:
	double locMin = m_minBeamEnergy;
	double locMax = m_minBeamEnergy + m_beamEnergyBinSize * (double)(locNBins);
	
	if(fabs(locMax-m_maxBeamEnergy) > 1.e-6) {
		printf("\n\nWarning: Beam energy bin edges of flux histogram are inconsistent with specificied bin widths.\n");
		printf("  Specificed energy range: %.3f-%.3f\n", m_minBeamEnergy, m_maxBeamEnergy);
		printf("  Flux weights histogram min-max: %.3f-%.3f\n\n", locMin, locMax);
	}
	
	h_fluxWeights = new TH1F("fluxWeights", "Fractional Photon Flux; E_{#gamma} [GeV]", locNBins, locMin, locMax);
	h_fluxWeights->GetYaxis()->SetTitle(Form("w_{i} = #Gamma(i) / #int_{%.2f GeV}^{%.2f GeV}#Gamma(E)dE", 
		m_minBeamEnergy, m_maxBeamEnergy));
	h_fluxWeights->GetYaxis()->SetTitleOffset(1.5);
	h_fluxWeights->GetYaxis()->CenterTitle(true);
	h_fluxWeights->SetDirectory(0);
	
	return;
}
