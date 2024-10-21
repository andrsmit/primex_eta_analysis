/*
ROOT Macro to calculate the integrated photon flux within a specified energy range for the PrimEx run periods.
Example usage:
	root -l -b -q 'GetIntegratedFlux_PrimEx(1, 8.0, 10.9)'

This will sum the total flux measured in the energy range 8-10.9 GeV for the first phase of PrimEx.
Note: To run, the relevant run lists (stored on the /work disk must be accessible).
*/

#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <stdlib.h> 
#include <fstream>

#include "TString.h"

double tagh_en[274], tagm_en[102];
double endpointEnergy, endpointEnergyCalib;

TString fluxDir = "/work/halld/home/andrsmit/primex_eta_analysis/photon_flux";
TString runListDirectory = "/work/halld/home/andrsmit/run_lists";
TString targetStr[2]     = {"full","empty"};

double targetThickness = 5.408e+23;
	// number of atoms per cm^2, assuming density = 0.1217 g/cm3 and length = 29.535cm

double cm2pb = 1.e-36;

int GetCounterEnergies(int phase, int run);
int GetIntegratedFlux(int phase, int run, double minE, double maxE, double &flux, double &fluxE);

void GetIntegratedFlux_PrimEx(int primexPhase=1, double minBeamEnergy=8.0, double maxBeamEnergy=12.0) 
{
	printf("\nSumming photon flux between %f-%f GeV for PrimEx-eta phase %d\n", minBeamEnergy, maxBeamEnergy, primexPhase);
	std::cout << "" << std::endl;
	
	//---------------------------------------------------------------//
	// Get run list:
	TString locRunList[2] = {"",""};
	switch(primexPhase) {
		case 1:
			for(int i=0; i<2; i++) {
				locRunList[i] = Form("%s/primex_phase1/he_%s_nobfield.txt", 
					runListDirectory.Data(), targetStr[i].Data());
			}
			break;
		case 2:
			for(int i=0; i<2; i++) {
				locRunList[i] = Form("%s/primex_phase2/he_%s_bfield.txt", 
					runListDirectory.Data(), targetStr[i].Data());
				//locRunList[i] = Form("%s/primex_phase2/he_%s_nobfield.txt", 
				//	runListDirectory.Data(), targetStr[i].Data());
			}
			break;
		case 3:
			for(int i=0; i<2; i++) {
				locRunList[i] = Form("%s/primex_phase3/he_%s_bfield.txt", 
					runListDirectory.Data(), targetStr[i].Data());
			}
			break;
		default:
			std::cout << "Invalid PrimEx run period specified. Please choose between 1, 2, or 3." << std::endl;
			exit(1);
	}
	
	//----------------------------------------------------------------//
	
	pair<double,double>  fullTargetFlux = {0.0, 0.0};
	pair<double,double> emptyTargetFlux = {0.0, 0.0};
	
	for(int itarget=0; itarget<2; itarget++) {
		printf("Processing %s target runs...\n",targetStr[itarget].Data());
		ifstream locInf(locRunList[itarget].Data());
		if(!locInf.is_open()) {
			std::cout << "Problem accessing run list." << std::endl;
			std::cout << "  filename: " << locRunList[itarget] << std::endl;
			exit(1);
		}
		
		double locIntegratedFlux  = 0.0;
		double locIntegratedFluxE = 0.0;
		
		int locRun;
		double locFlux, locFluxE;
		while(locInf >> locRun) {
			if(GetCounterEnergies(primexPhase, locRun)) {
				std::cout << "Problem loading tagger counters energy scale for run " << locRun << std::endl;
				continue;
			}
			//std::cout << "  " << locRun << std::endl;
			//std::cout << "    endpoint_energy = " << endpointEnergy << " GeV" << std::endl;
			if(GetIntegratedFlux(primexPhase, locRun, minBeamEnergy, maxBeamEnergy, locFlux, locFluxE)) {
				std::cout << "Problem getting flux for run " << locRun << std::endl;
				continue;
			}
			locIntegratedFlux  += locFlux;
			locIntegratedFluxE += locFluxE;
		}
		locInf.close();
		
		if(targetStr[itarget]=="full") {
			fullTargetFlux.first  = locIntegratedFlux;
			fullTargetFlux.second = sqrt(locIntegratedFluxE);
		} else if(targetStr[itarget]=="empty") {
			emptyTargetFlux.first  = locIntegratedFlux;
			emptyTargetFlux.second = sqrt(locIntegratedFluxE);
		}
	}
	
	double integratedLum = fullTargetFlux.first * targetThickness * cm2pb;
	
	printf("\n\n\n");
	printf("Integrated flux for full target runs: %e +/- %e\n", fullTargetFlux.first, fullTargetFlux.second);
	printf("Integrated luminosity: %f pb-1\n", integratedLum);
	printf("  Empty-target flux ratio (full/empty): %f\n", fullTargetFlux.first/emptyTargetFlux.first);
	printf("\n\n\n");
	
	return 0;
}

int GetCounterEnergies(int phase, int run) {
	
	//----------------------------------------//
	// Get endpoint_calib from run number:
	
	endpointEnergyCalib = 0.0;
	if(phase==1) {
		if(run<61434) endpointEnergyCalib = 11.6061;
		else endpointEnergyCalib = 11.1689;
	} else if(phase==2) {
		endpointEnergyCalib = 10.042;
	} else if(phase==3) {
		endpointEnergyCalib = 11.555;
	} else {
		return 1;
	}
	
	//----------------------------------------//
	// Get endpoint_energy from txtFile:
	
	endpointEnergy = 0.0;
	TString endpointFilename = fluxDir+Form("/endpointEnergies/phase%d/%d.txt", phase, run);
	if(gSystem->AccessPathName(endpointFilename.Data())) {
		return 1;
	}
	ifstream infEndpoint(endpointFilename.Data());
	
	string commentLine;
	getline(infEndpoint, commentLine);
	
	infEndpoint >> endpointEnergy;
	infEndpoint.close();
	
	//---------------------------------------//
	// Get energy scales from txtFile:
	
	TString tagh_xscale_fname = fluxDir+Form("/xscales/phase%d/primex_tagh.txt", phase);
	TString tagm_xscale_fname = fluxDir+Form("/xscales/phase%d/primex_tagm.txt", phase);
	
	ifstream infTAGH(tagh_xscale_fname.Data());
	ifstream infTAGM(tagm_xscale_fname.Data());
	
	// first check that input filenames exist:
	if(!infTAGH.good() || !infTAGM.good()) {
		return 1;
	}
	
	int a; double b, c;
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		infTAGH >> a >> b >> c;
		double deltaE = endpointEnergy - endpointEnergyCalib;
		double emin   = b * endpointEnergyCalib  +  deltaE;
		double emax   = c * endpointEnergyCalib  +  deltaE;
		tagh_en[tagh_counter-1] = 0.5 * (emin + emax);
	}
	infTAGH.close();
	
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		infTAGM >> a >> b >> c;
		double deltaE = endpointEnergy - endpointEnergyCalib;
		double emin   = b * endpointEnergyCalib  +  deltaE;
		double emax   = c * endpointEnergyCalib  +  deltaE;
		tagm_en[tagm_counter-1] = 0.5 * (emin + emax);
	}
	infTAGM.close();
	
	return 0;
}

int GetIntegratedFlux(int phase, int run, double minE, double maxE, double &flux, double &fluxE) {
	
	/*
	Function to get the integrated flux between 'minE' and 'maxE' for specified 'run'.
	This function will attempt to read the total flux for each counter stored in text files for each
	run in this directory: "/work/halld/home/andrsmit/primex_eta_analysis/photon_flux/txtFiles"
	For each run, the endpoint energy is needed to calculate the energy associated with each tagger counter.
	If the midpoint energy of the counter is between minE and maxE, the flux of photons measured
	by that counter is added into the total 'flux' for each run.
	*/
	
	TString tagh_flux_filename = fluxDir+Form("/txtFiles/phase%d/%d_tagh_ps_acc_cor.txt", phase, run);
	TString tagm_flux_filename = fluxDir+Form("/txtFiles/phase%d/%d_tagm_ps_acc_cor.txt", phase, run);
	
	ifstream infTAGH(tagh_flux_filename.Data());
	ifstream infTAGM(tagm_flux_filename.Data());
	
	// first check that input filenames exist:
	if(!infTAGH.good() || !infTAGM.good()) {
		return 1;
	}
	
	flux = 0.0, fluxE = 0.0;
	
	int a; double b, c;
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		infTAGH >> a >> b >> c;
		double locE = tagh_en[tagh_counter-1];
		if((minE<=locE) && (locE<maxE)) {
			flux  += b;
			fluxE += pow(c,2.0);
		}
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		infTAGM >> a >> b >> c;
		double locE = tagm_en[tagm_counter-1];
		if((minE<=locE) && (locE<maxE)) {
			flux  += b;
			fluxE += pow(c,2.0);
		}
	}
	
	return 0;
}
