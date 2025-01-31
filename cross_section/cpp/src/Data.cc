#include "EtaAnalyzer.h"
#include "MggFitter.h"
#include "CrossSection.h"

TString    dataDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/analyze_trees/rootFiles";
TString   etamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles";
TString omegamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/omega_gammapi0/analyze_trees/rootFiles";

int EtaAnalyzer::GetDataHistograms(int anaOption, TString histName)
{
	TString fieldString = "bfield";
	if(m_phase==1) fieldString = "nobfield";
	
	TString anaString = "";
	switch(anaOption) {
		case 1:
			anaString = "_matrix";
			break;
		case 2:
			anaString = "_FCAL";
			break;
		case 3:
			anaString = "_BCAL";
			break;
		case 4:
			anaString = "_BEAM";
			break;
		case 5:
			anaString = "_TOF";
			break;
		default:
			anaString = "";
			break;
	}
	
	TString fullTargetFileName  = Form("%s/phase%d/full_target_%s%s.root",  dataDirectory.Data(), m_phase, 
		fieldString.Data(), anaString.Data());
	
	TString emptyTargetFileName = Form("%s/phase%d/empty_target_%s%s.root", dataDirectory.Data(), m_phase, 
		fieldString.Data(), anaString.Data());
	
	if(gSystem->AccessPathName(fullTargetFileName.Data()) || gSystem->AccessPathName(emptyTargetFileName.Data()))
		return 1;
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (anaOption=1) is fundamentally different from the rest:
	
	if(anaOption==1) {
		
		//-------------------------------------//
		// Get the full target data:
		
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)fullTargetFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of data histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", minBeamEnergy-locMinimumEnergy, maxBeamEnergy-locMaximumEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_mggVsThetaFull = (TH2F*)h3Full->Project3D("yx")->Clone("mgg_vs_theta_full");
		h_mggVsThetaFull->SetDirectory(0);
		
		fullTargetFile->Close();
		
		//-------------------------------------//
		// Do the same for the empty target data:
		
		TFile *emptyTargetFile = new TFile(emptyTargetFileName.Data(), "READ");
		
		TH3F *h3Empty = (TH3F*)emptyTargetFile->Get("invmassMatrix")->Clone("h3Empty");
		
		locBinSize = h3Empty->GetZaxis()->GetBinWidth(1);
		ebin1 = h3Empty->GetZaxis()->FindBin(minBeamEnergy + 0.5*locBinSize);
		ebin2 = h3Empty->GetZaxis()->FindBin(maxBeamEnergy - 0.5*locBinSize);
		
		h3Empty->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_mggVsThetaEmpty = (TH2F*)h3Empty->Project3D("yx")->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		
		emptyTargetFile->Close();
	}
	else {
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		h_mggVsThetaFull = (TH2F*)fullTargetFile->Get(Form("%s",histName.Data()))->Clone("mgg_vs_theta_full");
		h_mggVsThetaFull->SetDirectory(0);
		fullTargetFile->Close();
		
		TFile *emptyTargetFile = new TFile(emptyTargetFileName.Data(), "READ");
		
		h_mggVsThetaEmpty = (TH2F*)emptyTargetFile->Get(Form("%s",histName.Data()))->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		emptyTargetFile->Close();
	}
	
	// scale empty target data by photon flux ratio:
	
	h_mggVsThetaEmpty->Scale(m_emptyTargetFluxRatio);
	
	return 0;
}

int EtaAnalyzer::LoadLineshapes(int anaOption, TString histName)
{
	if(m_fitOption_signal>=5) 
	{
		if(LoadEtaLineshape(anaOption, histName)) return 1;
	}
	
	if(m_fitOption_omega>1) 
	{
		if(LoadOmegaLineshape(anaOption, histName)) return 1;
	}
	
	return 0;
}
int EtaAnalyzer::LoadEtaLineshape(int anaOption, TString histName)
{
	TString anaString = "";
	switch(anaOption) {
		case 1:
			anaString = "_matrix";
			break;
		case 2:
			anaString = "_FCAL";
			break;
		case 3:
			anaString = "_BCAL";
			break;
		case 4:
			anaString = "_BEAM";
			break;
		case 5:
			anaString = "_TOF";
			break;
		default:
			anaString = "";
			break;
	}
	
	TString mcFileName  = Form("%s/phase%d/phase%d%s.root",  etamcDirectory.Data(), 1, 1, anaString.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (anaOption=1) is fundamentally different from the rest:
	
	if(anaOption==1) {
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)mcFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of eta mc histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", minBeamEnergy-locMinimumEnergy, maxBeamEnergy-locMaximumEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_etaLineshape = (TH2F*)h3Full->Project3D("yx")->Clone("etaLineshape");
		h_etaLineshape->SetDirectory(0);
		
		mcFile->Close();
	}
	else {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		h_etaLineshape = (TH2F*)mcFile->Get(Form("%s",histName.Data()))->Clone("etaLineshape");
		h_etaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	
	return 0;
}

int EtaAnalyzer::LoadOmegaLineshape(int anaOption, TString histName)
{
	TString anaString = "";
	switch(anaOption) {
		case 1:
			anaString = "_matrix";
			break;
		case 2:
			anaString = "_FCAL";
			break;
		case 3:
			anaString = "_BCAL";
			break;
		case 4:
			anaString = "_BEAM";
			break;
		case 5:
			anaString = "_TOF";
			break;
		default:
			anaString = "";
			break;
	}
	
	TString mcFileName  = Form("%s/phase%d%s.root",  omegamcDirectory.Data(), 1, anaString.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (anaOption=1) is fundamentally different from the rest:
	/*
	TFile *f_omega_lineshape = new TFile("../omega_lineshape_new.root", "READ");
	h_omegaLineshape = (TH2F*)f_omega_lineshape->Get("omega_lineshape_2d");
	h_omegaLineshape->SetDirectory(0);
	f_omega_lineshape->Close();
	*/
	if(anaOption==1) {
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)mcFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of omega mc histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", minBeamEnergy-locMinimumEnergy, maxBeamEnergy-locMaximumEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_omegaLineshape = (TH2F*)h3Full->Project3D("yx")->Clone("omegaLineshape");
		h_omegaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	else {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		h_omegaLineshape = (TH2F*)mcFile->Get(Form("%s",histName.Data()))->Clone("omegaLineshape");
		h_omegaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	
	return 0;
}
