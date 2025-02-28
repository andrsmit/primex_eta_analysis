#include "EtaAnalyzer.h"
#include "MggFitter.h"
#include "CrossSection.h"

TString    dataDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/analyze_trees/rootFiles";
TString   etamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles";
TString omegamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/omega_gammapi0/analyze_trees/rootFiles";
TString   bggenDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana/analyze_trees/rootFiles";

int EtaAnalyzer::LoadDataHistograms()
{
	printf("\nREADING DATA HISTGRAMS...\n");
	
	TString fieldString = "bfield";
	if(m_phase==1) fieldString = "nobfield";
	
	TString anaString = "";
	switch(m_analysisOption) {
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
	
	printf("   Full Target: %s\n",  fullTargetFileName.Data());
	printf("  Empty Target: %s\n", emptyTargetFileName.Data());
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (m_analysisOption=1) is fundamentally different from the rest:
	
	if(m_analysisOption==1) {
		
		//-------------------------------------//
		// Get the full target data:
		
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)fullTargetFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of data histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
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
		ebin1 = h3Empty->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		ebin2 = h3Empty->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		h3Empty->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_mggVsThetaEmpty = (TH2F*)h3Empty->Project3D("yx")->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		
		emptyTargetFile->Close();
	}
	else {
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		h_mggVsThetaFull = (TH2F*)fullTargetFile->Get(Form("%s",m_mggHistName.Data()))->Clone("mgg_vs_theta_full");
		h_mggVsThetaFull->SetDirectory(0);
		fullTargetFile->Close();
		
		TFile *emptyTargetFile = new TFile(emptyTargetFileName.Data(), "READ");
		
		h_mggVsThetaEmpty = (TH2F*)emptyTargetFile->Get(Form("%s",m_mggHistName.Data()))->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		emptyTargetFile->Close();
	}
	
	// scale empty target data by photon flux ratio:
	
	h_mggVsThetaEmpty->Scale(m_emptyTargetFluxRatio);
	
	return 0;
}

int EtaAnalyzer::LoadLineshapes()
{
	printf("\nREADING LINESHAPES...\n");
	
	if(m_fitOption_signal>=5) 
	{
		if(LoadEtaLineshape()) return 1;
		if(m_fitOption_signal>5) {
			if(LoadEtaPionLineshape()) return 1;
			if(LoadEtaPionFraction()) return 1;
		}
	}
	
	if(m_fitOption_omega>=2) 
	{
		if(LoadOmegaLineshape()) return 1;
	}
	
	// Read the lineshape of omega->pi0+gamma from the first FDC package:
	if(m_fitOption_empty==1 && m_emptyFitOption_fdc>=2)
	{
		if(LoadFDCOmegaLineshape()) return 1;
	}
	return 0;
}
int EtaAnalyzer::LoadEtaLineshape()
{
	if(m_fitOption_signal==7) {
		TString mcFileName  = Form("%s/phase%d/Helium.root",  bggenDirectory.Data(), m_phase);
		
		if(gSystem->AccessPathName(mcFileName.Data())) return 1;
		
		printf("  Eta lineshape from %s\n", mcFileName.Data());
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		h_etaLineshape = (TH2F*)mcFile->Get("mgg_const_bggen_signal")->Clone("etaLineshape");
		h_etaLineshape->SetDirectory(0);
		mcFile->Close();
		return 0;
	}
	
	TString anaString = "";
	switch(m_analysisOption) {
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
	printf("  Eta lineshape from %s\n", mcFileName.Data());
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (m_analysisOption=1) is fundamentally different from the rest:
	
	if(m_analysisOption==1) {
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)mcFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of eta mc histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_etaLineshape = (TH2F*)h3Full->Project3D("yx")->Clone("etaLineshape");
		h_etaLineshape->SetDirectory(0);
		
		mcFile->Close();
	}
	else {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		h_etaLineshape = (TH2F*)mcFile->Get(Form("%s",m_mggHistName.Data()))->Clone("etaLineshape");
		h_etaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	
	return 0;
}

int EtaAnalyzer::LoadEtaPionLineshape()
{
	TString mcFileName  = Form("%s/phase%d/Helium.root",  bggenDirectory.Data(), m_phase);
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Eta+Pion lineshape from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_etaPionLineshape = (TH2F*)mcFile->Get("mgg_const_bggen_etapion")->Clone("etapionLineshape");
	h_etaPionLineshape->SetDirectory(0);
	
	// add in contribution from eta+2pion background channels:
	
	TH2F *h_eta2PionLS = (TH2F*)mcFile->Get("mgg_const_bggen_eta2pion")->Clone("eta2pionLineshape");
	h_etaPionLineshape->Add(h_eta2PionLS);
	
	mcFile->Close();
	
	return 0;
}

int EtaAnalyzer::LoadEtaPionFraction()
{
	TString mcFileName  = Form("%s/phase%d/Helium.root",  bggenDirectory.Data(), m_phase);
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Eta+Pion fraction from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	TH2F *h2_exclusive = (TH2F*)mcFile->Get("mgg_const_bggen_signal");
	TH1F *h1_exclusive = (TH1F*)h2_exclusive->ProjectionX("h1_exclusive");
	h1_exclusive->Rebin(m_rebinsTheta);
	
	TH2F *h2_etapion   = (TH2F*)mcFile->Get("mgg_const_bggen_etapion");
	TH1F *h1_etapion   = (TH1F*)h2_etapion->ProjectionX("h1_etapion");
	TH2F *h2_eta2pion  = (TH2F*)mcFile->Get("mgg_const_bggen_eta2pion");
	TH1F *h1_eta2pion  = (TH1F*)h2_eta2pion->ProjectionX("h1_eta2pion");
	h1_etapion->Add(h1_eta2pion);
	h1_etapion->Rebin(m_rebinsTheta);
	
	h1_etapion->Divide(h1_exclusive);
	
	h_EtaPionFraction_bggen = (TH1F*)h1_etapion->Clone("etaPionFraction");
	h_EtaPionFraction_bggen->SetDirectory(0);
	mcFile->Close();
	
	return 0;
}

int EtaAnalyzer::LoadOmegaLineshape()
{
	// Use the lineshape from bggen for the omega+other backgrounds:
	
	TString mcFileName  = Form("%s/phase%d/Helium.root",  bggenDirectory.Data(), 1);
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Omega lineshape from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_omegaLineshape = (TH2F*)mcFile->Get("mgg_const_bggen_omega")->Clone("omegaLineshape");
	h_omegaLineshape->SetDirectory(0);
	mcFile->Close();
	return 0;
	
	/*
	TString anaString = "";
	switch(m_analysisOption) {
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
	
	//
	// TEMPORARY:
	//
	anaString = "_matrix";
	
	TString mcFileName  = Form("%s/phase%d%s.root",  omegamcDirectory.Data(), 1, anaString.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (m_analysisOption=1) is fundamentally different from the rest:
	//if(m_analysisOption==1) {
	if(1) {
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)mcFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of omega mc histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_omegaLineshape = (TH2F*)h3Full->Project3D("yx")->Clone("omegaLineshape");
		h_omegaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	else {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		h_omegaLineshape = (TH2F*)mcFile->Get(Form("%s",m_mggHistName.Data()))->Clone("omegaLineshape");
		h_omegaLineshape->SetDirectory(0);
		mcFile->Close();
	}
	return 0;
	*/
}

int EtaAnalyzer::LoadFDCOmegaLineshape()
{
	// Use the lineshape from bggen for the omega+other backgrounds:
	
	TString mcFileName  = Form("%s/phase3_1st-FDC-pack.root",  omegamcDirectory.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_fdcOmegaLineshape = (TH2F*)mcFile->Get("mgg_const_veto_5")->Clone("fdcOmegaLineshape");
	h_fdcOmegaLineshape->SetDirectory(0);
	mcFile->Close();
	return 0;
}
