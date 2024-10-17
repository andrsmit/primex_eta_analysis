#ifndef _ETAGG_ANA_
#define _ETAGG_ANA_

#define MAX_BEAM 400
#define MAX_FCAL 100
#define MAX_BCAL 100
#define MAX_TOF  100
#define MAX_SC   100
#define MAX_MC   5

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom3.h"

#include "particleType.h"

class EtaAna {
	private:
		
		int m_runNumber;
		TRandom3 *m_random;
		
		// Geometry:
		
		int    m_phaseVal;
		Particle_t m_Target;
		double m_targetLength, m_targetDensity, m_targetAtten;
		
		TVector3 m_fcalFace;
		TVector3 m_fcalCorrection;
		TVector3 m_vertex;
		
		// Constants:
		
		static constexpr double m_c = TMath::C() * 1.e-7; // [cm/ns]
		static constexpr double m_e = 0.510998928e-3;     // [GeV]
		static constexpr double m_fcalBlockSize = 4.0157; // [cm]
		
		// Variables:
		
		string m_outputFileName;
		
		TFile *m_inputFile;
		TTree *m_tree;
		
		int m_event;
		
		int  LoadTree();
		int  SetGeometry();
		void ReadEvent();
		int  CheckEventMultiplicities();
		void EtaggAnalysis();
		int  AcceptRejectEvent();
		
		TVector3 GetFCALPosition(int index);
		TVector3 GetBCALPosition(int index);
		
		int FCALFiducialCut(TVector3 pos, double cutLayer);
		
		void CheckTOFMatch(TVector3 pos, double &dxMin, double &dyMin, double &dtMin, double rfTimingCut);
		
		double GetEnergyAfterRecoil(double eb, double theta, double m0, double mp);
		double GetFCALEnergyResolution(double e);
		
		// TTree Variables:
		int    m_eventNum;
		double m_rfTime;
		int    m_nbeam;
		int    m_tagSystem[MAX_BEAM];
		int    m_tagCounter[MAX_BEAM];
		double m_beamE[MAX_BEAM];
		double m_beamT[MAX_BEAM];
		double m_accScaleFactor[MAX_BEAM];
		
		int    m_nfcal;
		double m_fcalE[MAX_FCAL];
		double m_fcalX[MAX_FCAL];
		double m_fcalY[MAX_FCAL];
		double m_fcalZ[MAX_FCAL];
		double m_fcalT[MAX_FCAL];
		int    m_fcalNblocks[MAX_FCAL];
		double m_fcalE1E9[MAX_FCAL];
		double m_fcalE9E25[MAX_FCAL];
		
		int    m_nbcal;
		double m_bcalE[MAX_BCAL];
		double m_bcalX[MAX_BCAL];
		double m_bcalY[MAX_BCAL];
		double m_bcalZ[MAX_BCAL];
		double m_bcalT[MAX_BCAL];
		
		int    m_ntof;
		double m_tofX[MAX_TOF];
		double m_tofY[MAX_TOF];
		double m_tofZ[MAX_TOF];
		double m_tofT[MAX_TOF];
		
		int    m_nsc;
		double m_scSector[MAX_SC];
		double m_scPhi[MAX_SC];
		double m_scdE[MAX_SC];
		double m_scT[MAX_SC];
		double m_scPulseHeight[MAX_SC];
		
		int    m_nmc;
		int    m_mcPDGType[MAX_MC];
		double m_mcX[MAX_MC];
		double m_mcY[MAX_MC];
		double m_mcZ[MAX_MC];
		double m_mcT[MAX_MC];
		double m_mcE[MAX_MC];
		double m_mcP[MAX_MC];
		double m_mcTheta[MAX_MC];
		double m_mcPhi[MAX_MC];
		
		// FCAL channel number:
		
		int m_fcalChannelNumber[59][59];
		int m_fcalRow[3481];
		int m_fcalColumn[3481];
		
		TVector2 m_fcalPositionOnFace[59][59];
		
		// Histograms:
		
		TH1F *h_mcVertex, *h_mcVertexAccepted;
		TH1F *h_mcReactionWeight;
		
		TH1F *h_fcalRFdt, *h_bcalRFdt, *h_tofRFdt, *h_scRFdt;
		TH1F *h_beamRFdt, *h_beamRFdt_cut;
		
		static const int m_nVetos = 8;
		TH2F *h_elasticity[m_nVetos];
		TH2F *h_mgg[m_nVetos];
		TH2F *h_mggConstr[m_nVetos];
		TH2F *h_mggConstr_coh[m_nVetos];
		TH2F *h_mggConstr_etaCut[m_nVetos];
		
		TH2F *h_xy1, *h_xy2;
		
		//----------------------------------------------------------------------------------------------//
		
	public:
		
		// Cuts:
		
		double m_FCALRFCut;
		double m_BCALRFCut;
		double m_BeamRFCut;
		double m_TOFRFCut;
		double m_FCALEnergyCut;
		double m_BCALEnergyCut;
		double m_minBeamEnergyCut, m_maxBeamEnergyCut;
		double m_FCALTOFCut;
		double m_SCDeltaPhiCut;
		double m_BCALDeltaPhiCut;
		double m_ElasMean_p0;
		double m_ElasMean_p1;
		double m_ElasWidth;
		double m_ElasSigmaCut;
		
		double beamBunchesMain = 1.0;
		double beamBunchesAcc  = 5.0; // how many bunches to use on each side for accidental subtraction
		
		EtaAna();
		~EtaAna(){};
		
		int GetPrimexPhase(int runNumber);
		
		// default analysis:
		void InitHistograms();
		void RunAnalysis(TString inputFileName);
		void ResetHistograms();
		void WriteHistograms();
		
		int SetRunNumber(int runNumber);
		void SetOutputFileName(string fileName);
};

#endif
