/**************************************************************************                                                                                                                           
Code copied from $HALLD_SIM_HOME/src/libraries/UTILITIES
  (with slight alteration)

**************************************************************************/

#ifndef _MYREADCONFIG_H_
#define _MYREADCONFIG_H_

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#include <string>
#include <math.h>
#include "TFile.h"
#include "TObject.h"
#include "TSystem.h"
#include "TLorentzVector.h"

const Int_t MAX_LINE = 1000;

class MyReadConfig {
	private:
		Int_t   nLine;
		TString strLine[MAX_LINE];
		TString strCaLibPath;
		
	protected:
	public:
		MyReadConfig();
		MyReadConfig(Char_t *);
		virtual ~MyReadConfig(){};
		
		TString ExtractName(TString);
		
		TString GetConfigName(TString);
		int ReadConfigFile(const Char_t*);
		
		Double_t* GetConfig1Par(TString);
		Double_t* GetConfig2Par(TString);
		Double_t* GetConfig3Par(TString);
		Double_t* GetConfig4Par(TString);
		Double_t* GetConfig5Par(TString);
		Double_t* GetConfig6Par(TString);
};

#endif
