#ifndef Analysis_H_
#define Analysis_H_


#include<iostream>
#include<vector>
#include <iterator>
#include "TFile.h"
#include"TF1.h"
#include"TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TH3F.h"
#include "TH3.h"
#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TVectorD.h"
#include "TBox.h"
#include "TArrow.h"
#include <functional>
#include <TApplication.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <THttpServer.h>
#include <TMemFile.h>
#include <TNtuple.h>
#include "TROOT.h"
#include<string>
#include<fstream>
#include<cmath>
#include <TLorentzVector.h>
#include<TChain.h>


#include "TopVar.hh"
namespace TopAnalysis{

class Analysis {
         
         public:  
          	   Analysis();
          	  ~Analysis();
          	   TH1F *ChargeAssy_vs_TopMassBin=new TH1F("massbin","Charge Assymetry vs Top-AntiTop Mass (t#bar{t}+#gamma)",7,300,1000);
          	   //TH1F *ChargeAssy_vs_TopPtBin=
          	   //TH1F *ChargeAssy_vs_PhoPtBin=
          	  
          	  
         private: 	
};
}
#endif
