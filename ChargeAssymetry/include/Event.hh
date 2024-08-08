#ifndef Event_H_
#define Event_H_


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

class Event {

		public:
			Event();
			Event(TChain *trev);
			~Event();
			void TopAntiTopPair(TLorentzVector &TopAntiTop,Double_t &RapiDiff);
			std::vector<ParticleType>ParticleStack;
			float SetChannelAndCuts(std::string SelChParticle_temp="Ele",int nJetCut_temp=4, int nBJetCut_temp=1, int nPhoCut_temp=1);
			std::string GetChannel(){return SelChParticle;}
			
			bool FindSelectedEvents(std::vector<double> &BinArray,std::vector<double> &sel_event_eta_plus,std::vector<double> &sel_event_eta_minus,std::string BinType="TopMass");
			void SetDebugMode(){debug=1;}
			void ChargeAssymetry(std::vector<double> &sel_event_eta_plus,std::vector<double> &sel_event_eta_minus,std::vector<double>&CAssym,std::vector<double> &BinError);

			
			TLorentzVector GetPhoton()
			{
			  Photon.SetPtEtaPhiM(this->phoEt->at(0), this->phoEta->at(0), this->phoPhi->at(0), 0.0);
			  return Photon;
			}
			
			TLorentzVector GetTop()
			{
        		 
        		  if (this->TopLep_charge > 0)
      			  {
        			Top.SetPtEtaPhiM(this->TopLep_pt, this->TopLep_eta, this->TopLep_phi, this->Mt_blgammaMET);

     		 	  }
      			  else if (this->TopLep_charge < 0)
      			  {

        			Top.SetPtEtaPhiM(this->TopHad_pt, this->TopHad_eta, this->TopHad_phi, this->M_bjj);
      			  }
        		 

			  return Top;
			}
			TLorentzVector GetAntiTop()
			{
			  if (this->TopLep_charge > 0)
      			  {
        			AntiTop.SetPtEtaPhiM(this->TopHad_pt, this->TopHad_eta, this->TopHad_phi, this->M_bjj);
     		 	  }
      			  else if (this->TopLep_charge < 0)
      			  {
        			AntiTop.SetPtEtaPhiM(this->TopLep_pt, this->TopLep_eta, this->TopLep_phi, this->Mt_blgammaMET);
      			  }
			  return AntiTop;
			}
			
			TLorentzVector GetMuon()
			{
			  Muon.SetPtEtaPhiM(this->muPt->at(0), this->muEta->at(0), this->muPhi->at(0), 0.1057);
			  return Muon;
			}
			
			TLorentzVector GetElectron()
			{
			  Electron.SetPtEtaPhiM(this->elePt->at(0), this->eleEta->at(0), this->elePhi->at(0), 0.000511);
			  return Electron;
			}
			
			TLorentzVector GetTopAntiTop()
			{
			
			return GetTop()+GetAntiTop();
			}
			
			// Data File Variables
			TChain *tr;
			std::vector<Float_t> *elePt = 0;
			std::vector<Float_t> *eleEta = 0;
			std::vector<Float_t> *elePhi = 0;
			
			std::vector<Float_t> *muPt = 0;
			std::vector<Float_t> *muEta = 0;
			std::vector<Float_t> *muPhi = 0;
			
			std::vector<Float_t> *jetPt = 0;
			
			std::vector<Float_t> *phoEt = 0;
			std::vector<Float_t> *phoEta = 0;
			std::vector<Float_t> *phoPhi = 0;

			Int_t nEle = 0;
			Int_t nMu = 0;
			Int_t nJet = 0;
			Int_t nBJet = 0;
			Int_t nPho = 0;
			Int_t lumis = 0;
			Float_t evtWeight = 0;
			Float_t btagWeight_1a = 0;
			Float_t prefireSF = 0;
			Float_t muEffWeight = 0;
			Float_t eleEffWeight = 0;
			Float_t TopHad_pt = 0;
			Float_t TopLep_pt = 0;
			Float_t TopHad_eta = 0;
			Float_t TopLep_eta = 0;
			Float_t TopLep_phi = 0;
			Float_t TopHad_phi = 0;
			Float_t Mt_blgammaMET = 0; // toplepmass
			Float_t M_bjj = 0;         // tophadmass
			Float_t TopLep_charge = 0;
			Float_t PUweight = 0;
			
			Float_t Rapidity_T_ele; // taking T for top quark
			Float_t Rapidity_t_ele; // taking t for antitop quark
			Float_t Rapidity_T_mu;
			Float_t Rapidity_t_mu;
			Float_t Rapidity_T1; // selection 1
			Float_t Rapidity_t1; // selection 1
			Float_t Rapidity_T2; // selection 2
			Float_t Rapidity_t2; // selection 2
			Float_t PhoRapidity;
			Float_t Mass_top;
			Float_t Mass_antitop;
			Float_t ttbar_mass;

			bool passPresel_Ele;
			bool passPresel_Mu;
			
			// Class variables
			
			
			
				
		private: 
			std::string SelChParticle="Ele";
			int nJetCut=4;
			int nBJetCut=1;
			int nPhoCut=1;  
			bool debug=0;
						
			TLorentzVector Top;
			TLorentzVector AntiTop; 
			TLorentzVector TopAntiTop; 
			TLorentzVector Photon;
			TLorentzVector Muon;  
			TLorentzVector Electron; 



};
}
#endif
