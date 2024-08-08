//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  5 13:38:28 2023 by ROOT version 6.12/07
// from TTree AnalysisTree/AnalysisTree
// found on file: /eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root
//////////////////////////////////////////////////////////

#ifndef genloop_h
#define genloop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class genloop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Float_t         PUweight;
   Float_t         PUweight_Up;
   Float_t         PUweight_Do;
   Float_t         prefireSF;
   Float_t         prefireSF_Up;
   Float_t         prefireSF_Do;
   vector<float>   *btagWeight;
   Float_t         btagWeight_1a;
   vector<float>   *btagWeight_b_Up;
   vector<float>   *btagWeight_b_Do;
   vector<float>   *btagWeight_l_Up;
   vector<float>   *btagWeight_l_Do;
   Float_t         btagWeight_1a_b_Up;
   Float_t         btagWeight_1a_b_Do;
   Float_t         btagWeight_1a_l_Up;
   Float_t         btagWeight_1a_l_Do;
   vector<float>   *btagSF;
   Float_t         muEffWeight;
   Float_t         muEffWeight_IdIso;
   Float_t         muEffWeight_Trig;
   Float_t         muEffWeight_Up;
   Float_t         muEffWeight_Do;
   Float_t         muEffWeight_IdIso_Up;
   Float_t         muEffWeight_IdIso_Do;
   Float_t         muEffWeight_Trig_Up;
   Float_t         muEffWeight_Trig_Do;
   Float_t         eleEffWeight;
   Float_t         eleEffWeight_IdReco;
   Float_t         eleEffWeight_Trig;
   Float_t         eleEffWeight_Up;
   Float_t         eleEffWeight_Do;
   Float_t         eleEffWeight_IdReco_Up;
   Float_t         eleEffWeight_IdReco_Do;
   Float_t         eleEffWeight_Trig_Up;
   Float_t         eleEffWeight_Trig_Do;
   vector<float>   *phoEffWeight;
   vector<float>   *phoEffWeight_Id;
   vector<float>   *phoEffWeight_eVeto;
   vector<float>   *loosePhoEffWeight;
   vector<float>   *loosePhoEffWeight_Id;
   vector<float>   *loosePhoEffWeight_eVeto;
   vector<float>   *phoNoIDEffWeight;
   vector<float>   *phoNoIDEffWeight_Id;
   vector<float>   *phoNoIDEffWeight_eVeto;
   vector<float>   *phoEffWeight_Up;
   vector<float>   *phoEffWeight_Do;
   vector<float>   *phoEffWeight_Id_Up;
   vector<float>   *phoEffWeight_Id_Do;
   vector<float>   *phoEffWeight_eVeto_Up;
   vector<float>   *phoEffWeight_eVeto_Do;
   vector<float>   *loosePhoEffWeight_Up;
   vector<float>   *loosePhoEffWeight_Do;
   vector<float>   *loosePhoEffWeight_Id_Up;
   vector<float>   *loosePhoEffWeight_Id_Do;
   vector<float>   *loosePhoEffWeight_eVeto_Up;
   vector<float>   *loosePhoEffWeight_eVeto_Do;
   vector<float>   *phoNoIDEffWeight_Up;
   vector<float>   *phoNoIDEffWeight_Do;
   vector<float>   *phoNoIDEffWeight_Id_Up;
   vector<float>   *phoNoIDEffWeight_Id_Do;
   vector<float>   *phoNoIDEffWeight_eVeto_Up;
   vector<float>   *phoNoIDEffWeight_eVeto_Do;
   Float_t         q2weight_Up;
   Float_t         q2weight_Do;
   Float_t         q2weight_nominal;
   vector<float>   *genScaleSystWeights;
   Float_t         pdfWeight;
   Float_t         pdfuncer;
   Float_t         pdfweight_Up;
   Float_t         pdfweight_Do;
   vector<float>   *pdfSystWeight;
   Float_t         ISRweight_Up;
   Float_t         ISRweight_Do;
   Float_t         FSRweight_Up;
   Float_t         FSRweight_Do;
   Float_t         evtWeight;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Float_t         genMET;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         nu_pz;
   Float_t         nu_pz_other;
   Float_t         WtransMass;
   Float_t         Mt_blgammaMET;
   Float_t         Mt_lgammaMET;
   Float_t         M_bjj;
   Float_t         M_bjjgamma;
   Float_t         M_jj;
   Bool_t          MassCuts;
   Float_t         TopHad_pt;
   Float_t         TopHad_eta;
   Float_t         TopHad_phi;
   Float_t         TopLep_pt;
   Float_t         TopLep_eta;
   Float_t         TopLep_phi;
   Float_t         TopLep_charge;
   Float_t         chi2;
   Float_t         DiphoMass;
   Float_t         DilepMass;
   Float_t         DilepDelR;
   Int_t           nPho;
   Int_t           nPhoBarrel;
   Int_t           nPhoEndcap;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoR9;
   vector<float>   *phoPhi;
   vector<bool>    *phoIsBarrel;
   vector<float>   *phoHoverE;
   vector<float>   *phoSIEIE;
   vector<float>   *phoPFChIso;
   vector<bool>    *phoTightID;
   vector<bool>    *phoMediumID;
   vector<int>     *phoGenMatchInd;
   vector<float>   *phoMassLepGamma;
   Int_t           nLoosePho;
   vector<float>   *loosePhoEt;
   vector<float>   *loosePhoEta;
   vector<float>   *loosePhoPhi;
   vector<bool>    *loosePhoIsBarrel;
   vector<float>   *loosePhoHoverE;
   vector<float>   *loosePhoSIEIE;
   vector<float>   *loosePhoPFChIso;
   vector<bool>    *loosePhoTightID;
   vector<bool>    *loosePhoMediumID;
   vector<bool>    *loosePhoLooseID;
   vector<float>   *loosePhoMVAId;
   vector<float>   *loosePhoMVAId17v1;
   vector<int>     *loosePhoGenMatchInd;
   vector<float>   *loosePhoMassLepGamma;
   vector<bool>    *loosePhoMediumIDFunction;
   vector<bool>    *loosePhoMediumIDPassHoverE;
   vector<bool>    *loosePhoMediumIDPassSIEIE;
   vector<bool>    *loosePhoMediumIDPassChIso;
   vector<bool>    *loosePhoMediumIDPassNeuIso;
   vector<bool>    *loosePhoMediumIDPassPhoIso;
   Int_t           nPhoNoID;
   vector<float>   *phoNoIDEt;
   vector<float>   *phoNoIDEta;
   vector<float>   *phoNoIDPhi;
   vector<bool>    *phoNoIDIsBarrel;
   vector<float>   *phoNoIDHoverE;
   vector<float>   *phoNoIDSIEIE;
   vector<float>   *phoNoIDPFChIso;
   vector<bool>    *phoNoIDTightID;
   vector<bool>    *phoNoIDMediumID;
   vector<bool>    *phoNoIDLooseID;
   vector<float>   *phoNoIDMVAId;
   vector<float>   *phoNoIDMVAId17v1;
   vector<int>     *phoNoIDGenMatchInd;
   vector<float>   *phoNoIDMassLepGamma;
   vector<bool>    *phoNoIDMediumIDFunction;
   vector<bool>    *phoNoIDMediumIDPassHoverE;
   vector<bool>    *phoNoIDMediumIDPassSIEIE;
   vector<bool>    *phoNoIDMediumIDPassChIso;
   vector<bool>    *phoNoIDMediumIDPassNeuIso;
   vector<bool>    *phoNoIDMediumIDPassPhoIso;
   Int_t           nEle;
   vector<float>   *elePt;
   vector<float>   *elePhi;
   vector<float>   *eleEta;
   vector<float>   *eleSCEta;
   vector<float>   *elePFRelIso;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<float>   *muPFRelIso;
   Int_t           nJet;
   Int_t           nfwdJet;
   Int_t           nBJet;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetDeepB;
   vector<int>     *jetGenJetIdx;
   vector<float>   *fwdJetPt;
   vector<float>   *fwdJetEta;
   vector<float>   *fwdJetPhi;
   vector<float>   *fwdJetMass;
   vector<float>   *dRPhotonJet;
   vector<float>   *dRPhotonLepton;
   vector<float>   *MPhotonLepton;
   vector<float>   *AnglePhotonLepton;
   Int_t           nGenPart;
   vector<float>   *genPt;
   vector<float>   *genEta;
   vector<float>   *genPhi;
   vector<float>   *genMass;
   vector<int>     *genStatus;
   vector<int>     *genStatusFlag;
   vector<int>     *genPDGID;
   vector<int>     *genMomIdx;
   Int_t           nGenJet;
   vector<float>   *genJetPt;
   vector<float>   *genJetEta;
   vector<float>   *genJetPhi;
   vector<float>   *genJetMass;
   vector<int>     *genJetPartonFlavour;
   Double_t        M3;
   Double_t        M3_gamma;
   Double_t        HT;
   Bool_t          passPresel_Ele;
   Bool_t          passPresel_Mu;
   Bool_t          passAll_Ele;
   Bool_t          passAll_Mu;
   Bool_t          inHEMVeto;
   vector<bool>    *photonIsGenuine;
   vector<bool>    *photonIsMisIDEle;
   vector<bool>    *photonIsHadronicPhoton;
   vector<bool>    *photonIsHadronicFake;
   vector<bool>    *loosePhotonIsGenuine;
   vector<bool>    *loosePhotonIsMisIDEle;
   vector<bool>    *loosePhotonIsHadronicPhoton;
   vector<bool>    *loosePhotonIsHadronicFake;
   vector<bool>    *photonNoIDIsGenuine;
   vector<bool>    *photonNoIDIsMisIDEle;
   vector<bool>    *photonNoIDIsHadronicPhoton;
   vector<bool>    *photonNoIDIsHadronicFake;
   vector<int>     *photonParentage;
   vector<int>     *photonParentPID;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_PUweight_Up;   //!
   TBranch        *b_PUweight_Do;   //!
   TBranch        *b_prefireSF;   //!
   TBranch        *b_prefireSF_Up;   //!
   TBranch        *b_prefireSF_Do;   //!
   TBranch        *b_btagWeight;   //!
   TBranch        *b_btagWeight_1a;   //!
   TBranch        *b_btagWeight_b_Up;   //!
   TBranch        *b_btagWeight_b_Do;   //!
   TBranch        *b_btagWeight_l_Up;   //!
   TBranch        *b_btagWeight_l_Do;   //!
   TBranch        *b_btagWeight_1a_b_Up;   //!
   TBranch        *b_btagWeight_1a_b_Do;   //!
   TBranch        *b_btagWeight_1a_l_Up;   //!
   TBranch        *b_btagWeight_1a_l_Do;   //!
   TBranch        *b_btagSF;   //!
   TBranch        *b_muEffWeight;   //!
   TBranch        *b_muEffWeight_IdIso;   //!
   TBranch        *b_muEffWeight_Trig;   //!
   TBranch        *b_muEffWeight_Up;   //!
   TBranch        *b_muEffWeight_Do;   //!
   TBranch        *b_muEffWeight_IdIso_Up;   //!
   TBranch        *b_muEffWeight_IdIso_Do;   //!
   TBranch        *b_muEffWeight_Trig_Up;   //!
   TBranch        *b_muEffWeight_Trig_Do;   //!
   TBranch        *b_eleEffWeight;   //!
   TBranch        *b_eleEffWeight_IdReco;   //!
   TBranch        *b_eleEffWeight_Trig;   //!
   TBranch        *b_eleEffWeight_Up;   //!
   TBranch        *b_eleEffWeight_Do;   //!
   TBranch        *b_eleEffWeight_IdReco_Up;   //!
   TBranch        *b_eleEffWeight_IdReco_Do;   //!
   TBranch        *b_eleEffWeight_Trig_Up;   //!
   TBranch        *b_eleEffWeight_Trig_Do;   //!
   TBranch        *b_phoEffWeight;   //!
   TBranch        *b_phoEffWeight_Id;   //!
   TBranch        *b_phoEffWeight_eVeto;   //!
   TBranch        *b_loosePhoEffWeight;   //!
   TBranch        *b_loosePhoEffWeight_Id;   //!
   TBranch        *b_loosePhoEffWeight_eVeto;   //!
   TBranch        *b_phoNoIDEffWeight;   //!
   TBranch        *b_phoNoIDEffWeight_Id;   //!
   TBranch        *b_phoNoIDEffWeight_eVeto;   //!
   TBranch        *b_phoEffWeight_Up;   //!
   TBranch        *b_phoEffWeight_Do;   //!
   TBranch        *b_phoEffWeight_Id_Up;   //!
   TBranch        *b_phoEffWeight_Id_Do;   //!
   TBranch        *b_phoEffWeight_eVeto_Up;   //!
   TBranch        *b_phoEffWeight_eVeto_Do;   //!
   TBranch        *b_loosePhoEffWeight_Up;   //!
   TBranch        *b_loosePhoEffWeight_Do;   //!
   TBranch        *b_loosePhoEffWeight_Id_Up;   //!
   TBranch        *b_loosePhoEffWeight_Id_Do;   //!
   TBranch        *b_loosePhoEffWeight_eVeto_Up;   //!
   TBranch        *b_loosePhoEffWeight_eVeto_Do;   //!
   TBranch        *b_phoNoIDEffWeight_Up;   //!
   TBranch        *b_phoNoIDEffWeight_Do;   //!
   TBranch        *b_phoNoIDEffWeight_Id_Up;   //!
   TBranch        *b_phoNoIDEffWeight_Id_Do;   //!
   TBranch        *b_phoNoIDEffWeight_eVeto_Up;   //!
   TBranch        *b_phoNoIDEffWeight_eVeto_Do;   //!
   TBranch        *b_q2weight_Up;   //!
   TBranch        *b_q2weight_Do;   //!
   TBranch        *b_q2weight_nominal;   //!
   TBranch        *b_genScaleSystWeights;   //!
   TBranch        *b_pdfWeight;   //!
   TBranch        *b_pdfuncer;   //!
   TBranch        *b_pdfweight_Up;   //!
   TBranch        *b_pdfweight_Do;   //!
   TBranch        *b_pdfSystWeight;   //!
   TBranch        *b_ISRweight_Up;   //!
   TBranch        *b_ISRweight_Do;   //!
   TBranch        *b_FSRweight_Up;   //!
   TBranch        *b_FSRweight_Do;   //!
   TBranch        *b_evtWeight;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_nu_pz;   //!
   TBranch        *b_nu_pz_other;   //!
   TBranch        *b_WtransMass;   //!
   TBranch        *b_Mt_blgammaMET;   //!
   TBranch        *b_Mt_lgammaMET;   //!
   TBranch        *b_M_bjj;   //!
   TBranch        *b_M_bjjgamma;   //!
   TBranch        *b_M_jj;   //!
   TBranch        *b_MassCuts;   //!
   TBranch        *b_TopHad_pt;   //!
   TBranch        *b_TopHad_eta;   //!
   TBranch        *b_TopHad_phi;   //!
   TBranch        *b_TopLep_pt;   //!
   TBranch        *b_TopLep_eta;   //!
   TBranch        *b_TopLep_phi;   //!
   TBranch        *b_TopLep_charge;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_DiphoMass;   //!
   TBranch        *b_DilepMass;   //!
   TBranch        *b_DilepDelR;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_nPhoBarrel;   //!
   TBranch        *b_nPhoEndcap;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoIsBarrel;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSIEIE;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoTightID;   //!
   TBranch        *b_phoMediumID;   //!
   TBranch        *b_phoGenMatchInd;   //!
   TBranch        *b_phoMassLepGamma;   //!
   TBranch        *b_nLoosePho;   //!
   TBranch        *b_loosePhoEt;   //!
   TBranch        *b_loosePhoEta;   //!
   TBranch        *b_loosePhoPhi;   //!
   TBranch        *b_loosePhoIsBarrel;   //!
   TBranch        *b_loosePhoHoverE;   //!
   TBranch        *b_loosePhoSIEIE;   //!
   TBranch        *b_loosePhoPFChIso;   //!
   TBranch        *b_loosePhoTightID;   //!
   TBranch        *b_loosePhoMediumID;   //!
   TBranch        *b_loosePhoLooseID;   //!
   TBranch        *b_loosePhoMVAId;   //!
   TBranch        *b_loosePhoMVAId17v1;   //!
   TBranch        *b_loosePhoGenMatchInd;   //!
   TBranch        *b_loosePhoMassLepGamma;   //!
   TBranch        *b_loosePhoMediumIDFunction;   //!
   TBranch        *b_loosePhoMediumIDPassHoverE;   //!
   TBranch        *b_loosePhoMediumIDPassSIEIE;   //!
   TBranch        *b_loosePhoMediumIDPassChIso;   //!
   TBranch        *b_loosePhoMediumIDPassNeuIso;   //!
   TBranch        *b_loosePhoMediumIDPassPhoIso;   //!
   TBranch        *b_nPhoNoID;   //!
   TBranch        *b_phoNoIDEt;   //!
   TBranch        *b_phoNoIDEta;   //!
   TBranch        *b_phoNoIDPhi;   //!
   TBranch        *b_phoNoIDIsBarrel;   //!
   TBranch        *b_phoNoIDHoverE;   //!
   TBranch        *b_phoNoIDSIEIE;   //!
   TBranch        *b_phoNoIDPFChIso;   //!
   TBranch        *b_phoNoIDTightID;   //!
   TBranch        *b_phoNoIDMediumID;   //!
   TBranch        *b_phoNoIDLooseID;   //!
   TBranch        *b_phoNoIDMVAId;   //!
   TBranch        *b_phoNoIDMVAId17v1;   //!
   TBranch        *b_phoNoIDGenMatchInd;   //!
   TBranch        *b_phoNoIDMassLepGamma;   //!
   TBranch        *b_phoNoIDMediumIDFunction;   //!
   TBranch        *b_phoNoIDMediumIDPassHoverE;   //!
   TBranch        *b_phoNoIDMediumIDPassSIEIE;   //!
   TBranch        *b_phoNoIDMediumIDPassChIso;   //!
   TBranch        *b_phoNoIDMediumIDPassNeuIso;   //!
   TBranch        *b_phoNoIDMediumIDPassPhoIso;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_elePFRelIso;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muPFRelIso;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_nfwdJet;   //!
   TBranch        *b_nBJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetDeepB;   //!
   TBranch        *b_jetGenJetIdx;   //!
   TBranch        *b_fwdJetPt;   //!
   TBranch        *b_fwdJetEta;   //!
   TBranch        *b_fwdJetPhi;   //!
   TBranch        *b_fwdJetMass;   //!
   TBranch        *b_dRPhotonJet;   //!
   TBranch        *b_dRPhotonLepton;   //!
   TBranch        *b_MPhotonLepton;   //!
   TBranch        *b_AnglePhotonLepton;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_genMass;   //!
   TBranch        *b_genStatus;   //!
   TBranch        *b_genStatusFlag;   //!
   TBranch        *b_genPDGID;   //!
   TBranch        *b_genMomIdx;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetMass;   //!
   TBranch        *b_genJetPartonFlavour;   //!
   TBranch        *b_M3;   //!
   TBranch        *b_M3_gamma;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_passPresel_Ele;   //!
   TBranch        *b_passPresel_Mu;   //!
   TBranch        *b_passAll_Ele;   //!
   TBranch        *b_passAll_Mu;   //!
   TBranch        *b_inHEMVeto;   //!
   TBranch        *b_photonIsGenuine;   //!
   TBranch        *b_photonIsMisIDEle;   //!
   TBranch        *b_photonIsHadronicPhoton;   //!
   TBranch        *b_photonIsHadronicFake;   //!
   TBranch        *b_loosePhotonIsGenuine;   //!
   TBranch        *b_loosePhotonIsMisIDEle;   //!
   TBranch        *b_loosePhotonIsHadronicPhoton;   //!
   TBranch        *b_loosePhotonIsHadronicFake;   //!
   TBranch        *b_photonNoIDIsGenuine;   //!
   TBranch        *b_photonNoIDIsMisIDEle;   //!
   TBranch        *b_photonNoIDIsHadronicPhoton;   //!
   TBranch        *b_photonNoIDIsHadronicFake;   //!
   TBranch        *b_photonParentage;   //!
   TBranch        *b_photonParentPID;   //!

   genloop(TTree *tree=0);
   virtual ~genloop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef genloop_cxx
genloop::genloop(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
      }
      f->GetObject("AnalysisTree",tree);

   }
   Init(tree);
}

genloop::~genloop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t genloop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t genloop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void genloop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   btagWeight = 0;
   btagWeight_b_Up = 0;
   btagWeight_b_Do = 0;
   btagWeight_l_Up = 0;
   btagWeight_l_Do = 0;
   btagSF = 0;
   phoEffWeight = 0;
   phoEffWeight_Id = 0;
   phoEffWeight_eVeto = 0;
   loosePhoEffWeight = 0;
   loosePhoEffWeight_Id = 0;
   loosePhoEffWeight_eVeto = 0;
   phoNoIDEffWeight = 0;
   phoNoIDEffWeight_Id = 0;
   phoNoIDEffWeight_eVeto = 0;
   phoEffWeight_Up = 0;
   phoEffWeight_Do = 0;
   phoEffWeight_Id_Up = 0;
   phoEffWeight_Id_Do = 0;
   phoEffWeight_eVeto_Up = 0;
   phoEffWeight_eVeto_Do = 0;
   loosePhoEffWeight_Up = 0;
   loosePhoEffWeight_Do = 0;
   loosePhoEffWeight_Id_Up = 0;
   loosePhoEffWeight_Id_Do = 0;
   loosePhoEffWeight_eVeto_Up = 0;
   loosePhoEffWeight_eVeto_Do = 0;
   phoNoIDEffWeight_Up = 0;
   phoNoIDEffWeight_Do = 0;
   phoNoIDEffWeight_Id_Up = 0;
   phoNoIDEffWeight_Id_Do = 0;
   phoNoIDEffWeight_eVeto_Up = 0;
   phoNoIDEffWeight_eVeto_Do = 0;
   genScaleSystWeights = 0;
   pdfSystWeight = 0;
   phoEt = 0;
   phoEta = 0;
   phoR9 = 0;
   phoPhi = 0;
   phoIsBarrel = 0;
   phoHoverE = 0;
   phoSIEIE = 0;
   phoPFChIso = 0;
   phoTightID = 0;
   phoMediumID = 0;
   phoGenMatchInd = 0;
   phoMassLepGamma = 0;
   loosePhoEt = 0;
   loosePhoEta = 0;
   loosePhoPhi = 0;
   loosePhoIsBarrel = 0;
   loosePhoHoverE = 0;
   loosePhoSIEIE = 0;
   loosePhoPFChIso = 0;
   loosePhoTightID = 0;
   loosePhoMediumID = 0;
   loosePhoLooseID = 0;
   loosePhoMVAId = 0;
   loosePhoMVAId17v1 = 0;
   loosePhoGenMatchInd = 0;
   loosePhoMassLepGamma = 0;
   loosePhoMediumIDFunction = 0;
   loosePhoMediumIDPassHoverE = 0;
   loosePhoMediumIDPassSIEIE = 0;
   loosePhoMediumIDPassChIso = 0;
   loosePhoMediumIDPassNeuIso = 0;
   loosePhoMediumIDPassPhoIso = 0;
   phoNoIDEt = 0;
   phoNoIDEta = 0;
   phoNoIDPhi = 0;
   phoNoIDIsBarrel = 0;
   phoNoIDHoverE = 0;
   phoNoIDSIEIE = 0;
   phoNoIDPFChIso = 0;
   phoNoIDTightID = 0;
   phoNoIDMediumID = 0;
   phoNoIDLooseID = 0;
   phoNoIDMVAId = 0;
   phoNoIDMVAId17v1 = 0;
   phoNoIDGenMatchInd = 0;
   phoNoIDMassLepGamma = 0;
   phoNoIDMediumIDFunction = 0;
   phoNoIDMediumIDPassHoverE = 0;
   phoNoIDMediumIDPassSIEIE = 0;
   phoNoIDMediumIDPassChIso = 0;
   phoNoIDMediumIDPassNeuIso = 0;
   phoNoIDMediumIDPassPhoIso = 0;
   elePt = 0;
   elePhi = 0;
   eleEta = 0;
   eleSCEta = 0;
   elePFRelIso = 0;
   muPt = 0;
   muEta = 0;
   muPhi = 0;
   muPFRelIso = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetDeepB = 0;
   jetGenJetIdx = 0;
   fwdJetPt = 0;
   fwdJetEta = 0;
   fwdJetPhi = 0;
   fwdJetMass = 0;
   dRPhotonJet = 0;
   dRPhotonLepton = 0;
   MPhotonLepton = 0;
   AnglePhotonLepton = 0;
   genPt = 0;
   genEta = 0;
   genPhi = 0;
   genMass = 0;
   genStatus = 0;
   genStatusFlag = 0;
   genPDGID = 0;
   genMomIdx = 0;
   genJetPt = 0;
   genJetEta = 0;
   genJetPhi = 0;
   genJetMass = 0;
   genJetPartonFlavour = 0;
   photonIsGenuine = 0;
   photonIsMisIDEle = 0;
   photonIsHadronicPhoton = 0;
   photonIsHadronicFake = 0;
   loosePhotonIsGenuine = 0;
   loosePhotonIsMisIDEle = 0;
   loosePhotonIsHadronicPhoton = 0;
   loosePhotonIsHadronicFake = 0;
   photonNoIDIsGenuine = 0;
   photonNoIDIsMisIDEle = 0;
   photonNoIDIsHadronicPhoton = 0;
   photonNoIDIsHadronicFake = 0;
   photonParentage = 0;
   photonParentPID = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweight_Up", &PUweight_Up, &b_PUweight_Up);
   fChain->SetBranchAddress("PUweight_Do", &PUweight_Do, &b_PUweight_Do);
   fChain->SetBranchAddress("prefireSF", &prefireSF, &b_prefireSF);
   fChain->SetBranchAddress("prefireSF_Up", &prefireSF_Up, &b_prefireSF_Up);
   fChain->SetBranchAddress("prefireSF_Do", &prefireSF_Do, &b_prefireSF_Do);
   fChain->SetBranchAddress("btagWeight", &btagWeight, &b_btagWeight);
   fChain->SetBranchAddress("btagWeight_1a", &btagWeight_1a, &b_btagWeight_1a);
   fChain->SetBranchAddress("btagWeight_b_Up", &btagWeight_b_Up, &b_btagWeight_b_Up);
   fChain->SetBranchAddress("btagWeight_b_Do", &btagWeight_b_Do, &b_btagWeight_b_Do);
   fChain->SetBranchAddress("btagWeight_l_Up", &btagWeight_l_Up, &b_btagWeight_l_Up);
   fChain->SetBranchAddress("btagWeight_l_Do", &btagWeight_l_Do, &b_btagWeight_l_Do);
   fChain->SetBranchAddress("btagWeight_1a_b_Up", &btagWeight_1a_b_Up, &b_btagWeight_1a_b_Up);
   fChain->SetBranchAddress("btagWeight_1a_b_Do", &btagWeight_1a_b_Do, &b_btagWeight_1a_b_Do);
   fChain->SetBranchAddress("btagWeight_1a_l_Up", &btagWeight_1a_l_Up, &b_btagWeight_1a_l_Up);
   fChain->SetBranchAddress("btagWeight_1a_l_Do", &btagWeight_1a_l_Do, &b_btagWeight_1a_l_Do);
   fChain->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
   fChain->SetBranchAddress("muEffWeight", &muEffWeight, &b_muEffWeight);
   fChain->SetBranchAddress("muEffWeight_IdIso", &muEffWeight_IdIso, &b_muEffWeight_IdIso);
   fChain->SetBranchAddress("muEffWeight_Trig", &muEffWeight_Trig, &b_muEffWeight_Trig);
   fChain->SetBranchAddress("muEffWeight_Up", &muEffWeight_Up, &b_muEffWeight_Up);
   fChain->SetBranchAddress("muEffWeight_Do", &muEffWeight_Do, &b_muEffWeight_Do);
   fChain->SetBranchAddress("muEffWeight_IdIso_Up", &muEffWeight_IdIso_Up, &b_muEffWeight_IdIso_Up);
   fChain->SetBranchAddress("muEffWeight_IdIso_Do", &muEffWeight_IdIso_Do, &b_muEffWeight_IdIso_Do);
   fChain->SetBranchAddress("muEffWeight_Trig_Up", &muEffWeight_Trig_Up, &b_muEffWeight_Trig_Up);
   fChain->SetBranchAddress("muEffWeight_Trig_Do", &muEffWeight_Trig_Do, &b_muEffWeight_Trig_Do);
   fChain->SetBranchAddress("eleEffWeight", &eleEffWeight, &b_eleEffWeight);
   fChain->SetBranchAddress("eleEffWeight_IdReco", &eleEffWeight_IdReco, &b_eleEffWeight_IdReco);
   fChain->SetBranchAddress("eleEffWeight_Trig", &eleEffWeight_Trig, &b_eleEffWeight_Trig);
   fChain->SetBranchAddress("eleEffWeight_Up", &eleEffWeight_Up, &b_eleEffWeight_Up);
   fChain->SetBranchAddress("eleEffWeight_Do", &eleEffWeight_Do, &b_eleEffWeight_Do);
   fChain->SetBranchAddress("eleEffWeight_IdReco_Up", &eleEffWeight_IdReco_Up, &b_eleEffWeight_IdReco_Up);
   fChain->SetBranchAddress("eleEffWeight_IdReco_Do", &eleEffWeight_IdReco_Do, &b_eleEffWeight_IdReco_Do);
   fChain->SetBranchAddress("eleEffWeight_Trig_Up", &eleEffWeight_Trig_Up, &b_eleEffWeight_Trig_Up);
   fChain->SetBranchAddress("eleEffWeight_Trig_Do", &eleEffWeight_Trig_Do, &b_eleEffWeight_Trig_Do);
   fChain->SetBranchAddress("phoEffWeight", &phoEffWeight, &b_phoEffWeight);
   fChain->SetBranchAddress("phoEffWeight_Id", &phoEffWeight_Id, &b_phoEffWeight_Id);
   fChain->SetBranchAddress("phoEffWeight_eVeto", &phoEffWeight_eVeto, &b_phoEffWeight_eVeto);
   fChain->SetBranchAddress("loosePhoEffWeight", &loosePhoEffWeight, &b_loosePhoEffWeight);
   fChain->SetBranchAddress("loosePhoEffWeight_Id", &loosePhoEffWeight_Id, &b_loosePhoEffWeight_Id);
   fChain->SetBranchAddress("loosePhoEffWeight_eVeto", &loosePhoEffWeight_eVeto, &b_loosePhoEffWeight_eVeto);
   fChain->SetBranchAddress("phoNoIDEffWeight", &phoNoIDEffWeight, &b_phoNoIDEffWeight);
   fChain->SetBranchAddress("phoNoIDEffWeight_Id", &phoNoIDEffWeight_Id, &b_phoNoIDEffWeight_Id);
   fChain->SetBranchAddress("phoNoIDEffWeight_eVeto", &phoNoIDEffWeight_eVeto, &b_phoNoIDEffWeight_eVeto);
   fChain->SetBranchAddress("phoEffWeight_Up", &phoEffWeight_Up, &b_phoEffWeight_Up);
   fChain->SetBranchAddress("phoEffWeight_Do", &phoEffWeight_Do, &b_phoEffWeight_Do);
   fChain->SetBranchAddress("phoEffWeight_Id_Up", &phoEffWeight_Id_Up, &b_phoEffWeight_Id_Up);
   fChain->SetBranchAddress("phoEffWeight_Id_Do", &phoEffWeight_Id_Do, &b_phoEffWeight_Id_Do);
   fChain->SetBranchAddress("phoEffWeight_eVeto_Up", &phoEffWeight_eVeto_Up, &b_phoEffWeight_eVeto_Up);
   fChain->SetBranchAddress("phoEffWeight_eVeto_Do", &phoEffWeight_eVeto_Do, &b_phoEffWeight_eVeto_Do);
   fChain->SetBranchAddress("loosePhoEffWeight_Up", &loosePhoEffWeight_Up, &b_loosePhoEffWeight_Up);
   fChain->SetBranchAddress("loosePhoEffWeight_Do", &loosePhoEffWeight_Do, &b_loosePhoEffWeight_Do);
   fChain->SetBranchAddress("loosePhoEffWeight_Id_Up", &loosePhoEffWeight_Id_Up, &b_loosePhoEffWeight_Id_Up);
   fChain->SetBranchAddress("loosePhoEffWeight_Id_Do", &loosePhoEffWeight_Id_Do, &b_loosePhoEffWeight_Id_Do);
   fChain->SetBranchAddress("loosePhoEffWeight_eVeto_Up", &loosePhoEffWeight_eVeto_Up, &b_loosePhoEffWeight_eVeto_Up);
   fChain->SetBranchAddress("loosePhoEffWeight_eVeto_Do", &loosePhoEffWeight_eVeto_Do, &b_loosePhoEffWeight_eVeto_Do);
   fChain->SetBranchAddress("phoNoIDEffWeight_Up", &phoNoIDEffWeight_Up, &b_phoNoIDEffWeight_Up);
   fChain->SetBranchAddress("phoNoIDEffWeight_Do", &phoNoIDEffWeight_Do, &b_phoNoIDEffWeight_Do);
   fChain->SetBranchAddress("phoNoIDEffWeight_Id_Up", &phoNoIDEffWeight_Id_Up, &b_phoNoIDEffWeight_Id_Up);
   fChain->SetBranchAddress("phoNoIDEffWeight_Id_Do", &phoNoIDEffWeight_Id_Do, &b_phoNoIDEffWeight_Id_Do);
   fChain->SetBranchAddress("phoNoIDEffWeight_eVeto_Up", &phoNoIDEffWeight_eVeto_Up, &b_phoNoIDEffWeight_eVeto_Up);
   fChain->SetBranchAddress("phoNoIDEffWeight_eVeto_Do", &phoNoIDEffWeight_eVeto_Do, &b_phoNoIDEffWeight_eVeto_Do);
   fChain->SetBranchAddress("q2weight_Up", &q2weight_Up, &b_q2weight_Up);
   fChain->SetBranchAddress("q2weight_Do", &q2weight_Do, &b_q2weight_Do);
   fChain->SetBranchAddress("q2weight_nominal", &q2weight_nominal, &b_q2weight_nominal);
   fChain->SetBranchAddress("genScaleSystWeights", &genScaleSystWeights, &b_genScaleSystWeights);
   fChain->SetBranchAddress("pdfWeight", &pdfWeight, &b_pdfWeight);
   fChain->SetBranchAddress("pdfuncer", &pdfuncer, &b_pdfuncer);
   fChain->SetBranchAddress("pdfweight_Up", &pdfweight_Up, &b_pdfweight_Up);
   fChain->SetBranchAddress("pdfweight_Do", &pdfweight_Do, &b_pdfweight_Do);
   fChain->SetBranchAddress("pdfSystWeight", &pdfSystWeight, &b_pdfSystWeight);
   fChain->SetBranchAddress("ISRweight_Up", &ISRweight_Up, &b_ISRweight_Up);
   fChain->SetBranchAddress("ISRweight_Do", &ISRweight_Do, &b_ISRweight_Do);
   fChain->SetBranchAddress("FSRweight_Up", &FSRweight_Up, &b_FSRweight_Up);
   fChain->SetBranchAddress("FSRweight_Do", &FSRweight_Do, &b_FSRweight_Do);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("nu_pz", &nu_pz, &b_nu_pz);
   fChain->SetBranchAddress("nu_pz_other", &nu_pz_other, &b_nu_pz_other);
   fChain->SetBranchAddress("WtransMass", &WtransMass, &b_WtransMass);
   fChain->SetBranchAddress("Mt_blgammaMET", &Mt_blgammaMET, &b_Mt_blgammaMET);
   fChain->SetBranchAddress("Mt_lgammaMET", &Mt_lgammaMET, &b_Mt_lgammaMET);
   fChain->SetBranchAddress("M_bjj", &M_bjj, &b_M_bjj);
   fChain->SetBranchAddress("M_bjjgamma", &M_bjjgamma, &b_M_bjjgamma);
   fChain->SetBranchAddress("M_jj", &M_jj, &b_M_jj);
   fChain->SetBranchAddress("MassCuts", &MassCuts, &b_MassCuts);
   fChain->SetBranchAddress("TopHad_pt", &TopHad_pt, &b_TopHad_pt);
   fChain->SetBranchAddress("TopHad_eta", &TopHad_eta, &b_TopHad_eta);
   fChain->SetBranchAddress("TopHad_phi", &TopHad_phi, &b_TopHad_phi);
   fChain->SetBranchAddress("TopLep_pt", &TopLep_pt, &b_TopLep_pt);
   fChain->SetBranchAddress("TopLep_eta", &TopLep_eta, &b_TopLep_eta);
   fChain->SetBranchAddress("TopLep_phi", &TopLep_phi, &b_TopLep_phi);
   fChain->SetBranchAddress("TopLep_charge", &TopLep_charge, &b_TopLep_charge);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("DiphoMass", &DiphoMass, &b_DiphoMass);
   fChain->SetBranchAddress("DilepMass", &DilepMass, &b_DilepMass);
   fChain->SetBranchAddress("DilepDelR", &DilepDelR, &b_DilepDelR);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("nPhoBarrel", &nPhoBarrel, &b_nPhoBarrel);
   fChain->SetBranchAddress("nPhoEndcap", &nPhoEndcap, &b_nPhoEndcap);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoIsBarrel", &phoIsBarrel, &b_phoIsBarrel);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSIEIE", &phoSIEIE, &b_phoSIEIE);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoTightID", &phoTightID, &b_phoTightID);
   fChain->SetBranchAddress("phoMediumID", &phoMediumID, &b_phoMediumID);
   fChain->SetBranchAddress("phoGenMatchInd", &phoGenMatchInd, &b_phoGenMatchInd);
   fChain->SetBranchAddress("phoMassLepGamma", &phoMassLepGamma, &b_phoMassLepGamma);
   fChain->SetBranchAddress("nLoosePho", &nLoosePho, &b_nLoosePho);
   fChain->SetBranchAddress("loosePhoEt", &loosePhoEt, &b_loosePhoEt);
   fChain->SetBranchAddress("loosePhoEta", &loosePhoEta, &b_loosePhoEta);
   fChain->SetBranchAddress("loosePhoPhi", &loosePhoPhi, &b_loosePhoPhi);
   fChain->SetBranchAddress("loosePhoIsBarrel", &loosePhoIsBarrel, &b_loosePhoIsBarrel);
   fChain->SetBranchAddress("loosePhoHoverE", &loosePhoHoverE, &b_loosePhoHoverE);
   fChain->SetBranchAddress("loosePhoSIEIE", &loosePhoSIEIE, &b_loosePhoSIEIE);
   fChain->SetBranchAddress("loosePhoPFChIso", &loosePhoPFChIso, &b_loosePhoPFChIso);
   fChain->SetBranchAddress("loosePhoTightID", &loosePhoTightID, &b_loosePhoTightID);
   fChain->SetBranchAddress("loosePhoMediumID", &loosePhoMediumID, &b_loosePhoMediumID);
   fChain->SetBranchAddress("loosePhoLooseID", &loosePhoLooseID, &b_loosePhoLooseID);
   fChain->SetBranchAddress("loosePhoMVAId", &loosePhoMVAId, &b_loosePhoMVAId);
   fChain->SetBranchAddress("loosePhoMVAId17v1", &loosePhoMVAId17v1, &b_loosePhoMVAId17v1);
   fChain->SetBranchAddress("loosePhoGenMatchInd", &loosePhoGenMatchInd, &b_loosePhoGenMatchInd);
   fChain->SetBranchAddress("loosePhoMassLepGamma", &loosePhoMassLepGamma, &b_loosePhoMassLepGamma);
   fChain->SetBranchAddress("loosePhoMediumIDFunction", &loosePhoMediumIDFunction, &b_loosePhoMediumIDFunction);
   fChain->SetBranchAddress("loosePhoMediumIDPassHoverE", &loosePhoMediumIDPassHoverE, &b_loosePhoMediumIDPassHoverE);
   fChain->SetBranchAddress("loosePhoMediumIDPassSIEIE", &loosePhoMediumIDPassSIEIE, &b_loosePhoMediumIDPassSIEIE);
   fChain->SetBranchAddress("loosePhoMediumIDPassChIso", &loosePhoMediumIDPassChIso, &b_loosePhoMediumIDPassChIso);
   fChain->SetBranchAddress("loosePhoMediumIDPassNeuIso", &loosePhoMediumIDPassNeuIso, &b_loosePhoMediumIDPassNeuIso);
   fChain->SetBranchAddress("loosePhoMediumIDPassPhoIso", &loosePhoMediumIDPassPhoIso, &b_loosePhoMediumIDPassPhoIso);
   fChain->SetBranchAddress("nPhoNoID", &nPhoNoID, &b_nPhoNoID);
   fChain->SetBranchAddress("phoNoIDEt", &phoNoIDEt, &b_phoNoIDEt);
   fChain->SetBranchAddress("phoNoIDEta", &phoNoIDEta, &b_phoNoIDEta);
   fChain->SetBranchAddress("phoNoIDPhi", &phoNoIDPhi, &b_phoNoIDPhi);
   fChain->SetBranchAddress("phoNoIDIsBarrel", &phoNoIDIsBarrel, &b_phoNoIDIsBarrel);
   fChain->SetBranchAddress("phoNoIDHoverE", &phoNoIDHoverE, &b_phoNoIDHoverE);
   fChain->SetBranchAddress("phoNoIDSIEIE", &phoNoIDSIEIE, &b_phoNoIDSIEIE);
   fChain->SetBranchAddress("phoNoIDPFChIso", &phoNoIDPFChIso, &b_phoNoIDPFChIso);
   fChain->SetBranchAddress("phoNoIDTightID", &phoNoIDTightID, &b_phoNoIDTightID);
   fChain->SetBranchAddress("phoNoIDMediumID", &phoNoIDMediumID, &b_phoNoIDMediumID);
   fChain->SetBranchAddress("phoNoIDLooseID", &phoNoIDLooseID, &b_phoNoIDLooseID);
   fChain->SetBranchAddress("phoNoIDMVAId", &phoNoIDMVAId, &b_phoNoIDMVAId);
   fChain->SetBranchAddress("phoNoIDMVAId17v1", &phoNoIDMVAId17v1, &b_phoNoIDMVAId17v1);
   fChain->SetBranchAddress("phoNoIDGenMatchInd", &phoNoIDGenMatchInd, &b_phoNoIDGenMatchInd);
   fChain->SetBranchAddress("phoNoIDMassLepGamma", &phoNoIDMassLepGamma, &b_phoNoIDMassLepGamma);
   fChain->SetBranchAddress("phoNoIDMediumIDFunction", &phoNoIDMediumIDFunction, &b_phoNoIDMediumIDFunction);
   fChain->SetBranchAddress("phoNoIDMediumIDPassHoverE", &phoNoIDMediumIDPassHoverE, &b_phoNoIDMediumIDPassHoverE);
   fChain->SetBranchAddress("phoNoIDMediumIDPassSIEIE", &phoNoIDMediumIDPassSIEIE, &b_phoNoIDMediumIDPassSIEIE);
   fChain->SetBranchAddress("phoNoIDMediumIDPassChIso", &phoNoIDMediumIDPassChIso, &b_phoNoIDMediumIDPassChIso);
   fChain->SetBranchAddress("phoNoIDMediumIDPassNeuIso", &phoNoIDMediumIDPassNeuIso, &b_phoNoIDMediumIDPassNeuIso);
   fChain->SetBranchAddress("phoNoIDMediumIDPassPhoIso", &phoNoIDMediumIDPassPhoIso, &b_phoNoIDMediumIDPassPhoIso);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("elePFRelIso", &elePFRelIso, &b_elePFRelIso);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muPFRelIso", &muPFRelIso, &b_muPFRelIso);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("nfwdJet", &nfwdJet, &b_nfwdJet);
   fChain->SetBranchAddress("nBJet", &nBJet, &b_nBJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetDeepB", &jetDeepB, &b_jetDeepB);
   fChain->SetBranchAddress("jetGenJetIdx", &jetGenJetIdx, &b_jetGenJetIdx);
   fChain->SetBranchAddress("fwdJetPt", &fwdJetPt, &b_fwdJetPt);
   fChain->SetBranchAddress("fwdJetEta", &fwdJetEta, &b_fwdJetEta);
   fChain->SetBranchAddress("fwdJetPhi", &fwdJetPhi, &b_fwdJetPhi);
   fChain->SetBranchAddress("fwdJetMass", &fwdJetMass, &b_fwdJetMass);
   fChain->SetBranchAddress("dRPhotonJet", &dRPhotonJet, &b_dRPhotonJet);
   fChain->SetBranchAddress("dRPhotonLepton", &dRPhotonLepton, &b_dRPhotonLepton);
   fChain->SetBranchAddress("MPhotonLepton", &MPhotonLepton, &b_MPhotonLepton);
   fChain->SetBranchAddress("AnglePhotonLepton", &AnglePhotonLepton, &b_AnglePhotonLepton);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("genPt", &genPt, &b_genPt);
   fChain->SetBranchAddress("genEta", &genEta, &b_genEta);
   fChain->SetBranchAddress("genPhi", &genPhi, &b_genPhi);
   fChain->SetBranchAddress("genMass", &genMass, &b_genMass);
   fChain->SetBranchAddress("genStatus", &genStatus, &b_genStatus);
   fChain->SetBranchAddress("genStatusFlag", &genStatusFlag, &b_genStatusFlag);
   fChain->SetBranchAddress("genPDGID", &genPDGID, &b_genPDGID);
   fChain->SetBranchAddress("genMomIdx", &genMomIdx, &b_genMomIdx);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetMass", &genJetMass, &b_genJetMass);
   fChain->SetBranchAddress("genJetPartonFlavour", &genJetPartonFlavour, &b_genJetPartonFlavour);
   fChain->SetBranchAddress("M3", &M3, &b_M3);
   fChain->SetBranchAddress("M3_gamma", &M3_gamma, &b_M3_gamma);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("passPresel_Ele", &passPresel_Ele, &b_passPresel_Ele);
   fChain->SetBranchAddress("passPresel_Mu", &passPresel_Mu, &b_passPresel_Mu);
   fChain->SetBranchAddress("passAll_Ele", &passAll_Ele, &b_passAll_Ele);
   fChain->SetBranchAddress("passAll_Mu", &passAll_Mu, &b_passAll_Mu);
   fChain->SetBranchAddress("inHEMVeto", &inHEMVeto, &b_inHEMVeto);
   fChain->SetBranchAddress("photonIsGenuine", &photonIsGenuine, &b_photonIsGenuine);
   fChain->SetBranchAddress("photonIsMisIDEle", &photonIsMisIDEle, &b_photonIsMisIDEle);
   fChain->SetBranchAddress("photonIsHadronicPhoton", &photonIsHadronicPhoton, &b_photonIsHadronicPhoton);
   fChain->SetBranchAddress("photonIsHadronicFake", &photonIsHadronicFake, &b_photonIsHadronicFake);
   fChain->SetBranchAddress("loosePhotonIsGenuine", &loosePhotonIsGenuine, &b_loosePhotonIsGenuine);
   fChain->SetBranchAddress("loosePhotonIsMisIDEle", &loosePhotonIsMisIDEle, &b_loosePhotonIsMisIDEle);
   fChain->SetBranchAddress("loosePhotonIsHadronicPhoton", &loosePhotonIsHadronicPhoton, &b_loosePhotonIsHadronicPhoton);
   fChain->SetBranchAddress("loosePhotonIsHadronicFake", &loosePhotonIsHadronicFake, &b_loosePhotonIsHadronicFake);
   fChain->SetBranchAddress("photonNoIDIsGenuine", &photonNoIDIsGenuine, &b_photonNoIDIsGenuine);
   fChain->SetBranchAddress("photonNoIDIsMisIDEle", &photonNoIDIsMisIDEle, &b_photonNoIDIsMisIDEle);
   fChain->SetBranchAddress("photonNoIDIsHadronicPhoton", &photonNoIDIsHadronicPhoton, &b_photonNoIDIsHadronicPhoton);
   fChain->SetBranchAddress("photonNoIDIsHadronicFake", &photonNoIDIsHadronicFake, &b_photonNoIDIsHadronicFake);
   fChain->SetBranchAddress("photonParentage", &photonParentage, &b_photonParentage);
   fChain->SetBranchAddress("photonParentPID", &photonParentPID, &b_photonParentPID);
   Notify();
}

Bool_t genloop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void genloop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t genloop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef genloop_cxx
