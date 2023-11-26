#ifndef GENRECOANALYZER_H
#define GENRECOANALYZER_H

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>

class GenRecoAnalyzer {
public:
    GenRecoAnalyzer(const char* inputFileName, const char* outputFileName);
    ~GenRecoAnalyzer();
    void RunAnalysis();

private:
    TChain* chain;
    TFile* outputFile;

    // Histograms
    TH1F* Rapidity_diff_Gen;
    TH1F* Rapidity_diff_Reco;
    TH2F* correlMatrix;
    TH2F* correlMatrix_reco;

    // Data members for tree branches
    std::vector<Float_t>* elePt;
    std::vector<Float_t>* muPt;  
    std::vector<Float_t>* jetPt;
    std::vector<Float_t>* phoEt; 
    std::vector<Float_t>* phoEta; 
    std::vector<Float_t>* phoPhi;
    std::vector<int>* genPDGID;
    std::vector<int>* genMomIdx;
    std::vector<Float_t>* genMass;
    std::vector<Float_t>* genPhi;
    std::vector<Float_t>* genEta;
    std::vector<Float_t>* genPt;  
    float evtWeight;
    int nEle;
    int nMu;
    int nGenJet;
    int nGenPart;
    int nJet;
    int nBJet;
    int nPho;
    float TopHad_pt;
    float TopLep_pt;
    float TopHad_eta;
    float TopLep_eta;
    float TopLep_phi;
    float TopHad_phi;
    float Mt_blgammaMET; 
    float M_bjj;         
    float TopLep_charge;
    float PUweight;
    float btagWeight_1a;
    float prefireSF;
    float muEffWeight;
    float eleEffWeight;
    bool passPresel_Ele;
    bool passPresel_Mu;
    int lumis;
    float weight_gen;
    float weight_reco;
    float rapidity_diff_gen;
    float rapidity_diff_reco;
    double weight;
    

    // Utility methods
    void InitializeHistograms();
    void ProcessEvent(int eventNumber);
    void SetupBranches();
    void SaveResults();

    // Helper methods for analysis
    void CalculateRapidityDifferences();
    void FillHistograms();
    void AnalyzeGenParticles();
    void AnalyzeRecoParticles();
    void CalculateCorrelations();

   
};

#endif // GENRECOANALYZER_H
