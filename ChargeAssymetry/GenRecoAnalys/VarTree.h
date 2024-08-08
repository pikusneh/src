#ifndef VARTREE_H
#define VARTREE_H

#include <vector>
#include "TChain.h"
#include <TLorentzVector.h>

// Declare your Event class here
class Event
{
public:
    Event(TChain *trev);
    ~Event();

    void Analyze();
    void LoopOverEvents();
    
    TChain *tr;
    std::vector<Float_t> *elePt = 0;
    std::vector<Float_t> *muPt = 0;
    std::vector<Float_t> *jetPt = 0;
    std::vector<Float_t> *phoEt = 0;
    std::vector<Float_t> *phoEta = 0;
    std::vector<Float_t> *phoPhi = 0;
    std::vector<int> *genPDGID = 0;
    std::vector<int> *genMomIdx = 0;
    std::vector<Float_t> *genMass = 0;
    std::vector<Float_t> *genPhi = 0;
    std::vector<Float_t> *genEta = 0;
    std::vector<Float_t> *genPt = 0;

    // Scalars
    float evtWeight = 0;
    int nEle = 0;
    int nMu = 0;
    int nGenJet = 0;
    int nJet = 0;
    int nBJet = 0;
    int nPho = 0;
    int lumis = 0;
    float btagWeight_1a = 0;
    float prefireSF = 0;
    float muEffWeight = 0;
    float eleEffWeight = 0;
    float TopHad_pt = 0;
    float TopLep_pt = 0;
    float TopHad_eta = 0;
    float TopLep_eta = 0;
    float TopLep_phi = 0;
    float TopHad_phi = 0;
    float Mt_blgammaMET = 0;
    float M_bjj = 0;
    float TopLep_charge = 0;
    float PUweight = 0;
    bool passPresel_Ele;
    bool passPresel_Mu;
    float Rapidity_T_ele = 0;
    float Rapidity_t_ele = 0;
    float Rapidity_T_mu = 0;
    float Rapidity_t_mu = 0;
    float PhoRapidity = 0;
    double weight = 0;
    double weight_gen = 0;
    float rapidity_diff_gen = 0;
    float event_weight_gen = 0;
    float rapidity_diff_reco = 0;
    float event_weight_reco = 0;

    // Add all other member variables here if necessary

private:
    
    void SetBranchStatuses();
    void SetBranchAddresses();
    void GetEntry(int i);
};

#endif //
