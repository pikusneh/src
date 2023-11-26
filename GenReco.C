#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
// #include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"

void GenReco()
{

    TH1F *Rapidity_diff_Gen = new TH1F("Rapidity_diff_Gen", "Rapidity_Diff_Gen", 20, -3, 3);
    TH2F *responseMatrix = new TH2F("responseMatrix", "Response Matrix", 20, -3, 3, 20, -3, 3);

    Rapidity_diff_Gen->SetOption("hist");

    TFile *fl = new TFile("GenLevel.root", "RECREATE");
    TChain *tr = new TChain("AnalysisTree");

    TFile *file1 = new TFile("output_GenReco.root", "RECREATE");
    TTree *tree1 = new TTree("tree1", "tree1");

    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");

    //  Define the variables
    std::vector<Float_t> *elePt = 0;
    std::vector<Float_t> *muPt = 0;
    std::vector<Float_t> *jetPt = 0;
    std::vector<Float_t> *phoEt = 0;
    std::vector<Float_t> *phoEta = 0;
    std::vector<Float_t> *phoPhi = 0;

    int numevents = tr->GetEntries();

    int nGenPart;
    std::vector<int> *genPDGID = 0;
    std::vector<int> *genMomIdx = 0;
    std::vector<Float_t> *genMass = 0;
    std::vector<Float_t> *genPhi = 0;
    std::vector<Float_t> *genEta = 0;
    std::vector<Float_t> *genPt = 0;
    float evtWeight = 0;
    int nEle = 0;
    int nMu = 0;
    Int_t nJet = 0;
    Int_t nBJet = 0;
    Int_t nPho = 0;
    Int_t lumis = 0;
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
    bool passPresel_Ele;
    bool passPresel_Mu;
    Float_t Rapidity_T_ele; // taking T for top quark
    Float_t Rapidity_t_ele; // taking t for antitop quark
    Float_t Rapidity_T_mu;
    Float_t Rapidity_t_mu;

    Float_t PhoRapidity;

    ///////////////////////////////////////////////////////////////Set Branch Status//////////////////////////////////////////////////////////////////

    tr->SetBranchStatus("*", 0);
    tr->SetBranchStatus("elePt", 1);
    tr->SetBranchStatus("muPt", 1);
    tr->SetBranchStatus("jetPt", 1);
    tr->SetBranchStatus("phoEt", 1);
    tr->SetBranchStatus("phoEta", 1);
    tr->SetBranchStatus("phoPhi", 1);
    tr->SetBranchStatus("nEle", 1);
    tr->SetBranchStatus("nMu", 1);

    tr->SetBranchStatus("nGenPart", 1);
    tr->SetBranchStatus("genPDGID", 1);
    tr->SetBranchStatus("genMomIdx", 1);
    tr->SetBranchStatus("genMass", 1);
    tr->SetBranchStatus("genPhi", 1);
    tr->SetBranchStatus("genEta", 1);
    tr->SetBranchStatus("genPt", 1);
    tr->SetBranchStatus("evtWeight", 1);
    tr->SetBranchStatus("nJet", 1);
    tr->SetBranchStatus("nBJet", 1);
    tr->SetBranchStatus("nPho", 1);
    tr->SetBranchStatus("TopHad_pt", 1);
    tr->SetBranchStatus("TopLep_pt", 1);
    tr->SetBranchStatus("TopHad_eta", 1);
    tr->SetBranchStatus("TopLep_eta", 1);
    tr->SetBranchStatus("TopHad_phi", 1);
    tr->SetBranchStatus("TopLep_phi", 1);
    tr->SetBranchStatus("M_bjj", 1);
    tr->SetBranchStatus("Mt_blgammaMET", 1);
    tr->SetBranchStatus("TopLep_charge", 1);
    tr->SetBranchStatus("passPresel_Ele", 1);
    tr->SetBranchStatus("passPresel_Mu", 1);
    tr->SetBranchStatus("lumis", 1);
    tr->SetBranchStatus("evtWeight", 1);
    tr->SetBranchStatus("btagWeight_1a", 1);
    tr->SetBranchStatus("prefireSF", 1);
    tr->SetBranchStatus("muEffWeight", 1);
    tr->SetBranchStatus("eleEffWeight", 1);
    tr->SetBranchStatus("PUweight", 1);
    cout << "numevents" << numevents << std::endl;

    /////////////////////////////////////////////////////////////Set Branch Address/////////////////////////////////////////////////////////////////
    tr->SetBranchAddress("elePt", &elePt);
    tr->SetBranchAddress("muPt", &muPt);
    tr->SetBranchAddress("jetPt", &jetPt);
    tr->SetBranchAddress("phoEt", &phoEt);
    tr->SetBranchAddress("phoEta", &phoEta);
    tr->SetBranchAddress("phoPhi", &phoPhi);
    tr->SetBranchAddress("nEle", &nEle);
    tr->SetBranchAddress("nMu", &nMu);
    tr->SetBranchAddress("nGenPart", &nGenPart);
    tr->SetBranchAddress("genPDGID", &genPDGID);
    tr->SetBranchAddress("genMomIdx", &genMomIdx);
    tr->SetBranchAddress("genMass", &genMass);
    tr->SetBranchAddress("genPhi", &genPhi);
    tr->SetBranchAddress("genEta", &genEta);
    tr->SetBranchAddress("genPt", &genPt);
    tr->SetBranchAddress("evtWeight", &evtWeight);
    tr->SetBranchAddress("nJet", &nJet);
    tr->SetBranchAddress("nBJet", &nBJet);
    tr->SetBranchAddress("nPho", &nPho);
    tr->SetBranchAddress("Mt_blgammaMET", &Mt_blgammaMET); // toplepmass
    tr->SetBranchAddress("M_bjj", &M_bjj);                 // tophadmass
    tr->SetBranchAddress("TopHad_pt", &TopHad_pt);
    tr->SetBranchAddress("TopLep_pt", &TopLep_pt);
    tr->SetBranchAddress("TopHad_eta", &TopHad_eta);
    tr->SetBranchAddress("TopLep_eta", &TopLep_eta);
    tr->SetBranchAddress("TopHad_phi", &TopHad_phi);
    tr->SetBranchAddress("TopLep_phi", &TopLep_phi);
    tr->SetBranchAddress("TopLep_charge", &TopLep_charge);
    tr->SetBranchAddress("passPresel_Ele", &passPresel_Ele);
    tr->SetBranchAddress("passPresel_Mu", &passPresel_Mu);
    tr->SetBranchAddress("lumis", &lumis);
    tr->SetBranchAddress("btagWeight_1a", &btagWeight_1a);
    tr->SetBranchAddress("evtWeight", &evtWeight);
    tr->SetBranchAddress("prefireSF", &prefireSF);
    tr->SetBranchAddress("muEffWeight", &muEffWeight);
    tr->SetBranchAddress("eleEffWeight", &eleEffWeight);
    tr->SetBranchAddress("PUweight", &PUweight);

    double weight = 0;
    double weight_gen = 0;
    float rapidity_diff_gen = 0;
    float event_weight_gen;
    float rapidity_diff_reco = 0;
    float event_weight_reco;
    tree1->Branch("rapidity_diff_gen", &rapidity_diff_gen, "Rapidity_diff_gen/F");
    tree1->Branch("event_weight_gen", &event_weight_gen, "event_weight_gen/F");
    tree1->Branch("rapidity_diff_reco", &rapidity_diff_reco, "Rapidity_diff_reco/F");
    tree1->Branch("event_weight_reco", &event_weight_reco, "event_weight_reco/F");

    // T is for top
    // t is for anti top

    float Rapidity_T_gen, Rapidity_t_gen;

    for (int j = 0; j < numevents; j++)
    {

        int topindex;
        int antitopindex;
        tr->GetEntry(j);
        weight_gen = evtWeight;
        weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;
        TLorentzVector Top_gen;
        TLorentzVector Anti_Top_gen;
        TLorentzVector Top;
        TLorentzVector AntiTop;
        if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1)
        {
            for (int i = 0; i < nGenPart; i++)
            {
                if (genPDGID->at(i) == 6)
                {
                    topindex = i;
                }
                if (genPDGID->at(i) == -6)
                {
                    antitopindex = i;
                }
            }

            Top_gen.SetPtEtaPhiM((*genPt)[topindex], (*genEta)[topindex], (*genPhi)[topindex], (*genMass)[topindex]);
            Anti_Top_gen.SetPtEtaPhiM((*genPt)[antitopindex], (*genEta)[antitopindex], (*genPhi)[antitopindex], (*genMass)[antitopindex]);
            TLorentzVector TopAntiTop_gen = Top_gen + Anti_Top_gen;

            ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
            Rapidity_T_gen = Top_gen.Rapidity();
            Rapidity_t_gen = Anti_Top_gen.Rapidity();
            float YT_gen = TMath::Abs(Rapidity_T_gen); // top rapidity
            float Yt_gen = TMath::Abs(Rapidity_t_gen); // antitop rapidity
            float rapidity_diff_gen = YT_gen - Yt_gen;
            event_weight_gen = weight_gen;
            Rapidity_diff_Gen->Fill(rapidity_diff_gen, weight);

            // Rapidity_diff_Gen->Scale(1.0 / Rapidity_diff_Gen->Integral());

            ////////////////////////////////////////////////////////////////////////////////reco level/////////////////////////////////////////////////////////////
            //////////////Tlorentz for photon///////////////////////
            TLorentzVector Photon;
            Photon.SetPtEtaPhiM(phoEt->at(0), phoEta->at(0), phoPhi->at(0), 0.0);
            PhoRapidity = Photon.Rapidity();
            float YPho = TMath::Abs(PhoRapidity);
            float PhoPt_ele = Photon.Pt();
            if (YPho >= 1)
            {
                if (TopLep_charge > 0)
                {
                    Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
                    AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
                }
                else if (TopLep_charge < 0)
                {
                    AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
                    Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
                }

                ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
                Rapidity_T_ele = Top.Rapidity();
                Rapidity_t_ele = AntiTop.Rapidity();
                // absolute value of rapidity
                float YT_ele = TMath::Abs(Rapidity_T_ele); // top
                float Yt_ele = TMath::Abs(Rapidity_t_ele); // antitop
                float rapidity_diff_ele = YT_ele - Yt_ele; // difference of rapidity; delta y

                rapidity_diff_reco = rapidity_diff_ele;
                event_weight_reco = weight;
                responseMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
                tree1->Fill();
            }
        }
    }
    TCanvas *c7 = new TCanvas();

    c7->cd();
    Rapidity_diff_Gen->Draw();
    c7->Update();
    c7->SaveAs("/eos/user/s/ssnehshu/GenLevel/Rapidity_Gen.png");
    Rapidity_diff_Gen->Write();

    TCanvas *canvas = new TCanvas("canvas", "Response Matrix Canvas", 800, 600);
    responseMatrix->Draw("COLZTEXT");
    responseMatrix->GetXaxis()->SetTitle("Gen Level");
    responseMatrix->GetYaxis()->SetTitle("Reco Level");
    responseMatrix->SetTitle("Response Matrix");
    canvas->Update();
    canvas->SaveAs("responseMatrix.png");
    responseMatrix->Write();
    tree1->Write();
    fl->Close();
    file1->Close();
}