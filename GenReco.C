#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
// #include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"

#endif

void GenReco()
{

    TH1F *Rapidity_diff_Gen = new TH1F("Rapidity_diff_Gen", "Rapidity_Diff_Gen", 20, -3, 3);
    TH1F *Rapidity_diff_Reco = new TH1F("Rapidity_diff_Reco", "Rapidity_Diff_Reco", 20, -3, 3);
    TH2F *correlMatrix = new TH2F("correlMatrix", "correl Matrix", 30, -3, 3, 30, -3, 3);
    TH2F *correlMatrix_reco = new TH2F("correlMatrix_reco", "correl Matrix_reco", 20, -3, 3, 20, -3, 3);

    Rapidity_diff_Gen->SetOption("hist");
    Rapidity_diff_Reco->SetOption("hist");
    TTree *tree1 = new TTree("tree1", "tree1");

    TFile *fl = new TFile("GenLevel.root", "RECREATE");
    TChain *tr = new TChain("AnalysisTree");

    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");

    // RooUnfoldResponse response(20, -3, 3, 20, -3, 3);

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
    Int_t nGenJet = 0;
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
    tr->SetBranchStatus("nGenJet", 1);
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
    cout << "nGenPart" << nGenPart << std::endl;

    /////////////////////////////////////////////////////////////Set Branch Address/////////////////////////////////////////////////////////////////
    tr->SetBranchAddress("elePt", &elePt);
    tr->SetBranchAddress("muPt", &muPt);
    tr->SetBranchAddress("jetPt", &jetPt);
    tr->SetBranchAddress("phoEt", &phoEt);
    tr->SetBranchAddress("phoEta", &phoEta);
    tr->SetBranchAddress("phoPhi", &phoPhi);
    tr->SetBranchAddress("nEle", &nEle);
    tr->SetBranchAddress("nMu", &nMu);
    tr->SetBranchAddress("nGenJet", &nGenJet);
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
    float event_weight_gen = 0;
    float rapidity_diff_reco = 0;
    float event_weight_reco = 0;
    tree1->Branch("rapidity_diff_gen", &rapidity_diff_gen);
    tree1->Branch("rapidity_diff_reco", &rapidity_diff_reco);
    tree1->Branch("event_weight_gen", &event_weight_gen);
    tree1->Branch("event_weight_reco", &event_weight_reco);

    // numevents = 10000;

    float Rapidity_T_gen, Rapidity_t_gen;

    for (int j = 0; j < numevents; j++)
    {
        tr->GetEntry(j);
        int topindex = -1;
        int antitopindex = -1;
        int photonindex = -1;
        int electronindex = -1;
        bool passPresel_Ele_gen = false;

        weight_gen = evtWeight;
        weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;
        TLorentzVector Top_gen;
        TLorentzVector Anti_Top_gen;
        TLorentzVector Top;
        TLorentzVector AntiTop;
        TLorentzVector Photon_gen;
        TLorentzVector Ele_gen;
        int totalPhotons = 0;
        int totalbjet = 0;

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
            if (genPDGID->at(i) == 11)
            {
                electronindex = i;
            }
            Ele_gen.SetPtEtaPhiM((*genPt)[electronindex], (*genEta)[electronindex], (*genPhi)[electronindex], (*genMass)[electronindex]);
            float EleRapidity_gen = Ele_gen.Rapidity();
            float YEle_gen = TMath::Abs(EleRapidity_gen);
            float Ele_Pt = Ele_gen.Pt();

            if (Ele_Pt >= 35.0 && YEle_gen <= 2.4)
            {
                passPresel_Ele_gen = true;
            }

            if (genPDGID->at(i) == 22)
            {
                photonindex = i;
                totalPhotons++;
                // std::cout << "photons" << totalPhotons << std::endl;
            }

            if ((genPDGID->at(i) == 5 || genPDGID->at(i) == -5) && genPt->at(i) >= 30.0 && fabs(genEta->at(i)) <= 2.4)
            {
                totalbjet++;
            }
        }

        // correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
        // }

        if (passPresel_Ele_gen == true && nGenJet >= 4 && totalbjet >= 1 && totalPhotons == 1)
        {
            Photon_gen.SetPtEtaPhiM((*genPt)[photonindex], (*genEta)[photonindex], (*genPhi)[photonindex], (*genMass)[photonindex]);
            float PhoRapidity_gen = Photon_gen.Rapidity();
            float YPho_gen = TMath::Abs(PhoRapidity_gen);
            float Pho_Pt = Photon_gen.Pt();
            if (YPho_gen >= 1)
            {
                Top_gen.SetPtEtaPhiM((*genPt)[topindex], (*genEta)[topindex], (*genPhi)[topindex], (*genMass)[topindex]);
                Anti_Top_gen.SetPtEtaPhiM((*genPt)[antitopindex], (*genEta)[antitopindex], (*genPhi)[antitopindex], (*genMass)[antitopindex]);
                TLorentzVector TopAntiTop_gen = Top_gen + Anti_Top_gen;

                ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
                Rapidity_T_gen = Top_gen.Rapidity();
                Rapidity_t_gen = Anti_Top_gen.Rapidity();
                float YT_gen = TMath::Abs(Rapidity_T_gen); // top rapidity
                float Yt_gen = TMath::Abs(Rapidity_t_gen); // antitop rapidity
                rapidity_diff_gen = YT_gen - Yt_gen;
                event_weight_gen = weight_gen;
                Rapidity_diff_Gen->Fill(rapidity_diff_gen, event_weight_gen);
                Rapidity_diff_Gen->Scale(1.0 / Rapidity_diff_Gen->Integral());
                // std::cout << "Gen Level Rapidity Difference: " << rapidity_diff_gen << std::endl;
            }
        }
        if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1)
        {
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
                event_weight_reco = weight;
                rapidity_diff_reco = rapidity_diff_ele;
                Rapidity_diff_Reco->Fill(rapidity_diff_ele, event_weight_reco);
                Rapidity_diff_Reco->Scale(1.0 / Rapidity_diff_Reco->Integral());
                correlMatrix_reco->Fill(rapidity_diff_reco, rapidity_diff_reco);
                if (passPresel_Ele_gen && nGenJet >= 4 && totalbjet >= 1 && totalPhotons == 1)
                {
                    // response.Fill(rapidity_diff_reco, rapidity_diff_gen, event_weight_reco);
                    correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
                }
                // std::cout << " Reco Level Rapidity Difference: " << rapidity_diff_reco << std::endl;
                // correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
                // response.Fill(rapidity_diff_reco, rapidity_diff_gen, event_weight_reco);

                tree1->Fill();
            }
        }

        // correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
        // }

        //     correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
        // }
        // correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
        // correlMatrix->Scale(1.0 / correlMatrix->Integral());

        // std::cout << " Reco Level Rapidity Difference2: " << rapidity_diff_reco << std::endl;
        // std::cout << "Gen Level Rapidity Difference2: " << rapidity_diff_gen << std::endl;
    }

    fl->Close();

    TFile *file1 = new TFile("output_GenReco.root", "RECREATE");

    TCanvas *c7 = new TCanvas();
    c7->cd();
    Rapidity_diff_Gen->Draw();
    c7->Update();
    c7->SaveAs("/eos/user/s/ssnehshu/GenLevel/Rapidity_Gen.png");
    Rapidity_diff_Gen->Write();

    TCanvas *canvas2 = new TCanvas();
    canvas2->cd();
    Rapidity_diff_Reco->Draw();
    canvas2->Update();
    canvas2->SaveAs("/eos/user/s/ssnehshu/GenLevel/Rapidity_Reco.png");
    Rapidity_diff_Reco->Write();

    TCanvas *canvas = new TCanvas("canvas", "correl Matrix Canvas", 800, 600);
    correlMatrix->GetXaxis()->SetTitle("Gen Level");
    correlMatrix->GetYaxis()->SetTitle("Reco Level");
    correlMatrix->SetTitle("correl Matrix");
    correlMatrix->Draw("COLZ TEXT");
    canvas->Write();

    TCanvas *canvas3 = new TCanvas("canvas3", "correl Matrix Canvas", 800, 600);

    correlMatrix_reco->GetXaxis()->SetTitle("Reco Level");
    correlMatrix_reco->GetYaxis()->SetTitle("Reco Level");
    correlMatrix_reco->SetTitle("correl Matrix_reco");
    correlMatrix_reco->Draw("COLZTEXT");
    canvas3->SaveAs("correlMatrix_reco.png");
    canvas3->Write();

    // response.Write();

    tree1->Write();

    file1->Close();
}
#ifndef __CINT__
int main()
{
    GenReco();
    return 0;
}
#endif
