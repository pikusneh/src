#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"

void gen()
{
    // gROOT->SetBatch(true);
    //  gROOT->Reset();
    TH1F *h1 = new TH1F("h1", "toppt", 100, 0, 400);
    TH1F *h2 = new TH1F("h2", "topantitop_pt", 100, 0, 400);
    TH1F *h3 = new TH1F("h3", "topantitop_mass", 100, 200, 1200);
    TH1F *h4 = new TH1F("h4", "Ac vs topantitop_mass", 4, 300, 1100);
    TH1F *h5 = new TH1F("h5", "Top Rapidity", 6, -3, 3);
    TH1F *h6 = new TH1F("h6", "AntiTop Rapidity", 6, -3, 3);
    TH1F *h7 = new TH1F("h7", "Rapidity_Gen", 4, -3, 3);

    h7->SetOption("hist");

    TFile *fl = new TFile("GenLevel.root", "RECREATE");
    TChain *tr = new TChain("AnalysisTree");

    TFile *file1 = new TFile("output_GenLevel.root", "RECREATE");
    TTree *tree1 = new TTree("tree1", "tree1");

    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
    tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");

    int numevents = tr->GetEntries();

    int nGenPart;
    std::vector<int> *genPDGID = 0;
    std::vector<int> *genMomIdx = 0;
    std::vector<Float_t> *genMass = 0;
    std::vector<Float_t> *genPhi = 0;
    std::vector<Float_t> *genEta = 0;
    std::vector<Float_t> *genPt = 0;
    float evtWeight = 0;
    double weight = 0;

    tr->SetBranchStatus("*", 0);

    tr->SetBranchStatus("nGenPart", 1);
    tr->SetBranchStatus("genPDGID", 1);
    tr->SetBranchStatus("genMomIdx", 1);
    tr->SetBranchStatus("genMass", 1);
    tr->SetBranchStatus("genPhi", 1);
    tr->SetBranchStatus("genEta", 1);
    tr->SetBranchStatus("genPt", 1);
    tr->SetBranchStatus("evtWeight", 1);

    tr->SetBranchAddress("nGenPart", &nGenPart);
    tr->SetBranchAddress("genPDGID", &genPDGID);
    tr->SetBranchAddress("genMomIdx", &genMomIdx);
    tr->SetBranchAddress("genMass", &genMass);
    tr->SetBranchAddress("genPhi", &genPhi);
    tr->SetBranchAddress("genEta", &genEta);
    tr->SetBranchAddress("genPt", &genPt);
    tr->SetBranchAddress("evtWeight", &evtWeight);

    // T is for top
    // t is for anti top

    float Mass_Top, Mass_AntiTop, Rapidity_T, Rapidity_t, Pt_top, Pt_antitop, Pt_topantitop;
    float NetMass;
    float NetMassRange[] = {300, 400, 500, 600, 700, 800, 900};
    int numMassRanges = 7;
    float sumArr[numMassRanges];
    float diffArr[numMassRanges];
    float deltaAcArr[numMassRanges];
    float N_plusArr[numMassRanges];
    float N_minusArr[numMassRanges];

    for (int j = 0; j < numevents; j++)
    {
        int topindex = -1;
        int antitopindex = -1;
        tr->GetEntry(j);
        weight = evtWeight;
        TLorentzVector Top;
        TLorentzVector AntiTop;

        // std::cout << "nGenPart: genpdgidSize" << nGenPart<<"     " <<genPDGID->size() <<std::endl;
        for (int i = 0; i < nGenPart; i++)
        {
            //        std::cout << "PDGID[" << i << "]: " << (*genPDGID)[i] << std::endl;
            if (genPDGID->at(i) == 6)
            {
                topindex = i;
            }
            if (genPDGID->at(i) == -6)
            {
                antitopindex = i;
            }
        }

        Top.SetPtEtaPhiM((*genPt)[topindex], (*genEta)[topindex], (*genPhi)[topindex], (*genMass)[topindex]);
        AntiTop.SetPtEtaPhiM((*genPt)[antitopindex], (*genEta)[antitopindex], (*genPhi)[antitopindex], (*genMass)[antitopindex]);
        TLorentzVector TopAntiTop = Top + AntiTop;

        Mass_Top = Top.M();
        Mass_AntiTop = AntiTop.M();
        Pt_top = Top.Pt();
        Pt_antitop = AntiTop.Pt();
        Pt_topantitop = TopAntiTop.Pt();
        NetMass = TopAntiTop.M();

        ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
        Rapidity_T = Top.Y();
        Rapidity_t = AntiTop.Y();
        float YT = TMath::Abs(Rapidity_T); // top rapidity
        float Yt = TMath::Abs(Rapidity_t); // antitop rapidity
        h5->Fill(YT, weight);
        h6->Fill(Yt, weight);
        float rapidity_diff = YT - Yt;
        h7->Fill(rapidity_diff, weight);
        // std::cout<<"rapidity_diff = "<< rapidity_diff <<std::endl;

        // Initialize arrays to store sums and differences for each mass range

        for (int k = 0; k < numMassRanges; k++)
        {
            sumArr[k] = 0;
            diffArr[k] = 0;
            deltaAcArr[k] = 0;
            N_plusArr[k] = 0;
            N_minusArr[k] = 0;
            if (NetMass >= NetMassRange[k] && NetMass < NetMassRange[k + 1])
            {
                if (rapidity_diff > 0)
                {
                    N_plusArr[k] += weight;
                }
                else if (rapidity_diff < 0)
                {
                    N_minusArr[k] += weight;
                }
            }
        }

        for (int l = 0; l < numMassRanges; l++)
        {
            sumArr[l] = N_plusArr[l] + N_minusArr[l];
            diffArr[l] = N_plusArr[l] - N_minusArr[l];

            if (sumArr[l] != 0)
            {
                float Ac = diffArr[l] / sumArr[l];
                float Delta_diff = sqrt(sumArr[l]);
                float Delta_sum = sqrt(sumArr[l]);
                deltaAcArr[l] = Ac * sqrt((Delta_diff / diffArr[l]) * (Delta_diff / diffArr[l]) + (Delta_sum / sumArr[l]) * (Delta_sum / sumArr[l]));

                int bin = h4->FindBin(NetMassRange[l]);
                h4->SetBinContent(bin, Ac);
                h4->SetBinError(bin, deltaAcArr[l]);
                std::cout << "Ac" << l + 1 << " = " << Ac << std::endl;
            }
        }
    }
    TCanvas *c4 = new TCanvas();
    c4->cd();
    h4->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
    h4->GetXaxis()->SetTitleSize(0.04);
    h4->GetXaxis()->SetLabelSize(0.03);
    h4->GetYaxis()->SetTitle("A_{c}");
    h4->GetYaxis()->SetTitleOffset(1.1);
    h4->GetYaxis()->SetTitleSize(0.04);
    h4->GetYaxis()->SetLabelSize(0.03);

    h4->SetMarkerSize(1);
    h4->SetMarkerColor(kRed);
    h4->SetLineWidth(1);

    h4->Draw("E1");
    c4->Update();
    c4->SaveAs("/eos/user/s/ssnehshu/GenLevel/A_c vs ttbarmass.png");
}