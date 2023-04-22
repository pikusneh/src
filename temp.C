#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
// use namespace std;
// gROOT.SetBatch(true);

void TT()
{
    int mode = 1;
    std::vector<TH1F *> hists;
    if (mode == 1)
    {
        const char *files[] = {"/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root"};
        Int_t colors[] = {kRed,
                          kCyan,
                          kBlue,
                          kYellow};
        // saving histogram in root file
        TFile *fl = new TFile("TTGamma.root", "RECREATE");
        char *event = "TTGamma";
    }
    else if (mode == 2)
    {
        const char *files[] = {"/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Dilepton_2016_AnalysisNtuple_1of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Dilepton_2016_AnalysisNtuple_2of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Dilepton_2016_AnalysisNtuple_3of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Dilepton_2016_AnalysisNtuple_4of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Dilepton_2016_AnalysisNtuple_5of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Hadronic_2016_AnalysisNtuple.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Semilept_2016_AnalysisNtuple_1of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Semilept_2016_AnalysisNtuple_2of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Semilept_2016_AnalysisNtuple_3of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Semilept_2016_AnalysisNtuple_4of5.root",
                               "/eos/user/s/ssnehshu/topquarksample/2016/TTbarPowheg_Semilept_2016_AnalysisNtuple_5of5.root"};
        Int_t colors[] = {kOrange, kGreen, kAzure, kViolet};
        TFile *fl = new TFile("TTbar.root", "RECREATE");
        char *event = "TTbar";
    }
    gROOT->SetBatch(true);
    // defining histograms

    TH1F *h1 = new TH1F("h1", "elept", 100, 0, 400);
    hists.push_back(h1);
    TH1F *h2 = new TH1F("h2", "muonpt", 100, 0, 400);
    hists.push_back(h2);
    TH1F *h3 = new TH1F("h3", "jetspt", 100, 0, 500);
    hists.push_back(h3);
    // TH1F *h4 = new TH1F("h4","tophadpt",100,-11000,2000);hists.push_back(h4);
    // TH1F *h5 = new TH1F("h5","topleppt",100,-11000,2000);hists.push_back(h5);
    TH1F *h6 = new TH1F("h6", "tophadeta", 100, -11000, 2000);
    hists.push_back(h6);
    TH1F *h7 = new TH1F("h7", "toplepeta", 100, -11000, 2000);
    hists.push_back(h7);
    TH1F *h8 = new TH1F("h8", "toplepcharge", 100, -11000, 2000);
    hists.push_back(h8);
    TH1F *h9 = new TH1F("h9", "top", 100, -3, 3);
    hists.push_back(h9);
    TH1F *h10 = new TH1F("h10", "antitop", 100, -3, 3);
    hists.push_back(h10);
    TH1F *h9ele = new TH1F("h9ele", "top", 50, -3, 3);
    hists.push_back(h9ele);
    TH1F *h10ele = new TH1F("h10ele", "antitop", 50, -3, 3);
    hists.push_back(h10ele);
    TH1F *h9mu = new TH1F("h9mu", "top", 50, -3, 3);
    hists.push_back(h9mu);
    TH1F *h10mu = new TH1F("h10mu", "antitop", 50, -3, 3);
    hists.push_back(h10mu);
    TH1F *h11phi = new TH1F("h11", "topphi", 50, -3, 3);
    hists.push_back(h11phi);
    TH1F *h12phi = new TH1F("h12", "antitopphi", 50, -3, 3);
    hists.push_back(h12phi);

    h9ele->SetOption("hist");
    h9mu->SetOption("hist");
    h10ele->SetOption("hist");
    h10mu->SetOption("hist");
    h11phi->SetOption("hist");
    h12phi->SetOption("hist");

    h9ele->SetLineColor(kBlack);
    h9mu->SetLineColor(kBlack);
    h10ele->SetLineColor(kBlack);
    h10mu->SetLineColor(kBlack);

    h9ele->SetFillColor(colors[0]);
    h9mu->SetFillColor(colors[1]);
    h10ele->SetFillColor(colors[2]);
    h10mu->SetFillColor(colors[3]);

    TChain *tr = new TChain("AnalysisTree");

    for (int i = 0; i < sizeof(files); ++i)
    {
        tr->Add(files[i]);
    }

    // tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
    // tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
    // tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");
    // TTree *tr = (TTree*)fl-> Get("AnalysisTree");
    Int_t numevents = tr->GetEntries();
    // numevents = 1000;
    // Define the variables
    std::vector<Float_t> *elePt = 0;
    std::vector<Float_t> *muPt = 0;
    std::vector<Float_t> *jetPt = 0;
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
    Float_t TopLep_mass = 0;
    Float_t TopHad_mass = 0;
    Float_t TopLep_charge = 0;
    Float_t PUweight = 0;
    bool passPresel_Ele;
    bool passPresel_Mu;

    // Set Branch Status

    tr->SetBranchStatus("*", 0);
    tr->SetBranchStatus("elePt", 1);
    tr->SetBranchStatus("muPt", 1);
    tr->SetBranchStatus("jetPt", 1);
    tr->SetBranchStatus("nEle", 1);
    tr->SetBranchStatus("nMu", 1);
    tr->SetBranchStatus("nJet", 1);
    tr->SetBranchStatus("nBJet", 1);
    tr->SetBranchStatus("nPho", 1);
    tr->SetBranchStatus("TopHad_pt", 1);
    tr->SetBranchStatus("TopLep_pt", 1);
    tr->SetBranchStatus("TopHad_eta", 1);
    tr->SetBranchStatus("TopLep_eta", 1);
    tr->SetBranchStatus("TopHad_phi", 1);
    tr->SetBranchStatus("TopLep_phi", 1);
    //  tr->SetBranchStatus("TopHad_mass", 1);
    //  tr->SetBranchStatus("TopLep_mass", 1);
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

    // Set Branch Address
    tr->SetBranchAddress("elePt", &elePt);
    tr->SetBranchAddress("muPt", &muPt);
    tr->SetBranchAddress("jetPt", &jetPt);
    tr->SetBranchAddress("nEle", &nEle);
    tr->SetBranchAddress("nMu", &nMu);
    tr->SetBranchAddress("nJet", &nJet);
    tr->SetBranchAddress("nBJet", &nBJet);
    tr->SetBranchAddress("nPho", &nPho);
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

    // Set the loop

    double weight = 0;

    for (Int_t i = 0; i < numevents; i++)
    {
        tr->GetEntry(i);
        // define weight
        weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;

        TLorentzVector TopLep;
        TLorentzVector TopHad;
        TLorentzVector AntiTopLep;
        TLorentzVector AntiTopHad;

        // Electron Channel
        if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1) // selection

        {

            h1->Fill(elePt->at(0), weight);

            TopLep.SetPtEtaPhiM(TopLep_pt, TopHad_eta, TopLep_phi, TopLep_mass);
            h9ele->Fill(TopLep.Eta(), weight); // top quark
            AntiTopHad.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, TopHad_mass);
            h10ele->Fill(AntiTopHad.Pt(), weight); // anti top quark
            h11phi->Fill(TopLep_phi, weight);
            //  }
        }

        // Muon channel
        if (passPresel_Mu && nJet >= 4 && nBJet >= 1 && nPho == 1)

        {

            h2->Fill(muPt->at(0), weight);

            TopHad.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, TopHad_mass);
            h9mu->Fill(TopHad.Eta(), weight); // top quark
            AntiTopLep.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, TopLep_mass);
            h10mu->Fill(AntiTopLep.Pt(), weight); // anti top quark
            h12phi->Fill(TopHad_phi, weight);
            //}
        }
        if (nJet > 0)
        {
            //   for (Int_t j=0; j<nJet;j++)
            //{
            h3->Fill(jetPt->at(0), weight);
            //}
        }

        if (TopLep_charge > 0)
        {
            h9->Fill(TopLep_eta, weight);  // top quark
            h10->Fill(TopHad_eta, weight); // anti top quark
        }

        if (TopLep_charge < 0)
        {
            h9->Fill(TopHad_eta, weight);  // top quark
            h10->Fill(TopLep_eta, weight); // antitop quark
        }

        //   h4->Fill(TopHad_pt);
        //  h5->Fill(TopLep_pt);
        //   h6->Fill(TopHad_eta);
        //   h7->Fill(TopLep_eta);
        h8->Fill(TopLep_charge, weight);
    }

    // h9->Write();
    // h10->Write();

    h9ele->Write();
    h10ele->Write();
    h9mu->Write();
    h10mu->Write();

    int num_of_hists = sizeof(hists);
    TCanvas *c[num_of_hists];
    for (Int_t i = 0; i < sizeof(hists); i++)
    {
        c[i] = new TCanvas();
        c[i]->cd();
        hists[i]->Draw();
        c[i]->Update();
        TString filename = Form("/eos/user/s/ssnehshu/temp/%s/h%d.png", event, i + 1);
        c[i]->SaveAs(filename);
    }

    fl->Close();
    std::cout << "numevents" << numevents << std::endl;

    // std::cout << "nEle" << nEle << std::endl;
}
