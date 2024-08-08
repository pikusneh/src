#include "GenRecoAnalyzer.h"
#include "TCanvas.h"
#include <iostream>


// Constructor
GenRecoAnalyzer::GenRecoAnalyzer(const char *inputFileName, const char *outputFileName)
{
    chain = new TChain("AnalysisTree");
    chain->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
    chain->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
    chain->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");

    outputFile = new TFile(outputFileName, "RECREATE");
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
    double weight = 0;

    InitializeHistograms();
    SetupBranches();
}

// Destructor
GenRecoAnalyzer::~GenRecoAnalyzer()
{
    delete chain;
    outputFile->Close();
    delete outputFile;
}

// Initialize histograms
void GenRecoAnalyzer::InitializeHistograms()
{
    Rapidity_diff_Gen = new TH1F("Rapidity_diff_Gen", "Rapidity_Diff_Gen", 20, -3, 3);
    Rapidity_diff_Reco = new TH1F("Rapidity_diff_Reco", "Rapidity_Diff_Reco", 20, -3, 3);
    correlMatrix = new TH2F("correlMatrix", "correl Matrix", 30, -3, 3, 30, -3, 3);
    correlMatrix_reco = new TH2F("correlMatrix_reco", "correl Matrix_reco", 20, -3, 3, 20, -3, 3);
}

void GenRecoAnalyzer::SetupBranches()
{
    // Example branch setup
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("elePt", 1);
    chain->SetBranchStatus("muPt", 1);
    chain->SetBranchStatus("jetPt", 1);
    chain->SetBranchStatus("phoEt", 1);
    chain->SetBranchStatus("phoEta", 1);
    chain->SetBranchStatus("phoPhi", 1);
    chain->SetBranchStatus("nEle", 1);
    chain->SetBranchStatus("nMu", 1);
    chain->SetBranchStatus("nGenJet", 1);
    chain->SetBranchStatus("nGenPart", 1);
    chain->SetBranchStatus("genPDGID", 1);
    chain->SetBranchStatus("genMomIdx", 1);
    chain->SetBranchStatus("genMass", 1);
    chain->SetBranchStatus("genPhi", 1);
    chain->SetBranchStatus("genEta", 1);
    chain->SetBranchStatus("genPt", 1);
    chain->SetBranchStatus("evtWeight", 1);
    chain->SetBranchStatus("nJet", 1);
    chain->SetBranchStatus("nBJet", 1);
    chain->SetBranchStatus("nPho", 1);
    chain->SetBranchStatus("TopHad_pt", 1);
    chain->SetBranchStatus("TopLep_pt", 1);
    chain->SetBranchStatus("TopHad_eta", 1);
    chain->SetBranchStatus("TopLep_eta", 1);
    chain->SetBranchStatus("TopHad_phi", 1);
    chain->SetBranchStatus("TopLep_phi", 1);
    chain->SetBranchStatus("M_bjj", 1);
    chain->SetBranchStatus("Mt_blgammaMET", 1);
    chain->SetBranchStatus("TopLep_charge", 1);
    chain->SetBranchStatus("passPresel_Ele", 1);
    chain->SetBranchStatus("passPresel_Mu", 1);
    chain->SetBranchStatus("lumis", 1);
    chain->SetBranchStatus("btagWeight_1a", 1);
    chain->SetBranchStatus("prefireSF", 1);
    chain->SetBranchStatus("muEffWeight", 1);
    chain->SetBranchStatus("eleEffWeight", 1);
    chain->SetBranchStatus("PUweight", 1);

    // Set branch addresses
    chain->SetBranchAddress("elePt", &elePt);
    chain->SetBranchAddress("muPt", &muPt);
    chain->SetBranchAddress("jetPt", &jetPt);
    chain->SetBranchAddress("phoEt", &phoEt);
    chain->SetBranchAddress("phoEta", &phoEta);
    chain->SetBranchAddress("phoPhi", &phoPhi);
    chain->SetBranchAddress("nEle", &nEle);
    chain->SetBranchAddress("nMu", &nMu);
    chain->SetBranchAddress("nGenJet", &nGenJet);
    chain->SetBranchAddress("nGenPart", &nGenPart);
    chain->SetBranchAddress("genPDGID", &genPDGID);
    chain->SetBranchAddress("genMomIdx", &genMomIdx);
    chain->SetBranchAddress("genMass", &genMass);
    chain->SetBranchAddress("genPhi", &genPhi);
    chain->SetBranchAddress("genEta", &genEta);
    chain->SetBranchAddress("genPt", &genPt);
    chain->SetBranchAddress("evtWeight", &evtWeight);
    chain->SetBranchAddress("nJet", &nJet);
    chain->SetBranchAddress("nBJet", &nBJet);
    chain->SetBranchAddress("nPho", &nPho);
    chain->SetBranchAddress("Mt_blgammaMET", &Mt_blgammaMET); // toplepmass
    chain->SetBranchAddress("M_bjj", &M_bjj);                 // tophadmass
    chain->SetBranchAddress("TopHad_pt", &TopHad_pt);
    chain->SetBranchAddress("TopLep_pt", &TopLep_pt);
    chain->SetBranchAddress("TopHad_eta", &TopHad_eta);
    chain->SetBranchAddress("TopLep_eta", &TopLep_eta);
    chain->SetBranchAddress("TopHad_phi", &TopHad_phi);
    chain->SetBranchAddress("TopLep_phi", &TopLep_phi);
    chain->SetBranchAddress("TopLep_charge", &TopLep_charge);
    chain->SetBranchAddress("passPresel_Ele", &passPresel_Ele);
    chain->SetBranchAddress("passPresel_Mu", &passPresel_Mu);
    chain->SetBranchAddress("lumis", &lumis);
    chain->SetBranchAddress("btagWeight_1a", &btagWeight_1a);
    chain->SetBranchAddress("evtWeight", &evtWeight);
    chain->SetBranchAddress("prefireSF", &prefireSF);
    chain->SetBranchAddress("muEffWeight", &muEffWeight);
    chain->SetBranchAddress("eleEffWeight", &eleEffWeight);
    chain->SetBranchAddress("PUweight", &PUweight);
}

// Main analysis loop
void GenRecoAnalyzer::RunAnalysis()
{
    int numEvents = chain->GetEntries();
    for (int i = 0; i < numEvents; ++i)
    {
        ProcessEvent(i);
    }

    SaveResults();
}

// Process a single event
void GenRecoAnalyzer::ProcessEvent(int eventNumber)
{
    // Load the data for the current event
    chain->GetEntry(eventNumber);
    weight_gen = evtWeight;
    weight_reco = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;
    AnalyzeGenParticles();
    AnalyzeRecoParticles();

    // // Calculate rapidity differences between particles
    // CalculateRapidityDifferences();
    FillHistograms();
}

void GenRecoAnalyzer::AnalyzeGenParticles()
{
    int topIndex = -1, antiTopIndex = -1, photonIndex = -1, electronIndex = -1;
    bool passPresel_Ele_gen = false;

    // Counters for photons and b-jets
    int totalPhotons = 0, totalbjet = 0;
    // Loop over all generator level particles
    for (int i = 0; i < nGenPart; ++i)
    {
        int pdgId = genPDGID->at(i);

        // Identify top, anti-top, electron, and photon indices
        if (pdgId == 6)
        {
            topIndex = i;
        }
        else if (pdgId == -6)
        {
            antiTopIndex = i;
        }
        else if (pdgId == 11)
        {
            electronIndex = i;
            // Calculate and check electron properties
            TLorentzVector ele_gen;
            ele_gen.SetPtEtaPhiM(genPt->at(i), genEta->at(i), genPhi->at(i), genMass->at(i));
            if (ele_gen.Pt() >= 35.0 && fabs(ele_gen.Rapidity()) <= 2.4)
            {
                passPresel_Ele_gen = true;
            }
        }
        else if (pdgId == 22)
        {
            photonIndex = i;
            totalPhotons++;
        }

        // Count b-jets
        if ((pdgId == 5 || pdgId == -5) && genPt->at(i) >= 30.0 && fabs(genEta->at(i)) <= 2.4)
        {
            totalbjet++;
        }
    }

    // Apply selection criteria
    if (passPresel_Ele_gen && nGenJet >= 4 && totalbjet >= 1 && totalPhotons == 1)
    {
        TLorentzVector top_gen, antiTop_gen, photon_gen;
        photon_gen.SetPtEtaPhiM(genPt->at(photonIndex), genEta->at(photonIndex), genPhi->at(photonIndex), genMass->at(photonIndex));

        if (fabs(photon_gen.Rapidity()) >= 1)
        {
            top_gen.SetPtEtaPhiM(genPt->at(topIndex), genEta->at(topIndex), genPhi->at(topIndex), genMass->at(topIndex));
            antiTop_gen.SetPtEtaPhiM(genPt->at(antiTopIndex), genEta->at(antiTopIndex), genPhi->at(antiTopIndex), genMass->at(antiTopIndex));

            // Calculate rapidity differences
            float rapidity_diff_gen = fabs(top_gen.Rapidity()) - fabs(antiTop_gen.Rapidity());

            // Fill the histogram
            Rapidity_diff_Gen->Fill(rapidity_diff_gen, weight_gen);
            std::cout << "Gen Level Rapidity Difference: " << rapidity_diff_gen << std::endl;
        }
    }
}

void GenRecoAnalyzer::AnalyzeRecoParticles()
{
    // Check for preselection criteria: number of jets, b-jets, and photons
    if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1)
    {
        // Initialize TLorentzVector for Photon, Top, and AntiTop
        TLorentzVector Photon, Top, AntiTop;

        // Set Photon properties
        Photon.SetPtEtaPhiM(phoEt->at(0), phoEta->at(0), phoPhi->at(0), 0.0);
        float PhoRapidity = Photon.Rapidity();
        float YPho = fabs(PhoRapidity);

        // Check Photon rapidity
        if (YPho >= 1)
        {
            // Assign Top and AntiTop based on charge
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

            // Calculate rapidity and rapidity differences
            float Rapidity_T_ele = Top.Rapidity();
            float Rapidity_t_ele = AntiTop.Rapidity();
            float YT_ele = fabs(Rapidity_T_ele); // absolute value of top rapidity
            float Yt_ele = fabs(Rapidity_t_ele); // absolute value of antitop rapidity
            float rapidity_diff_ele = YT_ele - Yt_ele;

            // Fill histograms and matrices
            Rapidity_diff_Reco->Fill(rapidity_diff_ele, weight_reco);
            correlMatrix_reco->Fill(rapidity_diff_ele, rapidity_diff_ele);
            std::cout << "Reco Level Rapidity Difference: " << rapidity_diff_ele << std::endl;

            
        }
    }
}

// void GenRecoAnalyzer::CalculateRapidityDifferences()
// {
//     // Implementation of rapidity difference calculation
// }

void GenRecoAnalyzer::FillHistograms()
{
    // Fill 1D histograms with the calculated rapidity differences
    if (Rapidity_diff_Gen && Rapidity_diff_Reco)
    {
        Rapidity_diff_Gen->Fill(rapidity_diff_gen, weight_gen);
        Rapidity_diff_Reco->Fill(rapidity_diff_reco, weight_reco);
    }
    else
    {
        std::cerr << "1D Histograms not initialized properly!" << std::endl;
    }

    // Fill the 2D histogram with correlation data, if applicable
    if (correlMatrix)
    {
        correlMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, weight * evtWeight);
    }
    else
    {
        std::cerr << "2D Histogram (correlMatrix) not initialized properly!" << std::endl;
    }
}

// Helper method: Calculate correlations
void GenRecoAnalyzer::CalculateCorrelations()
{
    // Implementation of correlation calculation
}

// Save the results to the output file
void GenRecoAnalyzer::SaveResults()
{
    outputFile->cd();

     // Save histograms as canvas objects
    TCanvas* canvas1 = new TCanvas("canvas1", "Rapidity Difference Gen", 600, 400);
    Rapidity_diff_Gen->Draw();
    canvas1->Write();

    TCanvas* canvas2 = new TCanvas("canvas2", "Rapidity Difference Reco", 600, 400);
    Rapidity_diff_Reco->Draw();
    canvas2->Write();

    TCanvas* canvas3 = new TCanvas("canvas3", "Correlation Matrix", 600, 400);
    correlMatrix->Draw("COLZ"); // Draw as a color-filled plot
    canvas3->Write();

    TCanvas* canvas4 = new TCanvas("canvas4", "Correlation Matrix Reco", 600, 400);
    correlMatrix_reco->Draw("COLZ");
    canvas4->Write();

    Rapidity_diff_Gen->Write();
    Rapidity_diff_Reco->Write();
    correlMatrix->Write();
    correlMatrix_reco->Write();
    // Write other results as needed
}

// Additional helper method implementations...
