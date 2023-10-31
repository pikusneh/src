#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"

void response()
{
    TFile *file1 = new TFile("output_GenLevel.root", "READ");
    TFile *file2 = new TFile("output_TTGamma.root", "READ");

    TTree *tree_gen = (TTree *)file1->Get("tree_gen");
    TTree *tree_reco = (TTree *)file2->Get("tree_reco");

    Float_t rapidity_diff_reco, event_weight_reco;
    Float_t rapidity_diff_gen, event_weight_gen;

    tree_reco->SetBranchAddress("rapidity_diff_reco", &rapidity_diff_reco);
    tree_reco->SetBranchAddress("event_weight_reco", &event_weight_reco);
    tree_gen->SetBranchAddress("rapidity_diff_gen", &rapidity_diff_gen);
    tree_gen->SetBranchAddress("event_weight_gen", &event_weight_gen);

    Int_t nbins = 20;
    Double_t xmin = -3.0, xmax = 3.0;

    TH2F *responseMatrix = new TH2F("responseMatrix", "Response Matrix", nbins, xmin, xmax, nbins, xmin, xmax);
    TH1F *histReco = new TH1F("histReco", "Reco level", nbins, xmin, xmax);
    TH1F *histGen = new TH1F("histGen", "Gen level", nbins, xmin, xmax);

    histReco->SetOption("hist");
    histGen->SetOption("hist");

    Int_t nentriesReco = tree_reco->GetEntries();
    for (Int_t i = 0; i < nentriesReco; i++)
    {
        tree_reco->GetEntry(i);
        histReco->Fill(rapidity_diff_reco, event_weight_reco);
    }

    Int_t nentriesGen = tree_gen->GetEntries();
    for (Int_t i = 0; i < nentriesGen; i++)
    {
        tree_gen->GetEntry(i);
        histGen->Fill(rapidity_diff_gen, event_weight_gen);
        // responseMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
    }
    for (Int_t i = 0; i < nentriesReco; i++)
    {
        tree_reco->GetEntry(i);
        tree_gen->GetEntry(i);
        responseMatrix->Fill(rapidity_diff_gen, rapidity_diff_reco, event_weight_reco * event_weight_gen);
    }

    TCanvas *canvas = new TCanvas("canvas", "Response Matrix Canvas", 800, 600);
    responseMatrix->Draw("COLZTEXT");
    responseMatrix->GetXaxis()->SetTitle("Gen Level");
    responseMatrix->GetYaxis()->SetTitle("Reco Level");
    responseMatrix->SetTitle("Response Matrix");

    canvas->SaveAs("responseMatrix.png");

    TFile *outputFile = new TFile("outputFile.root", "RECREATE");
    responseMatrix->Write();
    histReco->Write();
    histGen->Write();
    canvas->Write();
    file1->Close();
    file2->Close();
    outputFile->Close();
    // delete file1;
    // delete file2;
    // delete outputFile;
    // delete responseMatrix;
    // delete histReco;
    // delete histGen;
    // delete canvas;
}
