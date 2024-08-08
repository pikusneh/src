

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// #include "RooUnfoldSvd.h"
// #include "RooUnfoldTUnfold.h"
// #include "RooUnfoldIds.h"
#endif

const Double_t cutdummy = -99999.0;
#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include <TH2.h>

#include <iostream>
using std::cout;
using std::endl;

void CalculateResponseMatrix()
{
    // Open the ROOT files containing the histograms
    TFile *file1 = new TFile("output_GenLevel.root", "READ");
    TFile *file2 = new TFile("output_TTGamma.root", "READ");

    // Get the histograms from the files
    TH1F *hist1 = (TH1F *)file1->Get("h7");
    TH1F *hist2 = (TH1F *)file2->Get("h25");

    TH1F *histgen = new TH1F("histgen", "Gen Level", 4, -3, 3);
    TH1F *histreco = new TH1F("histreco", "Reco Level", 4, -3, 3);
    TH2F *histresponse = new TH2F("histresponse", "histresponse", 4, -3, 3, 4, -3, 3);

    // Fill histgen and histreco with the content of hist1 and hist2, respectively
    histgen->Add(hist1);
    histreco->Add(hist2);

    Int_t nbins = 4;
    // Loop over Gen level histogram
    for (Int_t i = 1; i <= nbins; i++)
    {
        // Get the number of events in the Gen level bin and fill histgen
        Float_t genLevelEvents = histgen->GetXaxis()->GetBinCenter(i);
        // histgen->Fill(genLevelEvents);

        // Loop over Reco level histogram
        for (Int_t j = 1; j <= nbins; j++)
        {
            // Get the number of events in the Reco level bin and fill histreco
            Float_t recoLevelEvents = histreco->GetXaxis()->GetBinCenter(j);
            //  histreco->Fill(recoLevelEvents);

            Float_t genContent = histgen->GetBinContent(i);
            Float_t recoContent = histreco->GetBinContent(j);

            // Calculate the response matrix element
           // Float_t response = histgen->GetBinContent(i) / histreco->GetBinContent(j);
            //  Float_t response = histreco->GetBinContent(j) / histgen->GetBinContent(i);
             histresponse->SetBinContent(j, i, genContent / recoContent);
            // Fill histresponse with the response matrix element
           // histresponse->SetBinContent(i, j, response);
        }
    }

    // // Create a RooUnfoldResponse object
    // RooUnfoldResponse response(histgen, histreco, histresponse);

    // cout << "==================================== TEST =====================================" << endl;
    // // TFile *subtracted_file = new TFile("subtracted_histogram.root", "READ");
    // TFile *file3 = new TFile("subtracted_histogram.root", "READ");
    // TH1F *h11 = (TH1F *)file3->Get("h11");
    // TFile *file4 = new TFile("output_TTGamma.root", "READ");
    // TH1F *h25 = (TH1F *)file4->Get("h25");
    // cout << "==================================== UNFOLD ===================================" << endl;
    // RooUnfoldBayes unfold(&response, h25, 4); // h11-> test input

    // TH1D *hUnfold = (TH1D *)unfold.Hunfold();
    // hUnfold->SetOption("hist");

    // TCanvas *c1 = new TCanvas("canvas", "canvas");

    // unfold.PrintTable(cout, hist2); //(cout,test truth)
    // hUnfold->Draw();
    // h25->SetLineColor(8);
    // h25->Draw("SAME");

    // Print the content of histgen and histreco
    // histgen->Print();
    // histreco->Print();

    // Save histresponse to a ROOT file
    TFile *outputFile = new TFile("responseMatrix.root", "RECREATE");
    histresponse->Write();
    outputFile->Close();

    // Clean up
    // delete histgen;
    // delete histreco;
    // delete histresponse;
}
