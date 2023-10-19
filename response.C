#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include <TH2.h>

void response()
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
       

        // Loop over Reco level histogram
        for (Int_t j = 1; j <= nbins; j++)
        {
            // Get the number of events in the Reco level bin and fill histreco
            Float_t recoLevelEvents = histreco->GetXaxis()->GetBinCenter(j);
            

            Float_t genContent = histgen->GetBinContent(i);
            Float_t recoContent = histreco->GetBinContent(j);

            // Calculate the response matrix element
        
            histresponse->SetBinContent(j, i, genContent / recoContent);
        
        }
    }

    // Save histresponse to a ROOT file
    TFile *outputFile = new TFile("response.root", "RECREATE");
    histresponse->Write();
    outputFile->Close();
}