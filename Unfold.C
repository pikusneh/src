
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

void Unfold()
{

    TFile *file1 = new TFile("output_GenLevel.root", "READ");
    TH1F *hist_gen = (TH1F *)file1->Get("Rapidity_Diff_Gen");
    TFile *file2 = new TFile("output_TTGamma.root", "READ");
    TH1F *hist_reco = (TH1F *)file2->Get("h25");

    // Create the response matrix
    RooUnfoldResponse response(hist_gen->GetNbinsX(), hist_gen->GetXaxis()->GetXmin(), hist_gen->GetXaxis()->GetXmax());

    // Fill the response matrix
    for (Int_t i = 1; i <= hist_gen->GetNbinsX(); i++)
    {
        Double_t xt = hist_gen->GetBinCenter(i);
        if (xt >= hist_reco->GetXaxis()->GetXmin() && xt <= hist_reco->GetXaxis()->GetXmax())
        {
            Int_t bin = hist_reco->FindBin(xt);
            response.Fill(hist_reco->GetBinContent(bin), xt);
        }
        else
        {
            response.Miss(xt);
        }
    }

    // Get the 2-dimensional histogram from the response matrix
    TH2 *hist_recoD = response.Hresponse();

    // Optional: Create a projection of the 2D histogram to obtain a 1D histogram
    TH1D *hist = hist_recoD->ProjectionX();

    // Optional: Set the title and axis labels for the histogram
    hist->SetTitle("Response Histogram");
    hist->GetXaxis()->SetTitle("Measured Value");
    hist->GetYaxis()->SetTitle("True Value");

    // Print the response matrix
    response.Print();
    // Save response matrix to ROOT file
   
    // unable to invert the response matrix.

    // try inverting the response matrix and then unfold. it should come similar to gen level.

    // cout << "==================================== TEST =====================================" << endl;
    // // TFile *subtracted_file = new TFile("subtracted_histogram.root", "READ");
    // TFile *file = new TFile("subtracted_histogram.root", "READ");
    // TH1F *h11 = (TH1F *)file->Get("h11");

    // cout << "==================================== UNFOLD ===================================" << endl;
    // RooUnfoldBayes unfold(&response, h11, 4);

    // TH1D *hUnfold = (TH1D *)unfold.Hunfold();
    // hUnfold->SetOption("hist");

    // TCanvas *c1 = new TCanvas("canvas", "canvas");

    // unfold.PrintTable(cout, hist_reco);
    // hUnfold->Draw();
    // h11->SetLineColor(8);
    // h11->Draw("SAME");

    // hTrue->Draw("SAME");

    // Save the histogram and the response matrix to a ROOT file
    TFile *outputFile = new TFile("unfolding_output.root", "RECREATE");
    // hUnfold->Write();
    // h11->Write();
    hist_recoD->Write();
    outputFile->Close();

    // c1->SaveAs("RooUnfoldExample.png");
}

#ifndef __CINT__
int main()
{
    Unfold();
    return 0;
} // Main program when run stand-alone
#endif
