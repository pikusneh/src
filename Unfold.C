
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
    TH1F *hist1 = (TH1F *)file1->Get("h7");
    TFile *file2 = new TFile("output_TTGamma.root", "READ");
    TH1F *hist2 = (TH1F *)file2->Get("h25");

    // Create the response matrix
    RooUnfoldResponse response(hist1->GetNbinsX(), hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());

    // Fill the response matrix
    for (Int_t i = 1; i <= hist1->GetNbinsX(); i++)
    {
        Double_t xt = hist1->GetBinCenter(i);
        if (xt >= hist2->GetXaxis()->GetXmin() && xt <= hist2->GetXaxis()->GetXmax())
        {
            Int_t bin = hist2->FindBin(xt);
            response.Fill(hist2->GetBinContent(bin), xt);
        }
        else
        {
            response.Miss(xt);
        }
    }
    
    // Get the 2-dimensional histogram from the response matrix
    TH2 *hist2D = response.Hresponse();

    // Optional: Create a projection of the 2D histogram to obtain a 1D histogram
    TH1D *hist = hist2D->ProjectionX();

    // Optional: Set the title and axis labels for the histogram
    hist->SetTitle("Response Histogram");
    hist->GetXaxis()->SetTitle("Measured Value");
    hist->GetYaxis()->SetTitle("True Value");

    // Print the response matrix
    response.Print();

    // unable to invert the response matrix.

    // try inverting the response matrix and then unfold. it should come similar to gen level.

    cout << "==================================== TEST =====================================" << endl;
    // TFile *subtracted_file = new TFile("subtracted_histogram.root", "READ");
    TFile *file = new TFile("subtracted_histogram.root", "READ");
    TH1F *h11 = (TH1F *)file->Get("h11");

    cout << "==================================== UNFOLD ===================================" << endl;
    RooUnfoldBayes unfold(&response, h11, 4);

    TH1D *hUnfold = (TH1D *)unfold.Hunfold();
    hUnfold->SetOption("hist");

    TCanvas *c1 = new TCanvas("canvas", "canvas");

    unfold.PrintTable(cout, hist2);
    hUnfold->Draw();
    h11->SetLineColor(8);
    h11->Draw("SAME");

    // hTrue->Draw("SAME");

    // Save the histogram and the response matrix to a ROOT file
    TFile *outputFile = new TFile("unfolding_output.root", "RECREATE");
    hUnfold->Write();
    h11->Write();
    hist2D->Write();
    outputFile->Close();

    c1->SaveAs("RooUnfoldExample.png");
}

#ifndef __CINT__
int main()
{
    Unfold();
    return 0;
} // Main program when run stand-alone
#endif
