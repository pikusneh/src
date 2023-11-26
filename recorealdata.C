
#include <TH1F.h>
#include <TFile.h>

void recorealdata()
{

    TFile *fileEle = new TFile("Ele.root", "READ");

    TFile *fileCombined = new TFile("combined_histogram.root", "READ");

    TH1F *histEle = (TH1F *)fileEle->Get("h11");
    histEle->SetOption("hist");

    TH1F *histCombined = (TH1F *)fileCombined->Get("h25");

    // Subtract the number of events in h25 of combinedhistogram.root from h11 of Ele.root
    histEle->Add(histCombined, -1);

    // Save the resulting histogram to a new file
    TFile *outputFile = new TFile("subtracted_histogram.root", "RECREATE");
    histEle->Write();
    outputFile->Close();

    // Close the files
    fileEle->Close();
    fileCombined->Close();
}
