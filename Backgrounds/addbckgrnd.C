#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1.h"

float addbckgrnd
() {
  // List of root file names
  std::vector<std::string> fileNames = {
      "output_diboson.root",
      "output_qcd_ele.root",
      "output_qcd_mu.root",
      "output_singletop.root",
      "output_tjets.root",
      "output_wgamma.root",
      "output_wjets.root",
      "output_zgamma.root",
      "output_zjets.root",
      "output_TTbar.root"
  };

  // Create a vector to store histograms
  std::vector<TH1*> histograms;

  // Open the root files and get the histograms
  for (const auto& fileName : fileNames) {
    TFile* rootFile = new TFile(fileName.c_str(), "READ");
    TH1* histogram = (TH1*)rootFile->Get("h25");
    histograms.push_back(histogram);
  }

  // Check if the histograms were loaded successfully
  if (histograms.empty()) {
    std::cerr << "Error: No histograms were loaded." << std::endl;
    return 1;
  }

  // Create a combined histogram
  TH1* combinedHistogram = (TH1*)histograms[0]->Clone();
  combinedHistogram->Reset();

  // Add the histograms
  for (size_t i = 0; i < histograms.size(); ++i) {
    combinedHistogram->Add(histograms[i]);
  }

  // Do something with the combined histogram
  combinedHistogram->Draw();

  // Clean up memory
  for (auto& histogram : histograms) {
    delete histogram;
  }
  // Create a new ROOT file to save the combined histogram
TFile* outputFile = new TFile("combined_histogram.root", "RECREATE");

// Write the combined histogram to the output file
combinedHistogram->Write();

// Close the output file
outputFile->Close();


  return 0;
}
