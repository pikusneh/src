#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>

void gauss() {
    // Create a ROOT file to store the histogram
    TFile* file = new TFile("gaussian_histogram.root", "RECREATE");

    // Create a random number generator
    TRandom3 random;

    // Define Gaussian parameters
    double mean = 0.0;
    double sigma = 1.0;
    int numBins = 100;
    double xmin = -5.0;
    double xmax = 5.0;

    // Create a histogram
    TH1F* hist = new TH1F("gaussian_hist", "Gaussian Histogram", numBins, xmin, xmax);

    // Fill the histogram with random Gaussian-distributed values
    for (int i = 0; i < 10000; i++) {
        double value = random.Gaus(mean, sigma);
        hist->Fill(value);
    }

    // Write the histogram to the ROOT file
    hist->Write();

    // Close the ROOT file
    file->Close();

    std::cout << "Histogram saved to gaussian_histogram.root" << std::endl;

}
