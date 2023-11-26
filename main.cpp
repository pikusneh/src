#include "GenRecoAnalyzer.h"

int main() {
    GenRecoAnalyzer analyzer("inputfile.txt", "outputfile.root");
    analyzer.RunAnalysis();
    return 0;
}
