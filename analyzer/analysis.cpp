#include <iostream>
#include <vector>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"

int main()
{
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

    std::vector<TString> input_files{ "nano_0.root" };
    TString tree_name = "Events";

    Analyzer ana(tree_name, input_files);
    ana.Analyze();

    return 0;
}