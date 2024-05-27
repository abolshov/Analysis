#include <iostream>
#include <vector>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"

int main()
{
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

    std::vector<TString> input_files{ "GluGluToRadionToHHTo2B2WToLNu2J_M-450.root", 
                                      "GluGluToRadionToHHTo2B2WToLNu2J_M-500.root" };
    TString tree_name = "Events";

    Analyzer ana(tree_name, input_files);
    ana.Analyze();

    return 0;
}