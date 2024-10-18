#include <iostream>
#include <vector>
#include <map>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"
#include "Constants.hpp"

int main()
{
    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

    TString tree_name = "Events";

    std::map<TString, Channel> input_file_map = { { "nano_0.root", Channel::SL } };

    Analyzer ana(tree_name, input_file_map);
    ana.Analyze();

    return 0;
}