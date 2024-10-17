#include <iostream>
#include <vector>
#include <cassert>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"
#include "Constants.hpp"

int main()
{
    static_assert(static_cast<int>(Obj::count) == 6);

    gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

    std::vector<TString> input_files{ "nano_0.root" };
    TString tree_name = "Events";

    Analyzer ana(tree_name, input_files);
    ana.Analyze();

    return 0;
}