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
    TString pdf_file_name = "pdf.root";

    std::map<TString, Channel> input_file_map = { { "nano_0.root", Channel::SL } };

    Mode mode = Mode::Validation;

    Analyzer ana(tree_name, input_file_map, pdf_file_name, mode);
    // ana.ProcessFile("nano_0.root", Channel::SL);

    return 0;
}