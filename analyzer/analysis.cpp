#include <iostream>
#include <vector>
#include <map>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"
#include "Constants.hpp"

int main()
{
    // gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

    TString tree_name = "Events";
    TString pdf_file_name = "pdf_sl.root";
    // TString pdf_file_name = "pdf_dl.root";

    std::map<TString, Channel> input_file_map = { { "nano_sl_M800.root", Channel::SL },
                                                  { "nano_dl_M800.root", Channel::DL } };

    Analyzer ana(tree_name, input_file_map, pdf_file_name);
    ana.ProcessFile("nano_sl_M800.root", Channel::SL);
    // ana.ProcessFile("nano_dl_M800.root", Channel::DL);

    return 0;
}