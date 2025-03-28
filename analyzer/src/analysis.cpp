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
    std::cout << std::setprecision(3);
    TString tree_name = "Events";
    
    std::map<TString, Channel> input_file_map = { { "nano_sl_M800.root", Channel::SL },
                                                  { "nano_dl_M800.root", Channel::DL } };

    std::map<Channel, TString> pdf_file_map = { { Channel::SL, "pdf_sl.root" },
                                                { Channel::DL, "pdf_dl.root" } };

    Analyzer ana(tree_name, input_file_map, pdf_file_map);
    ana.ProcessFile("nano_sl_M800.root", Channel::SL);
    // ana.ProcessFile("nano_dl_M800.root", Channel::DL);

    return 0;
}