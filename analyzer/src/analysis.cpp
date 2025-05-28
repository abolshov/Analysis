#include <iostream>
#include <vector>
#include <map>

#include "TString.h"
#include "TROOT.h"

#include "Analyzer.hpp"
#include "Constants.hpp"
#include "Definitions.hpp"
#include "Sample.hpp"

int main()
{
    #if (!defined(DEV) && !defined(DEBUG))
        std::cout << "Release mode: ignoring all ROOT warnings\n";
        gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");
    #endif
    std::cout << std::setprecision(3);

    TH1::AddDirectory(false);
    TH2::AddDirectory(false);

    std::map<Channel, TString> pdf_file_map = { { Channel::SL, "pdf_sl.root" },
                                                { Channel::DL, "pdf_dl.root" } };

    Dataset_t dataset = { {"nano_sl_M1000.root", "Events", "GGF", Channel::SL, 1000} };

    Analyzer ana(pdf_file_map);
    ana.ProcessDataset(dataset);

    return 0;
}