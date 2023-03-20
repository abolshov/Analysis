#include <iostream>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include <math.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TSpline.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TStyle.h"
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <functional>
#include "TLegend.h"
#include <numeric>
#include "TRandom3.h"
#include <cmath>
#include <regex>
#include <fstream>
#include <cstdio>

int main()
{
    TFile *myFile = TFile::Open("out_radion_2L_2016_m400.root");
    // TFile *myFile = TFile::Open("out_ttbar_2L_2016.root");
    TDirectory *dir = (TDirectory*)myFile;
    TTree *myTree = (TTree*)dir->Get("Double_Tree");

    std::ofstream file;
    file.open("data_prep.cpp");
    file << std::endl;
    file << "int main() {\n" 
         << "\tTFile *myFile = TFile::Open(\"out_radion_2L_2016_m400.root\");\n"
         << "\tTDirectory *dir = (TDirectory*)myFile;\n"
         << "\tTTree *myTree = (TTree*)dir->Get(\"Double_Tree\");\n"
         << std::endl;

    int nEntries = myTree->GetEntries();

    std::cout << "nEntries = " << nEntries << std::endl;

    auto branches = myTree->GetListOfBranches();
    int nBranches = myTree->GetNbranches();

    std::vector<std::string> col_jet;
    std::vector<std::string> col_met;
    std::vector<std::string> col_lep;

    std::vector<std::string> branch_names;

    std::regex pat_jet("_jet[0-9]_[pE]");
    std::regex pat_met("met_[pE]");
    std::regex pat_lep("lep[0-9]_[pE]");

    // char formatted[1000];
    // std::string line = "function(%s);";
    // std::sprintf(formatted, line.c_str(), "variable");
    // std::cout << "formatted = " << formatted << std::endl;

    std::cout << "nBranches = " << nBranches << std::endl;
    // myTree->Print();
    for (size_t i = 0; i < nBranches; ++i)
    {
        std::string bName = branches->At(i)->GetName();
        size_t found = bName.find("jet");
        std::cout << branches->At(i)->GetName() << std::endl;
        if (found != std::string::npos)
        {
            std::cout << bName << /* ": " << myTree->GetBranch(bName.c_str())->GetEntries() << */ std::endl;
        }
        std::string decl("\tfloat %s;\n");
        std::string set_addr("\tmyTree->SetBranchAddress(\"%s\", &%s);\n");
        char buf[1000];
        if (std::regex_search(bName, pat_lep))
        {
            col_lep.push_back(bName);
            std::cout << "Branch " << bName << " pushed back to col_lep" << std::endl;
            std::sprintf(buf, decl.c_str(), bName.c_str());
            file << buf;
            std::sprintf(buf, set_addr.c_str(), bName.c_str(), bName.c_str());
            file << buf;
            file << std::endl;
            branch_names.push_back(bName);
        }
        else if (std::regex_search(bName, pat_jet))
        {
            col_jet.push_back(bName);
            std::cout << "Branch " << bName << " pushed back to col_jet" << std::endl;
            std::sprintf(buf, decl.c_str(), bName.c_str());
            file << buf;
            std::sprintf(buf, set_addr.c_str(), bName.c_str(), bName.c_str());
            file << buf;
            file << std::endl;
            branch_names.push_back(bName);
        }
        else if (std::regex_search(bName, pat_met))
        {
            col_jet.push_back(bName);
            std::cout << "Branch " << bName << " pushed back to col_met" << std::endl;
            std::sprintf(buf, decl.c_str(), bName.c_str());
            file << buf;
            std::sprintf(buf, set_addr.c_str(), bName.c_str(), bName.c_str());
            file << buf;
            file << std::endl;
            branch_names.push_back(bName);
        }
    }

    std::cout << "# branches = " << branch_names.size() << std::endl;

    for (auto const& name: branch_names)
    {
        std::cout << " << " << name << " << \",\"" << std::endl;
    }
    file << "\treturn 0;\n" << "}\n";
    file.close();

    return 0;
}