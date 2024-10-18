#include "Analyzer.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map): m_file_map(input_file_map), m_tree_name(tree_name) {}

void Analyzer::Analyze()
{
    for (auto const& [file_name, channel]: m_file_map)
    {
        TFile* file = TFile::Open(file_name);
        TTree* tree = static_cast<TTree*>(file->Get(m_tree_name));

        TString channel_name = channel == Channel::SL ? "Single lepton" : "Double Lepton";
        m_index = channel == Channel::SL ? GenTruthIdxMapSL : GenTruthIdxMapDL;
        std::cout << "File: " << file_name << ", channel: " << channel_name << "\n";

        Event event(tree, channel);
        AnalyzeEvent(event, tree);
        file->Close();
    }
}

void Analyzer::AnalyzeEvent(Event const& event, TTree* tree)
{
    for (long long i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        size_t met_idx = m_index.at("met");
        std::cout << event.gen_truth.pt[met_idx] << " " << event.gen_truth.mass[met_idx] << "\n";

        std::cout << event.recojet.nRecoJet << "\n";

        // break;
    }
}