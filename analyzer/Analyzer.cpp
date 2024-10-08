#include "Analyzer.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::vector<TString> const& input_files)
: m_chain(std::make_unique<TChain>(tree_name)),
  m_input_files(input_files)
{
    for (auto const& input_file: m_input_files)
    {
        m_chain->Add(input_file);
    }
    m_event = std::make_unique<Event>(m_chain.get()); 
}

void Analyzer::Analyze()
{
    std::cout << "Files in TChain:\n";
    for (auto const& input_file: m_input_files)
    {
        std::cout << "\t" << input_file << "\n";
    }

    auto nEntries = m_chain->GetEntries();
    std::cout << "chain size: " << nEntries << "\n";

    for (size_t i = 0; i < 1; ++i)
    {
        m_chain->GetEntry(i);

        std::cout << "gen=" << m_event->genjet.nGenJet << ", reco=" << m_event->recojet.nRecoJet << "\n";
    }
}