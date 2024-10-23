#include "Analyzer.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)
,   m_pdf_file_name(pdf_file_name)     
{}

HistVec1d_t Analyzer::PDFs1dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec1d_t res;
    auto file = std::make_unique<TFile>(file_name, "READ");
    for (auto const& pdf_name: pdf_list)
    {
        res.push_back(std::unique_ptr<TH1F>(static_cast<TH1F*>(file->Get(pdf_name))));
    }
    file->Close();
    return res;
}

HistVec2d_t Analyzer::PDFs2dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec2d_t res;
    auto file = std::make_unique<TFile>(file_name, "READ");
    for (auto const& pdf_name: pdf_list)
    {
        res.push_back(std::unique_ptr<TH2F>(static_cast<TH2F*>(file->Get(pdf_name))));
    }
    file->Close();
    return res;
}


void Analyzer::Analyze()
{
    for (auto const& [file_name, channel]: m_file_map)
    {
        TFile* file = TFile::Open(file_name);
        TTree* tree = static_cast<TTree*>(file->Get(m_tree_name));

        TString channel_name = channel == Channel::SL ? "Single lepton" : "Double Lepton";
        m_index = channel == Channel::SL ? GenTruthIdxMapSL : GenTruthIdxMapDL;
        std::cout << "File: " << file_name << ", channel: " << channel_name << "\n";

        HistVec1d_t pdf_1d = PDFs1dFromFile(m_pdf_file_name, pdf1d_names);
        HistVec2d_t pdf_2d = PDFs2dFromFile(m_pdf_file_name, pdf2d_names);

        Event event(tree, channel);
        AnalyzeEvent(event, tree);

        ValidatorInput vi(event, std::move(pdf_1d), std::move(pdf_2d));

        file->Close();
    }
}

void Analyzer::AnalyzeEvent(Event const& event, TTree* tree)
{
    for (long long i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        size_t met_idx = m_index.at("met");
        std::cout << "gen met: " << event.gen_truth.pt[met_idx] << " " << event.gen_truth.mass[met_idx] << "\n";
        std::cout << "reco jet: " << event.reco_jet.nRecoJet << "\n";
        std::cout << "gen jet: " << event.gen_jet.nGenJet << "\n";
        std::cout << "nu: " << event.nu.pdgId[0] << " " << event.nu.pdgId[1] << "\n";
        std::cout << "reco lep: " << event.reco_lep.lep_type[0] << " " << event.reco_lep.lep_iso[0] << " " << event.reco_lep.lep_iso[1] << "\n";
        std::cout << "lep1: " << event.reco_lep.pt[0] << " " << event.reco_lep.eta[0] << " " << event.reco_lep.phi[0] << " " << event.reco_lep.mass[0] << "\n";
        std::cout << "lep2: " << event.reco_lep.pt[1] << " " << event.reco_lep.eta[1] << " " << event.reco_lep.phi[1] << " " << event.reco_lep.mass[1] << "\n";
        std::cout << "\n";
    }
}