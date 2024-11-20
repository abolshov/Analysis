#include "Analyzer.hpp"
#include "Definitions.hpp"

#include <iostream>

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)     
{
    if (mode == Mode::Validation)
    {
        std::cout << "Constructing analyzer in validation mode\n";
    }
}

HistVec1d_t Analyzer::PDFs1dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec1d_t res;
    auto file = std::unique_ptr<TFile>(TFile::Open(file_name, "READ"));
    for (auto const& pdf_name: pdf_list)
    {
        std::unique_ptr<TH1F> hist(file->Get<TH1F>(pdf_name));
        hist->SetDirectory(nullptr); // unregister hist from root's garbage collector to avoid double free
        res.push_back(std::move(hist));
    }
    file->Close();
    return res;
}

HistVec2d_t Analyzer::PDFs2dFromFile(TString const& file_name, std::vector<TString> const& pdf_list) const
{
    HistVec2d_t res;
    auto file = std::unique_ptr<TFile>(TFile::Open(file_name, "READ"));
    for (auto const& pdf_name: pdf_list)
    {
        std::unique_ptr<TH2F> hist(file->Get<TH2F>(pdf_name));
        hist->SetDirectory(nullptr);
        res.push_back(std::move(hist));
    }
    file->Close();
    return res;
}

void Analyzer::ProcessFile(TString const& name, Channel ch)
{
    TFile* file = TFile::Open(name);
    TTree* tree = static_cast<TTree*>(file->Get(m_tree_name));

    m_storage.ConnectTree(tree, ch);
    
    ULong64_t n_events = tree->GetEntries();
    for (ULong64_t i = 0; i < n_events; ++i)
    {
        tree->GetEntry(i);
        auto res = GetPNetRes();
        for (auto r: res)
        {
            std::cout << r << " ";
        } 
        std::cout << "\n";
    }
}

VecLVF_t Analyzer::GetRecoJetP4()
{
    VecLVF_t res;
    for (int i = 0; i < m_storage.n_reco_jet; ++i)
    {
        res.emplace_back(m_storage.reco_jet_pt[i], m_storage.reco_jet_eta[i], m_storage.reco_jet_phi[i], m_storage.reco_jet_mass[i]);
    }
    return res;
}

VecLVF_t Analyzer::GetRecoLepP4(Channel ch)
{
    VecLVF_t res;
    res.emplace_back(m_storage.reco_lep_pt[0], m_storage.reco_lep_eta[0], m_storage.reco_lep_phi[0], m_storage.reco_lep_mass[0]);
    if (ch == Channel::DL)
    {
        res.emplace_back(m_storage.reco_lep_pt[1], m_storage.reco_lep_eta[1], m_storage.reco_lep_phi[1], m_storage.reco_lep_mass[1]);
    }
    return res;
}

std::vector<Float_t> Analyzer::GetPNetRes()
{
    std::vector<Float_t> res;
    for (int i = 0; i < m_storage.n_reco_jet; ++i)
    {
        res.push_back(m_storage.reco_jet_pt[i]*m_storage.reco_jet_corr[i]*m_storage.reco_jet_res[i]);
    }
    return res;
}

// VecLVF_t Analyzer::GetGenJetP4()
// VecLVF_t Analyzer::GetGenLepP4()