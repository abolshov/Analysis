#include "Analyzer.hpp"
#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "MatchingTools.hpp"
#include "SelectionUtils.hpp"

#include <iostream>

#include "TH1.h"

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, std::map<Channel, TString> const& pdf_file_map)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)  
,   m_estimator(pdf_file_map.at(Channel::SL), pdf_file_map.at(Channel::DL), "")
{
    TH1::AddDirectory(false);
    m_hm.Add("hme_mass", "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2500}, 100);
    m_hm.Add("hme_score", "HME score", {"score", "Count"}, {0, 4}, 50);
    m_hm.Add("hme_width", "HME width", {"width, [GeV]", "Count"}, {0, 1000}, 100);
    m_hm.Add("hme_peak", "HME peak value", {"peak value", "Count"}, {0, 200}, 100);
    m_hm.Add("hme_integral", "HME integral", {"integral", "Count"}, {0, 10000}, 100);
}

void Analyzer::ProcessFile(TString const& name, Channel ch)
{
    TFile* file = TFile::Open(name);
    TTree* tree = static_cast<TTree*>(file->Get<TTree>(m_tree_name));

    m_storage.ConnectTree(tree, ch);
    ULong64_t n_events = tree->GetEntries();
    for (ULong64_t evt = 0; evt < n_events; ++evt)
    {
        ProcessEvent(evt, tree, ch);
    }
    std::cout << "counter=" << counter << "\n";
    m_hm.Draw();
    file->Close();
}

void Analyzer::ProcessEvent(ULong64_t evt, TTree* tree, Channel ch)
{
    tree->GetEntry(evt);

    ULong64_t evt_id = m_storage.event_id;
    if (evt_id % 2 != 1)
    {
        return;
    }

    VecLVF_t jets = GetRecoJetP4(m_storage);
    VecLVF_t leptons = GetRecoLepP4(m_storage, ch);
    LorentzVectorF_t met = GetRecoMET(m_storage);

    if (!IsRecoverable(m_storage, ch))
    {
        return;
    }

    if (!IsFiducial(m_storage, jets, ch))
    {
        return;
    }

    ++counter;

    auto hme = m_estimator.EstimateMass(jets, leptons, met, evt_id, ch);
    if (hme)
    {
        auto out_array = hme.value();
        m_hm.Fill("hme_mass", out_array[static_cast<size_t>(EstimOut::mass)]);
        m_hm.Fill("hme_score", out_array[static_cast<size_t>(EstimOut::score)]);
        m_hm.Fill("hme_width", out_array[static_cast<size_t>(EstimOut::width)]);
        m_hm.Fill("hme_peak", out_array[static_cast<size_t>(EstimOut::peak_value)]);
        m_hm.Fill("hme_integral", out_array[static_cast<size_t>(EstimOut::integral)]);
    }
}