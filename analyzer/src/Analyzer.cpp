#include "Analyzer.hpp"
#include "Definitions.hpp"
#include "EstimatorUtils.hpp"

#ifdef DEV
    #include "MatchingTools.hpp"
    #include "SelectionUtils.hpp"
#endif 

#include <iostream>
#include <chrono>

#include "TH1.h"

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, std::map<Channel, TString> const& pdf_file_map)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)  
,   m_estimator(pdf_file_map.at(Channel::SL), pdf_file_map.at(Channel::DL), AggregationMode::Combination)
{
    TH1::AddDirectory(false);
    m_hm.Add("hme_mass", "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2500}, 100);
    m_hm.Add("hme_width", "HME width", {"width, [GeV]", "Count"}, {0, 1000}, 100);
    m_hm.Add("hme_peak", "HME peak value", {"peak value", "Count"}, {0, 200}, 100);
    m_hm.Add("hme_integral", "HME integral", {"integral", "Count"}, {0, 10000}, 100);
}

void Analyzer::ProcessFile(TString const& name, Channel ch)
{
    std::cout << "Start processing " << name << "\n";
    auto process_file_start{std::chrono::steady_clock::now()};
    TFile* file = TFile::Open(name);
    TTree* tree = static_cast<TTree*>(file->Get<TTree>(m_tree_name));

    #ifdef DEV
        if (m_record_iterations)
        {
            TString dbg_file_name = ch == Channel::SL ? "hme_iter_sl.root" : "hme_iter_dl.root";
            m_estimator.OpenDbgFile(dbg_file_name, ch);
        }
    #endif

    m_storage.ConnectTree(tree, ch);
    ULong64_t n_events = tree->GetEntries();
    std::chrono::duration<double> average_duration{};
    for (ULong64_t evt = 0; evt < n_events; ++evt)
    {
        auto process_event_start{std::chrono::steady_clock::now()};
        ProcessEvent(evt, tree, ch);
        auto process_event_end{std::chrono::steady_clock::now()};
        std::chrono::duration<double> process_event_duration{process_event_end - process_event_start};
        average_duration += process_event_duration;
        if (evt % 5000 == 0 && evt > 0)
        {
            std::cout << "===> Processed " << 100.0*evt/n_events << "%, time elapsed " 
                      << std::chrono::duration<double>(average_duration).count() << " s, average time per event: " 
                      << std::chrono::duration<double, std::milli>(average_duration).count()/counter << " ms\n";
        }
    }
    std::cout << "counter=" << counter << "\n";
    m_hm.Draw();
    file->Close();
    auto process_file_end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds{process_file_end - process_file_start};
    std::cout << "Processing " << name << " took " << elapsed_seconds.count() << " seconds\n";
    std::cout << "Average per event: " << std::chrono::duration<double, std::milli>(average_duration).count()/counter << " ms\n";
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
    
    #ifdef DEV 
        if (!IsRecoverable(m_storage, ch))
        {
            return;
        }

        if (!IsFiducial(m_storage, jets, ch))
        {
            return;
        }
    #endif

    ++counter;

    auto hme = m_estimator.EstimateMass(jets, leptons, met, evt_id, ch);
    if (hme)
    {
        auto out_array = hme.value();
        m_hm.Fill("hme_mass", out_array[static_cast<size_t>(EstimOut::mass)]);
        m_hm.Fill("hme_width", out_array[static_cast<size_t>(EstimOut::width)]);
        m_hm.Fill("hme_peak", out_array[static_cast<size_t>(EstimOut::peak_value)]);
        m_hm.Fill("hme_integral", out_array[static_cast<size_t>(EstimOut::integral)]);
    }
}