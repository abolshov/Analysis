#include "Analyzer.hpp"
#include "Definitions.hpp"
#include "EstimatorUtils.hpp"
#include "MatchingTools.hpp"
#include "SelectionUtils.hpp"

#include <iostream>

#include "TH1.h"

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)  
,   m_estimator(pdf_file_name, "")
,   m_hm()   
{
    TH1::AddDirectory(false);
    m_hm.Add("hme_mass", "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2500}, 100);
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

    if (m_storage.eventId % 2 != 1)
    {
        return;
    }

    VecLVF_t jets = GetRecoJetP4(m_storage);
    VecLVF_t leptons = GetRecoLepP4(m_storage, ch);
    LorentzVectorF_t met = GetRecoMET(m_storage);
    std::vector<Float_t> jet_resolutions = GetPNetRes(m_storage);

    if (!IsRecoverable(m_storage, ch))
    {
        return;
    }

    if (!IsFiducial(m_storage, jets, ch))
    {
        return;
    }

    ++counter;

    #ifdef DEBUG 
        VecLVF_t gen_leptons = GetGenLepP4(m_storage, ch);
        VecLVF_t gen_quarks = GetGenQuarksP4(m_storage, ch);
        VecLVF_t gen_nu = GetGenQuarksP4(m_storage, ch);
        LorentzVectorF_t gen_met = GetGenMET(m_storage);

        gen_truth_buf << "Event " << evt << ", eventId=" << m_storage.eventId << "\nMC truth values:\n";
        LogP4(gen_truth_buf, gen_quarks[static_cast<size_t>(Quark::b1)], "bq1");
        LogP4(gen_truth_buf, gen_quarks[static_cast<size_t>(Quark::b2)], "bq2");

        if (ch == Channel::SL)
        {
            LogP4(gen_truth_buf, gen_quarks[static_cast<size_t>(Quark::q1)], "lq1");
            LogP4(gen_truth_buf, gen_quarks[static_cast<size_t>(Quark::q2)], "lq2");
            LogP4(gen_truth_buf, gen_leptons[static_cast<size_t>(Lep::lep1)], "lep");
            LogP4(gen_truth_buf, gen_nu[static_cast<size_t>(Nu::nu1)], "nu");
        }
        else if (ch == Channel::DL)
        {
            LogP4(gen_truth_buf, gen_leptons[static_cast<size_t>(Lep::lep1)], "lep1");
            LogP4(gen_truth_buf, gen_nu[static_cast<size_t>(Nu::nu1)], "nu1");
            LogP4(gen_truth_buf, gen_leptons[static_cast<size_t>(Lep::lep2)], "lep2");
            LogP4(gen_truth_buf, gen_nu[static_cast<size_t>(Nu::nu2)], "nu2");
        }
        else
        {
            assert(false);
        }

        LogP4(gen_truth_buf, gen_met, "met");

        if (ch == Channel::SL)
        {
            gen_truth_buf << "\tmbb=" << (gen_quarks[static_cast<size_t>(Quark::b1)] + gen_quarks[static_cast<size_t>(Quark::b2)]).M() << "\n"
                          << "\tmWhad=" << (gen_quarks[static_cast<size_t>(Quark::q1)] + gen_quarks[static_cast<size_t>(Quark::q2)]).M() << "\n";
        }
    #endif

    TString chosen_comb = "";
    // auto hme = m_estimator.EstimateMass(jets, leptons, jet_resolutions, met, evt, chosen_comb); // sl
    auto hme = m_estimator.EstimateMass(jets, leptons, met, evt, chosen_comb);
    // auto hme = m_estimator.EstimateMass(jets, leptons, met, evt, chosen_comb); // dl
    if (hme)
    {
        m_hm.Fill("hme_mass", hme.value());
    }

    #ifdef DEBUG
        gen_truth_buf.str("");
    #endif
}