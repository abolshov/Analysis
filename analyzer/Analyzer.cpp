#include "Analyzer.hpp"
#include "Definitions.hpp"
#include "EstimatorUtils.hpp"

#include <iostream>

#include "TH1.h"

Analyzer::Analyzer(TString const& tree_name, std::map<TString, Channel> const& input_file_map, TString const& pdf_file_name, Mode mode)
:   m_file_map(input_file_map)
,   m_tree_name(tree_name)  
,   m_estimator(pdf_file_name)
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
    for (ULong64_t evt = 0; evt < 5; ++evt)
    {
        ProcessEvent(evt, tree, ch);
    }
    m_hm.Draw();
    file->Close();
}

void Analyzer::ProcessEvent(ULong64_t evt, TTree* tree, Channel ch)
{
    tree->GetEntry(evt);

    if (m_storage.eventId % 2 != 0)
    {
        return;
    }

    if (m_storage.n_reco_jet < 4)
    {
        return;
    }

    VecLVF_t jets = GetRecoJetP4(m_storage);
    VecLVF_t leptons = GetRecoLepP4(m_storage, ch);
    LorentzVectorF_t met = GetRecoMET(m_storage);
    std::vector<Float_t> jet_resolutions = GetPNetRes(m_storage);

    #ifdef DEBUG 
    VecLVF_t gen_leptons = GetGenLepP4(m_storage, ch);
    VecLVF_t gen_quarks = GetGenQuarksP4(m_storage, ch);
    VecLVF_t gen_nu = GetGenQuarksP4(m_storage, ch);
    LorentzVectorF_t gen_met = GetGenMET(m_storage);

    gen_truth_buf << "Event " << evt << "\n"
                  << "MC truth values:\n"
                  << "\tbq1=(" << gen_quarks[static_cast<size_t>(Quark::b1)].Pt() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b1)].Eta() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b1)].Phi() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b1)].M() << ")\n"
                  << "\tbq2=(" << gen_quarks[static_cast<size_t>(Quark::b2)].Pt() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b2)].Eta() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b2)].Phi() << ", " 
                               << gen_quarks[static_cast<size_t>(Quark::b2)].M() << ")\n";
    if (ch == Channel::SL)
    {
        gen_truth_buf << "\tlq1=(" << gen_quarks[static_cast<size_t>(Quark::q1)].Pt() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q1)].Eta() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q1)].Phi() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q1)].M() << ")\n"
                      << "\tlq2=(" << gen_quarks[static_cast<size_t>(Quark::q2)].Pt() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q2)].Eta() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q2)].Phi() << ", " 
                                   << gen_quarks[static_cast<size_t>(Quark::q2)].M() << ")\n"
                      << "\tlep=(" << gen_leptons[static_cast<size_t>(Lep::lep1)].Pt() << ", " 
                                   << gen_leptons[static_cast<size_t>(Lep::lep1)].Eta() << ", " 
                                   << gen_leptons[static_cast<size_t>(Lep::lep1)].Phi() << ", " 
                                   << gen_leptons[static_cast<size_t>(Lep::lep1)].M() << ")\n"
                      << "\tnu=("  << gen_nu[static_cast<size_t>(Nu::nu1)].Pt() << ", " 
                                   << gen_nu[static_cast<size_t>(Nu::nu1)].Eta() << ", " 
                                   << gen_nu[static_cast<size_t>(Nu::nu1)].Phi() << ", " 
                                   << gen_nu[static_cast<size_t>(Nu::nu1)].M() << ")\n";
    }
    else if (ch == Channel::DL)
    {
        gen_truth_buf << "\tlep1=(" << gen_leptons[static_cast<size_t>(Lep::lep1)].Pt() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep1)].Eta() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep1)].Phi() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep1)].M() << ")\n"
                      << "\tnu1=("  << gen_nu[static_cast<size_t>(Nu::nu1)].Pt() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu1)].Eta() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu1)].Phi() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu1)].M() << ")\n"
                      << "\tlep2=(" << gen_leptons[static_cast<size_t>(Lep::lep2)].Pt() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep2)].Eta() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep2)].Phi() << ", " 
                                    << gen_leptons[static_cast<size_t>(Lep::lep2)].M() << ")\n"
                      << "\tnu2=("  << gen_nu[static_cast<size_t>(Nu::nu2)].Pt() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu2)].Eta() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu2)].Phi() << ", " 
                                    << gen_nu[static_cast<size_t>(Nu::nu2)].M() << ")\n";
    }

    gen_truth_buf << "\tmet=(" << gen_met.Pt() << ", " << gen_met.Eta() << ", " << gen_met.Phi() << ", " << gen_met.M() << ")\n"
                  << "\tmbb=" << (gen_quarks[static_cast<size_t>(Quark::b1)] + gen_quarks[static_cast<size_t>(Quark::b2)]).M() << "\n"
                  << "\tmWhad=" << (gen_quarks[static_cast<size_t>(Quark::q1)] + gen_quarks[static_cast<size_t>(Quark::q2)]).M() << "\n";
    #endif


    auto hme = m_estimator.EstimateMass(jets, leptons, jet_resolutions, met, evt);
    if (hme)
    {
        m_hm.Fill("hme_mass", hme.value());
    }

    #ifdef DEBUG
    gen_truth_buf.str("");
    #endif
}