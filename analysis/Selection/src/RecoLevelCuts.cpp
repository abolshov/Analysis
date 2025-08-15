#include "RecoLevelCuts.hpp"
#include "Utils.hpp"

AK4JetAcceptCut::AK4JetAcceptCut(Float_t pt, Float_t eta, Int_t num_jets)
:   Specification<Event>("ak4_accept_cut", StrCat("pt > ", std::to_string(pt), " && abs(eta) < ", std::to_string(eta), " && n >= ", std::to_string(num_jets)))
,   m_pt(pt)
,   m_eta(eta)
,   m_num_jets(num_jets)
{}

bool AK4JetAcceptCut::IsSatisfied(Event const& event)
{
    Int_t nj = event.n_reco_jet;
    if (nj < m_num_jets)
    {
        return false;
    }

    Int_t num_jets_in_accept = 0;
    for (Int_t i = 0; i < nj; ++i)
    {
        if (event.reco_jet_pt[i] > m_pt && std::abs(event.reco_jet_eta[i]) < m_eta)
        {
            ++num_jets_in_accept;
        }
    }

    return num_jets_in_accept >= m_num_jets;
}

AK8JetAcceptCut::AK8JetAcceptCut(Float_t pt, Float_t eta, Int_t num_jets)
:   Specification<Event>("ak8_accept_cut", StrCat("pt > ", std::to_string(pt), " && abs(eta) < ", std::to_string(eta), " && n >= ", std::to_string(num_jets)))
,   m_pt(pt)
,   m_eta(eta)
,   m_num_jets(num_jets)
{}

bool AK8JetAcceptCut::IsSatisfied(Event const& event)
{
    Int_t nfj = event.n_reco_fatjet;
    if (nfj < m_num_jets)
    {
        return false;
    }

    Int_t num_fatjets_in_accept = 0;
    for (Int_t i = 0; i < nfj; ++i)
    {
        if (event.reco_jet_pt[i] > m_pt && std::abs(event.reco_jet_eta[i]) < m_eta)
        {
            ++num_fatjets_in_accept;
        }
    }

    return num_fatjets_in_accept >= m_num_jets;
}