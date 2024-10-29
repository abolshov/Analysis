#include "Validator.hpp"

void Validator::FillVariables(Event const& event)
{
    recoMET = LorentzVectorF_t(event.reco_met_pt, 0.0, event.reco_met_phi, 0.0);

    for (int i = 0; i < event.reco_jet.nRecoJet; ++i)
    {
        reco_jet_p4.emplace_back(event.reco_jet.pt[i], event.reco_jet.eta[i], event.reco_jet.phi[i], event.reco_jet.mass[i]);

        jet_PNet_resolutions.push_back(event.reco_jet.pt[i]*event.reco_jet.PNetRegPtRawRes[i]);
        jet_PNet_corrections.push_back(event.reco_jet.PNetRegPtRawCorr[i]);
    }

    // gen truth inputs
    for (auto const& [obj_name, idx]: event.m_index)
    {
        gen_truth_p4.emplace_back(event.gen_truth.pt[idx], event.gen_truth.eta[idx], event.gen_truth.phi[idx], event.gen_truth.mass[idx]);
    }

    // true nu inputs
    for (auto const& [obj_name, idx]: event.m_nu_index)
    {
        nu.emplace_back(event.nu.pt[idx], event.nu.eta[idx], event.nu.phi[idx], event.nu.mass[idx]);
    }
}

void Validator::Print()
{
    for (int i = 0; i < 10; ++i)
    {
        std::cout << pdf_1d[0]->GetRandom() << "\n";
    }
}