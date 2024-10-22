#include "Input.hpp"

EstimatorInput::EstimatorInput(Event const& event, HistVec1d_t&& vec_pdf_1d, HistVec2d_t&& vec_pdf_2d)
:   pdf1d(std::move(vec_pdf_1d))
,   pdf2d(std::move(vec_pdf_2d))
{
    for (size_t i = 0; i < event.m_index.size(); ++i)
    {
        p4.emplace_back(event.gen_truth.pt[i], event.gen_truth.eta[i], event.gen_truth.phi[i], event.gen_truth.mass[i]);
    }
}

ValidatorInput::ValidatorInput(Event const& event, HistVec1d_t&& pdf_1d, HistVec2d_t&& pdf_2d)
:   pdf1d(std::move(pdf_1d))
,   pdf2d(std::move(pdf_2d))
,   recoMET(event.reco_met_pt, 0.0, event.reco_met_phi, 0.0)
{
    for (int i = 0; i < event.reco_jet.nRecoJet; ++i)
    {
        reco_jet_p4.emplace_back(event.reco_jet.pt[i], event.reco_jet.eta[i], event.reco_jet.phi[i], event.reco_jet.mass[i]);

        jet_PNet_resolutions.push_back(event.reco_jet.PNetRegPtRawCorr[i]);
        jet_PNet_corrections.push_back(event.reco_jet.pt[i]*event.reco_jet.PNetRegPtRawRes[i]);
    }

    // gen truth inputs
    for (size_t i = 0; i < event.m_index.size(); ++i)
    {
        gen_truth_p4.emplace_back(event.gen_truth.pt[i], event.gen_truth.eta[i], event.gen_truth.phi[i], event.gen_truth.mass[i]);
    }

    // true nu inputs
    for (size_t i = 0; i < event.m_nu_index.size(); ++i)
    {
        nu.emplace_back(event.nu.pt[i], event.nu.eta[i], event.nu.phi[i], event.nu.mass[i]);
    }
}