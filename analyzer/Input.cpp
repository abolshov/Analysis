#include "Input.hpp"

EstimatorInput::EstimatorInput(std::vector<TLorentzVector>&& p4, std::string const& pdf_file)
:   p4(std::move(p4))
{
    auto file_pdf = std::make_unique<TFile>(pdf_file.c_str(), "READ");

    size_t n_pdf1d = static_cast<size_t>(PDF1dSLRes::count);
    pdf1d.reserve(n_pdf1d);
    for (size_t i = 0; i < n_pdf1d; ++i)
    {
        auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get(pdf1d_names.at(i))));
        pdf1d.push_back(std::move(pdf));
    }

    size_t n_pdf2d = static_cast<size_t>(PDF2dSLRes::count);
    pdf2d.reserve(n_pdf2d);
    for (size_t i = 0; i < n_pdf2d; ++i)
    {
        auto pdf = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get(pdf2d_names.at(i))));
        pdf2d.push_back(std::move(pdf));
    }

    file_pdf->Close();
}

EstimatorInput::EstimatorInput(std::vector<TLorentzVector>&& p4, std::vector<UHist1D>&& vec_pdf_1d, std::vector<UHist2D>&& vec_pdf_2d)
:   p4(std::move(p4))
,   pdf1d(std::move(vec_pdf_1d))
,   pdf2d(std::move(vec_pdf_2d))
{}

ValidatorInput::ValidatorInput(Event const& event, Vec1D_t&& pdf_1d, Vec2D_t&& pdf_2d)
:   pdf1d(std::move(pdf_1d))
,   pdf2d(std::move(pdf_2d))
{
    for (int i = 0; i < event.reco_jet.nRecoJet; ++i)
    {
        TLorentzVector p;
        p.SetPtEtaPhiM(event.reco_jet.pt[i], event.reco_jet.eta[i], event.reco_jet.phi[i], event.reco_jet.mass[i]);
        reco_jet_p4.push_back(p);

        jet_PNet_resolutions.push_back(event.reco_jet.PNetRegPtRawCorr[i]);
        jet_PNet_corrections.push_back(event.reco_jet.pt[i]*event.reco_jet.PNetRegPtRawRes[i]);
    }

    // gen truth inputs
    for (size_t i = 0; i < event.m_index.size(); ++i)
    {
        TLorentzVector p;
        p.SetPtEtaPhiM(event.gen_truth.pt[i], event.gen_truth.eta[i], event.gen_truth.phi[i], event.gen_truth.mass[i]);
        gen_truth_p4.push_back(p);
    }

    // true nu inputs
    for (size_t i = 0; i < event.m_nu_index.size(); ++i)
    {
        TLorentzVector p;
        p.SetPtEtaPhiM(event.nu.pt[i], event.nu.eta[i], event.nu.phi[i], event.nu.mass[i]);
        nu.push_back(p);
    }

    recoMET.SetPtEtaPhiM(event.reco_met_pt, 0.0, event.reco_met_phi, 0.0);
}