#include "Validator.hpp"
#include "EstimatorTools.hpp"

#include <iostream>

Validator::Validator()
{
    hm.Add("corr_ratio_px_1", "Ratio of corrected Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_py_1", "Ratio of corrected Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_pz_1", "Ratio of corrected Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_E_1", "Ratio of corrected E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 10}, 100);

    hm.Add("corr_ratio_px_2", "Ratio of corrected Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_py_2", "Ratio of corrected Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_pz_2", "Ratio of corrected Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("corr_ratio_E_2", "Ratio of corrected E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 10}, 100);

    hm.Add("raw_ratio_px_1", "Ratio of Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_py_1", "Ratio of Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_pz_1", "Ratio of Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_E_1", "Ratio of E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 10}, 100);

    hm.Add("raw_ratio_px_2", "Ratio of Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_py_2", "Ratio of Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_pz_2", "Ratio of Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-10, 10}, 100);
    hm.Add("raw_ratio_E_2", "Ratio of E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 10}, 100);
}

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

// void Validator::ResetHistograms()
// {
//     hist_corr_ratio_px_1->Reset("ICESM");
//     hist_corr_ratio_py_1->Reset("ICESM");
//     hist_corr_ratio_pz_1->Reset("ICESM");
//     hist_corr_ratio_E_1->Reset("ICESM");

//     hist_corr_ratio_px_2->Reset("ICESM");
//     hist_corr_ratio_py_2->Reset("ICESM");
//     hist_corr_ratio_pz_2->Reset("ICESM");
//     hist_corr_ratio_E_2->Reset("ICESM");
// }

void Validator::Clear()
{
    gen_truth_p4.clear();
    reco_jet_p4.clear(); 
    nu.clear();
    jet_PNet_resolutions.clear(); 
    jet_PNet_corrections.clear();
}

void Validator::CompareJetsToQuarksB()
{
    // reorder b1 and b2 to make quark with larger pt b1
    if (gen_truth_p4[truth_index["b1"]].Pt() < gen_truth_p4[truth_index["b2"]].Pt())
    {
        std::swap(gen_truth_p4[truth_index["b1"]], gen_truth_p4[truth_index["b2"]]);
    }

    LorentzVectorF_t const& bq1 = gen_truth_p4[truth_index["b1"]];
    LorentzVectorF_t const& bq2 = gen_truth_p4[truth_index["b2"]];

    UHist1d_t& pdf = pdf_1d[pdf1d_index["pdf_b1"]];

    // histograms for accumulating effective values for one event
    auto hist_ratio_px_1 = std::make_unique<TH1F>("hist_ratio_px_1", "hist_ratio_px_1", 100, -10, 10);
    auto hist_ratio_py_1 = std::make_unique<TH1F>("hist_ratio_py_1", "hist_ratio_py_1", 100, -10, 10);
    auto hist_ratio_pz_1 = std::make_unique<TH1F>("hist_ratio_pz_1", "hist_ratio_pz_1", 100, -10, 10);
    auto hist_ratio_E_1 = std::make_unique<TH1F>("hist_ratio_E_1", "hist_ratio_E_1", 100, 0, 10);

    auto hist_ratio_px_2 = std::make_unique<TH1F>("hist_ratio_px_2", "hist_ratio_px_2", 100, -10, 10);
    auto hist_ratio_py_2 = std::make_unique<TH1F>("hist_ratio_py_2", "hist_ratio_py_2", 100, -10, 10);
    auto hist_ratio_pz_2 = std::make_unique<TH1F>("hist_ratio_pz_2", "hist_ratio_pz_2", 100, -10, 10);
    auto hist_ratio_E_2 = std::make_unique<TH1F>("hist_ratio_E_2", "hist_ratio_E_2", 100, 0, 10);

    LorentzVectorF_t raw_bj1, raw_bj2;
    if (reco_jet_p4[0].Pt() > reco_jet_p4[1].Pt())
    {
        raw_bj1 = reco_jet_p4[0];
        raw_bj2 = reco_jet_p4[1];
    }
    else
    {
        raw_bj2 = reco_jet_p4[0];
        raw_bj1 = reco_jet_p4[1];
    }

    hm.Fill("raw_ratio_px_1", raw_bj1.Px()/bq1.Px());
    hm.Fill("raw_ratio_py_1", raw_bj1.Py()/bq1.Py());
    hm.Fill("raw_ratio_pz_1", raw_bj1.Pz()/bq1.Pz());
    hm.Fill("raw_ratio_E_1", raw_bj1.E()/bq1.E());

    hm.Fill("raw_ratio_px_2", raw_bj2.Px()/bq2.Px());
    hm.Fill("raw_ratio_py_2", raw_bj2.Py()/bq2.Py());
    hm.Fill("raw_ratio_pz_2", raw_bj2.Pz()/bq2.Pz());
    hm.Fill("raw_ratio_E_2", raw_bj2.E()/bq2.E());

    for (int i = 0; i < N_ITER; ++i)
    {
        LorentzVectorF_t bj1, bj2;
        if (reco_jet_p4[0].Pt() > reco_jet_p4[1].Pt())
        {
            bj1 = reco_jet_p4[0];
            bj2 = reco_jet_p4[1];
        }
        else
        {
            bj2 = reco_jet_p4[0];
            bj1 = reco_jet_p4[1];
        }

        auto resc_pair = JetRescFact(bj1, bj2, pdf, HIGGS_MASS);
        if (!resc_pair)
        {
            continue;
        }

        auto [c1, c2] = resc_pair.value();
        bj1 *= c1;
        bj2 *= c2;

        double ratio_px_1 = bj1.Px()/bq1.Px();
        double ratio_py_1 = bj1.Py()/bq1.Py();
        double ratio_pz_1 = bj1.Pz()/bq1.Pz();
        double ratio_E_1 = bj1.E()/bq1.E();

        hist_ratio_px_1->Fill(ratio_px_1);
        hist_ratio_py_1->Fill(ratio_py_1);
        hist_ratio_pz_1->Fill(ratio_pz_1);
        hist_ratio_E_1->Fill(ratio_E_1);

        double ratio_px_2 = bj2.Px()/bq2.Px();
        double ratio_py_2 = bj2.Py()/bq2.Py();
        double ratio_pz_2 = bj2.Pz()/bq2.Pz();
        double ratio_E_2 = bj2.E()/bq2.E();

        hist_ratio_px_2->Fill(ratio_px_2);
        hist_ratio_py_2->Fill(ratio_py_2);
        hist_ratio_pz_2->Fill(ratio_pz_2);
        hist_ratio_E_2->Fill(ratio_E_2);
    }

    auto GetEffValue = [](UHist1d_t& h){ return h->GetEntries() ? h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) : 0.0; };
    
    double eff_ratio_px_1 = GetEffValue(hist_ratio_px_1);
    double eff_ratio_py_1 = GetEffValue(hist_ratio_py_1);
    double eff_ratio_pz_1 = GetEffValue(hist_ratio_pz_1);
    double eff_ratio_E_1 = GetEffValue(hist_ratio_E_1);

    double eff_ratio_px_2 = GetEffValue(hist_ratio_px_2);
    double eff_ratio_py_2 = GetEffValue(hist_ratio_py_2);
    double eff_ratio_pz_2 = GetEffValue(hist_ratio_pz_2);
    double eff_ratio_E_2 = GetEffValue(hist_ratio_E_2);

    hm.Fill("corr_ratio_px_1", eff_ratio_px_1);
    hm.Fill("corr_ratio_py_1", eff_ratio_py_1);
    hm.Fill("corr_ratio_pz_1", eff_ratio_pz_1);
    hm.Fill("corr_ratio_E_1", eff_ratio_E_1);

    hm.Fill("corr_ratio_px_2", eff_ratio_px_2);
    hm.Fill("corr_ratio_py_2", eff_ratio_py_2);
    hm.Fill("corr_ratio_pz_2", eff_ratio_pz_2);
    hm.Fill("corr_ratio_E_2", eff_ratio_E_2);

    Clear();
}