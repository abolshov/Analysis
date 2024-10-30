#include "Validator.hpp"
#include "EstimatorTools.hpp"

#include <iostream>
#include <stdexcept>

Validator::Validator()
{
    // momentum ratios of quarks and jets
    hm.Add("corr_ratio_px_1", "Ratio of corrected Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_py_1", "Ratio of corrected Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_pz_1", "Ratio of corrected Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_E_1", "Ratio of corrected E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("corr_ratio_px_2", "Ratio of corrected Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_py_2", "Ratio of corrected Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_pz_2", "Ratio of corrected Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("corr_ratio_E_2", "Ratio of corrected E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("raw_ratio_px_1", "Ratio of Px of jet 1 to Px of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_py_1", "Ratio of Py of jet 1 to Py of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_pz_1", "Ratio of Pz of jet 1 to Pz of quark 1", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_E_1", "Ratio of E of jet 1 to E of quark 1", {"Ratio", "Count"}, {0, 25}, 100);

    hm.Add("raw_ratio_px_2", "Ratio of Px of jet 2 to Px of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_py_2", "Ratio of Py of jet 2 to Py of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_pz_2", "Ratio of Pz of jet 2 to Pz of quark 2", {"Ratio", "Count"}, {-25, 25}, 100);
    hm.Add("raw_ratio_E_2", "Ratio of E of jet 2 to E of quark 2", {"Ratio", "Count"}, {0, 25}, 100);

    // momenta of quarks and jets
    hm.Add("corr_px_1", "Corrected Px of jet 1", {"Px, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_py_1", "Corrected Py of jet 1", {"Py, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_pz_1", "Corrected Pz of jet 1", {"Pz, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_E_1", "Corrected E of jet 1", {"E, [GeV]", "Count"}, {0.0, 1000.0}, 100);

    hm.Add("corr_px_2", "Corrected Px of jet 2", {"Px, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_py_2", "Corrected Py of jet 2", {"Py, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_pz_2", "Corrected Pz of jet 2", {"Pz, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("corr_E_2", "Corrected E of jet 2", {"E, [GeV]", "Count"}, {0.0, 1000.0}, 100);

    hm.Add("true_px_1", "True Px of quark 1", {"Px, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_py_1", "True Py of quark 1", {"Py, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_pz_1", "True Pz of quark 1", {"Pz, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_E_1", "True E of quark 1", {"E, [GeV]", "Count"}, {0, 1000.0}, 100);

    hm.Add("true_px_2", "True Px of quark 2", {"Px, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_py_2", "True Py of quark 2", {"Py, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_pz_2", "True Pz of quark 2", {"Pz, [GeV]", "Count"}, {-1000.0, 1000.0}, 200);
    hm.Add("true_E_2", "True E of quark 2", {"E, [GeV]", "Count"}, {0, 1000.0}, 100);
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

void Validator::Clear()
{
    gen_truth_p4.clear();
    reco_jet_p4.clear(); 
    nu.clear();
    jet_PNet_resolutions.clear(); 
    jet_PNet_corrections.clear();
}

void Validator::CompareJetsToQuarksB(bool use_2d_pdf, bool skip_failures)
{
    // reorder b1 and b2 to make quark with larger pt b1
    if (gen_truth_p4[truth_index["b1"]].Pt() < gen_truth_p4[truth_index["b2"]].Pt())
    {
        std::swap(gen_truth_p4[truth_index["b1"]], gen_truth_p4[truth_index["b2"]]);
    }

    LorentzVectorF_t const& bq1 = gen_truth_p4[truth_index["b1"]];
    LorentzVectorF_t const& bq2 = gen_truth_p4[truth_index["b2"]];

    hm.Fill("true_px_1", bq1.Px());
    hm.Fill("true_py_1", bq1.Px());
    hm.Fill("true_pz_1", bq1.Py());
    hm.Fill("true_E_1", bq1.E());

    hm.Fill("true_px_2", bq2.Px());
    hm.Fill("true_py_2", bq2.Px());
    hm.Fill("true_pz_2", bq2.Py());
    hm.Fill("true_E_2", bq2.E());

    UHist1d_t& pdf1d = pdf_1d[pdf1d_index["pdf_b1"]];
    UHist2d_t& pdf2d = pdf_2d[pdf2d_index["pdf_b1b2"]];

    // histograms for accumulating effective ratio values for one event
    auto hist_ratio_px_1 = std::make_unique<TH1F>("hist_ratio_px_1", "hist_ratio_px_1", 100, -25, 25);
    auto hist_ratio_py_1 = std::make_unique<TH1F>("hist_ratio_py_1", "hist_ratio_py_1", 100, -25, 25);
    auto hist_ratio_pz_1 = std::make_unique<TH1F>("hist_ratio_pz_1", "hist_ratio_pz_1", 100, -25, 25);
    auto hist_ratio_E_1 = std::make_unique<TH1F>("hist_ratio_E_1", "hist_ratio_E_1", 100, 0, 25);

    auto hist_ratio_px_2 = std::make_unique<TH1F>("hist_ratio_px_2", "hist_ratio_px_2", 100, -25, 25);
    auto hist_ratio_py_2 = std::make_unique<TH1F>("hist_ratio_py_2", "hist_ratio_py_2", 100, -25, 25);
    auto hist_ratio_pz_2 = std::make_unique<TH1F>("hist_ratio_pz_2", "hist_ratio_pz_2", 100, -25, 25);
    auto hist_ratio_E_2 = std::make_unique<TH1F>("hist_ratio_E_2", "hist_ratio_E_2", 100, 0, 25);

    // histograms for accumulating effective momentum values for one event
    auto hist_px_1 = std::make_unique<TH1F>("hist_px_1", "hist_px_1", 200, -1000.0, 1000.0);
    auto hist_py_1 = std::make_unique<TH1F>("hist_py_1", "hist_py_1", 200, -1000.0, 1000.0);
    auto hist_pz_1 = std::make_unique<TH1F>("hist_pz_1", "hist_pz_1", 200, -1000.0, 1000.0);
    auto hist_E_1 = std::make_unique<TH1F>("hist_E_1", "hist_E_1", 100, 0.0, 1000.0);

    auto hist_px_2 = std::make_unique<TH1F>("hist_px_2", "hist_px_2", 200, -1000.0, 1000.0);
    auto hist_py_2 = std::make_unique<TH1F>("hist_py_2", "hist_py_2", 200, -1000.0, 1000.0);
    auto hist_pz_2 = std::make_unique<TH1F>("hist_pz_2", "hist_pz_2", 200, -1000.0, 1000.0);
    auto hist_E_2 = std::make_unique<TH1F>("hist_E_2", "hist_E_2", 100, 0.0, 1000.0);

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

        double c1 = 1.0;
        double c2 = 1.0;

        if (use_2d_pdf)
        {
            pdf2d->GetRandom2(c1, c2);
        }
        else 
        {
            auto resc_pair = JetRescFact(bj1, bj2, pdf1d, HIGGS_MASS);
            if (resc_pair)
            {
                auto [resc1, resc2] = resc_pair.value();
                c1 = resc1;
                c2 = resc2;
            }
            else 
            {
                if (skip_failures)
                {
                    continue;
                }    
            }
        }

        bj1 *= c1;
        bj2 *= c2;

        hist_px_1->Fill(bj1.Px());
        hist_py_1->Fill(bj1.Py());
        hist_pz_1->Fill(bj1.Pz());
        hist_E_1->Fill(bj1.E());

        hist_px_2->Fill(bj2.Px());
        hist_py_2->Fill(bj2.Py());
        hist_pz_2->Fill(bj2.Pz());
        hist_E_2->Fill(bj2.E());

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

    double eff_px_1 = GetEffValue(hist_px_1);
    double eff_py_1 = GetEffValue(hist_py_1);
    double eff_pz_1 = GetEffValue(hist_pz_1);
    double eff_E_1 = GetEffValue(hist_E_1);

    double eff_px_2 = GetEffValue(hist_px_2);
    double eff_py_2 = GetEffValue(hist_py_2);
    double eff_pz_2 = GetEffValue(hist_pz_2);
    double eff_E_2 = GetEffValue(hist_E_2);

    hm.Fill("corr_px_1", eff_px_1);
    hm.Fill("corr_py_1", eff_py_1);
    hm.Fill("corr_pz_1", eff_pz_1);
    hm.Fill("corr_E_1", eff_E_1);

    hm.Fill("corr_px_2", eff_px_2);
    hm.Fill("corr_py_2", eff_py_2);
    hm.Fill("corr_pz_2", eff_pz_2);
    hm.Fill("corr_E_2", eff_E_2);

    Clear();
}